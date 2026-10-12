#!/usr/bin/env python3
"""Run equivalent dense/CSR matvec implementations, sequentially, on Linux.

Timing and process-tree PSS are separate passes. No performance assertion is
used as a CI gate. See README.md for dependencies and measurement boundaries.
"""

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import signal
import statistics
import subprocess
import sys
import time

from reference import reference

HERE = Path(__file__).resolve().parent
IMPLEMENTATIONS = [
    ("cpp", "serial", False),
    ("cpp", "std-thread", True),
    ("cpp", "omp-static", True),
    ("cpp", "omp-dynamic", True),
    ("cpp", "tbb", True),
    ("rust", "serial", False),
    ("rust", "std-scoped", True),
    ("rust", "rayon", True),
    ("python", "serial", False),
    ("python", "thread", True),
    ("python", "process-copy", True),
    ("python", "process-shared", True),
    ("python", "numpy-serial", False),
    ("python", "numpy-thread", True),
    ("python", "numba-serial", False),
    ("python", "numba-parallel", True),
]


def read_text(path):
    try:
        return Path(path).read_text().strip()
    except (OSError, UnicodeError):
        return None


def tree_pids(pid):
    pending, seen = [pid], set()
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        # A subprocess may be created by an executor's management thread.
        for task in Path(f"/proc/{current}/task").glob('*/children'):
            children = read_text(task)
            if children:
                pending.extend(map(int, children.split()))
    return seen


def tree_pss(pid):
    """Contemporaneous live-process sum, not a sum of individual maxima."""
    total, count = 0, 0
    for current in tree_pids(pid):
        rollup = read_text(f"/proc/{current}/smaps_rollup")
        if rollup:
            values = dict(line.split(":", 1) for line in rollup.splitlines() if ":" in line)
            if "Pss" in values:
                total += int(values["Pss"].split()[0]) * 1024
                count += 1
    return total, count


def memory_budget_ok(pid, limit_mib):
    # RSS is deliberately conservative here: shared pages can be counted twice.
    # It is a stop condition, never the reported space metric.
    total = 0
    for current in tree_pids(pid):
        status = read_text(f'/proc/{current}/status') or ''
        total += next((int(line.split()[1]) * 1024 for line in status.splitlines() if line.startswith('VmRSS:')), 0)
    return total <= limit_mib * 2**20


def available_memory_mib():
    return next(int(line.split()[1]) / 1024 for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:'))


def quantiles(values):
    ordered = sorted(values)

    def percentile(fraction):
        index = (len(ordered) - 1) * fraction
        lower = int(index)
        upper = min(lower + 1, len(ordered) - 1)
        return ordered[lower] + (ordered[upper] - ordered[lower]) * (index - lower)

    return {"median": statistics.median(values), "q1": percentile(0.25),
            "q3": percentile(0.75), "min": min(values), "max": max(values),
            "p95": percentile(0.95) if len(values) >= 20 else None,
            "p99": percentile(0.99) if len(values) >= 100 else None, "count": len(values)}


def execute(command, env, memory, log_prefix, timeout, interval_ms, memory_limit_mib):
    peak_pss, max_processes, samples = 0, 0, 0
    start = time.perf_counter()
    last_budget_check = start
    with log_prefix.with_suffix(".out").open("w") as stdout, log_prefix.with_suffix(".err").open("w") as stderr:
        proc = subprocess.Popen(command, env=env, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            while proc.poll() is None:
                if time.perf_counter() - start > timeout:
                    raise TimeoutError(f"Timed out: {command[0]}")
                if time.perf_counter() - last_budget_check >= 1:
                    if not memory_budget_ok(proc.pid, memory_limit_mib) or available_memory_mib() < 2048:
                        raise MemoryError('Benchmark stopped at memory budget; no incomplete result is published')
                    last_budget_check = time.perf_counter()
                if memory:
                    pss, count = tree_pss(proc.pid)
                    peak_pss = max(peak_pss, pss)
                    max_processes = max(max_processes, count)
                    samples += 1
                time.sleep(interval_ms / 1000 if memory else 0.001)
            wall_ms = (time.perf_counter() - start) * 1000
        except BaseException:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            raise
    error = log_prefix.with_suffix(".err").read_text()
    output = log_prefix.with_suffix(".out").read_text()
    if proc.returncode:
        raise RuntimeError(f"Exit {proc.returncode}: {error}\n{output}")
    payload = json.loads(output.strip().splitlines()[-1])
    if memory and not peak_pss:
        raise RuntimeError("PSS unavailable: do not substitute parent-only RSS")
    return {"measurement": payload, "process_wall_ms": wall_ms,
            "sampled_peak_tree_pss_bytes": peak_pss if memory else None,
            "max_observed_processes": max_processes if memory else None,
            "pss_samples": samples if memory else None, "stderr": error}


def environment(args):
    topology = subprocess.check_output(["lscpu", "-J"], text=True)
    return {
        "measured_at": dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat(),
        "system": platform.system(), "kernel": platform.release(), "machine": platform.machine(),
        "cpu": json.loads(topology),
        "cpu_topology": subprocess.check_output(['lscpu', '-p=CPU,CORE,SOCKET,MAXMHZ'], text=True),
        "affinity": args.cpus, "monitor_affinity": args.monitor_cpus,
        "memory_total_kib": next((int(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemTotal:')), None),
        "loadavg_start": list(os.getloadavg()),
        "cgroup_cpu_max": read_text('/sys/fs/cgroup/cpu.max'),
        "cgroup_cpuset": read_text('/sys/fs/cgroup/cpuset.cpus.effective'),
        "cgroup_membership": read_text('/proc/self/cgroup'),
        "governor": read_text('/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor'),
        "seed": args.seed, "fresh_processes": args.launches, "warm_samples_per_process": args.repeats,
        "profile": args.profile, "stress_seconds": args.seconds,
        "memory_guard_tree_rss_mib": args.memory_limit_mib, "memory_guard_interval_ms": 1000,
        "memory_passes_per_case": 1, "pss_interval_ms": args.pss_interval_ms,
        "memory_final_hold_ms": 100, "timing_memory_separate": True,
        "cold_cache": False,
        "thread_env": {"OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "NUMBA_THREADING_LAYER": "omp"},
        "sources_sha256": {str(path.relative_to(HERE)): hashlib.sha256(path.read_bytes()).hexdigest()
                           for path in sorted(HERE.rglob('*')) if path.is_file()
                           and path.suffix in {'.py', '.cpp', '.rs', '.toml', '.lock', '.sh', '.txt'}
                           and '__pycache__' not in str(path)},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpp', required=True, type=Path)
    parser.add_argument('--rust', required=True, type=Path)
    parser.add_argument('--python', default=sys.executable, type=Path)
    parser.add_argument('--python-script', default=HERE / 'python' / 'bench.py', type=Path)
    parser.add_argument('--cpus', default='0,2,4,6')
    parser.add_argument('--monitor-cpus', default='8,9,10,11,12,13,14,15')
    parser.add_argument('--workers', default='1,2,4')
    parser.add_argument('--profile', choices=['baseline', 'stress'], default='baseline')
    parser.add_argument('--seconds', type=float, default=0)
    parser.add_argument('--launches', type=int, default=3)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--seed', type=int, default=20261012)
    parser.add_argument('--pss-interval-ms', type=float, default=10)
    parser.add_argument('--timeout', type=float, default=120)
    parser.add_argument('--memory-limit-mib', type=int, default=4096)
    parser.add_argument('--output', type=Path, default=Path('test-results/parallel-libraries'))
    parser.add_argument('--only', help='Comma-separated implementation IDs, e.g. cpp:omp-static,python:thread')
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--resume', action='store_true', help='Keep completed raw records with identical experiment settings')
    args = parser.parse_args()
    if min(args.launches, args.repeats) < 1:
        parser.error('launches and repeats must be positive')
    if not math.isfinite(args.seconds) or args.seconds < 0:
        parser.error('seconds must be finite and non-negative')
    if args.profile == 'stress' and not args.seconds:
        parser.error('stress requires --seconds greater than zero')
    os.sched_setaffinity(0, {int(cpu) for cpu in args.monitor_cpus.split(',')})
    args.output.mkdir(parents=True, exist_ok=True)
    raw_dir = args.output / 'raw'
    raw_dir.mkdir(exist_ok=True)
    workloads = [('dense', 2048, 2048, 64), ('sparse-uniform', 32768, 32768, 64), ('sparse-skew', 32768, 32768, 64)]
    if args.profile == 'stress':
        workloads = [('dense', 4096, 4096, 64), ('sparse-uniform', 131072, 131072, 64), ('sparse-skew', 131072, 131072, 64)]
    if args.smoke:
        workloads = [('dense', 17, 19, 64), ('sparse-uniform', 16, 1024, 64), ('sparse-skew', 16, 1024, 64)]
    print('Computing independent integer references...', flush=True)
    references = {name: reference(name, rows, cols, nnz) for name, rows, cols, nnz in workloads}
    implementations = [impl for impl in IMPLEMENTATIONS
                       if (args.profile != 'stress' or impl[2])
                       and (not args.only or f'{impl[0]}:{impl[1]}' in args.only.split(','))]
    if not implementations:
        parser.error('no implementation selected')
    cases = []
    for language, backend, parallel in implementations:
        for workers in (list(map(int, args.workers.split(','))) if parallel else [1]):
            for workload, rows, cols, nnz in workloads:
                cases.append({'id': f'{language}:{backend}/{workload}/{workers}', 'language': language,
                              'backend': backend, 'workers': workers, 'workload': workload,
                              'rows': rows, 'cols': cols, 'nnz_per_row': nnz})
    metadata = environment(args)
    executables = {name: {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'name': path.name}
                   for name, path in [('cpp', args.cpp), ('rust', args.rust), ('python', args.python), ('python_script', args.python_script)]}
    metadata['executables'] = executables
    runtime_versions = subprocess.check_output([str(args.python.absolute()), '-c',
        "import sys,numpy,scipy,numba,llvmlite; print(sys.version); print(numpy.__version__,scipy.__version__,numba.__version__,llvmlite.__version__)"], text=True)
    metadata['python_runtime_versions'] = runtime_versions
    settings = {'cases': cases, 'references': references, 'launches': args.launches, 'repeats': args.repeats,
                'cpus': args.cpus, 'monitor_cpus': args.monitor_cpus, 'seed': args.seed,
                'pss_interval_ms': args.pss_interval_ms, 'executables': executables,
                'seconds': args.seconds, 'profile': args.profile, 'memory_limit_mib': args.memory_limit_mib,
                'python_runtime_versions': runtime_versions, 'cpu_topology': metadata['cpu_topology'],
                'sources_sha256': metadata['sources_sha256']}
    signature = hashlib.sha256(json.dumps(settings, sort_keys=True).encode()).hexdigest()
    report_path = args.output / 'results.json'
    report = {'schema': 1, 'signature': signature, 'environment': metadata, 'references': references, 'records': []}
    if args.resume and report_path.exists():
        report = json.loads(report_path.read_text())
        if report['signature'] != signature:
            raise ValueError('Source/settings changed; use a new output directory rather than mixing experiments')
    finished = {(entry['case']['id'], entry['pass'], entry['launch']) for entry in report['records']}
    schedule = []
    rng = random.Random(args.seed)
    for launch in range(args.launches):
        shuffled = list(cases)
        rng.shuffle(shuffled)
        schedule.extend((case, 'timing', launch) for case in shuffled)
    shuffled = list(cases)
    rng.shuffle(shuffled)
    schedule.extend((case, 'memory', 0) for case in shuffled)
    for index, (case, kind, launch) in enumerate(schedule):
        if (case['id'], kind, launch) in finished:
            continue
        print(f"[{index + 1}/{len(schedule)}] {kind} {case['id']} launch={launch}", flush=True)
        if available_memory_mib() < 2048:
            raise MemoryError('Less than 2 GiB memory available before benchmark start')
        executable = {'cpp': [str(args.cpp.resolve())], 'rust': [str(args.rust.resolve())],
                      'python': [str(args.python.absolute()), str(args.python_script.resolve())]}[case['language']]
        command = ['taskset', '-c', args.cpus, *executable, '--backend', case['backend'],
                   '--workers', str(case['workers']), '--workload', case['workload'],
                   '--rows', str(case['rows']), '--cols', str(case['cols']), '--nnz-per-row', str(case['nnz_per_row']),
                   '--repeats', str(args.repeats if kind == 'timing' else 2), '--seconds', str(args.seconds),
                   '--hold-ms', '0' if kind == 'timing' else '100']
        env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'BLIS_NUM_THREADS': '1',
               'OMP_NUM_THREADS': str(case['workers']), 'OMP_DYNAMIC': 'FALSE',
               'NUMBA_NUM_THREADS': str(case['workers']), 'NUMBA_THREADING_LAYER': 'omp', 'PYTHONHASHSEED': '0'}
        result = execute(command, env, kind == 'memory', raw_dir / f'{index:04d}', args.timeout, args.pss_interval_ms, args.memory_limit_mib)
        observed = result['measurement']
        expected = references[case['workload']]
        for key in ['checksum', 'output_hash', 'logical_data_bytes']:
            if str(observed[key]) != str(expected[key]):
                raise AssertionError(f"{case['id']}: {key}: {observed[key]} != {expected[key]}")
        durations = observed['samples_ms']
        expected_count = args.repeats if kind == 'timing' else 2
        if (len(durations) < expected_count or (not args.seconds and len(durations) != expected_count)
                or any(not math.isfinite(value) or value <= 0 for value in durations)):
            raise ValueError(f'Invalid timing samples: {durations}')
        if args.seconds and observed['elapsed_warm_ms'] < args.seconds * 1000:
            raise ValueError('Sustained run ended too early')
        report['records'].append({'case': case, 'pass': kind, 'launch': launch, **result})
        report_path.write_text(json.dumps(report, indent=2) + '\n')
    report['summary'] = []
    for case in cases:
        records = [entry for entry in report['records'] if entry['case']['id'] == case['id']]
        timing = [entry for entry in records if entry['pass'] == 'timing']
        memory = [entry for entry in records if entry['pass'] == 'memory'][0]
        report['summary'].append({**case,
            'warm_ms': quantiles([sample for entry in timing for sample in entry['measurement']['samples_ms']]),
            'first_ms': quantiles([entry['measurement']['first_ms'] for entry in timing]),
            'setup_ms': quantiles([entry['measurement']['setup_ms'] for entry in timing]),
            'process_wall_ms': quantiles([entry['process_wall_ms'] for entry in timing]),
            'peak_tree_pss_mib': memory['sampled_peak_tree_pss_bytes'] / 2**20,
            'max_processes': memory['max_observed_processes'],
            'throughput_calls_per_second': quantiles([len(entry['measurement']['samples_ms']) / (entry['measurement']['elapsed_warm_ms'] / 1000) for entry in timing]),
            'measurement_details': timing[0]['measurement']})
    report['environment']['loadavg_end'] = list(os.getloadavg())
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    print(f"Verified {len(cases)} configurations; results: {report_path}", flush=True)


if __name__ == '__main__':
    main()
