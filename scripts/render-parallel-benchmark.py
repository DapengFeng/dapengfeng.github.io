#!/usr/bin/env python3
"""Render measured parallel-matvec reports, never run benchmark processes.

Use the pinned benchmark Python with Matplotlib available. Both input reports
must be complete: this script rejects missing cases, unvalidated outputs and
source hashes that differ from the measured files. It copies full evidence,
exports all summary rows, and builds a reproducible source archive.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import zipfile


REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / "benchmarks" / "parallel-libraries"
BACKGROUND = "#11202a"
FOREGROUND = "#e8edef"
MUTED = "#a9bbc5"
GRID = "#38505e"
CONTEXT_NOTE = ""
COLORS = {"cpp": "#9bd8f0", "rust": "#e9b47a", "python": "#c0e49b"}
METHODS = [
    ("cpp", "serial", "C++ · serial", False),
    ("cpp", "std-thread", "C++ · std::thread", True),
    ("cpp", "omp-static", "C++ · OpenMP static", True),
    ("cpp", "omp-dynamic", "C++ · OpenMP dynamic", True),
    ("cpp", "tbb", "C++ · oneTBB", True),
    ("rust", "serial", "Rust · serial", False),
    ("rust", "std-scoped", "Rust · scoped threads", True),
    ("rust", "rayon", "Rust · Rayon", True),
    ("python", "serial", "Python · scalar serial", False),
    ("python", "thread", "Python · thread pool", True),
    ("python", "process-copy", "Python · processes / copy", True),
    ("python", "process-shared", "Python · processes / shared", True),
    ("python", "numpy-serial", "Python · NumPy serial", False),
    ("python", "numpy-thread", "Python · NumPy + threads", True),
    ("python", "numba-serial", "Python · Numba serial", False),
    ("python", "numba-parallel", "Python · Numba prange", True),
]
LABELS = {(language, backend): label for language, backend, label, _ in METHODS}
SERIAL = {(language, backend) for language, backend, _label, parallel in METHODS if not parallel}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_positive(value, description):
    if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{description} must be finite and positive: {value!r}")


def read_report(path):
    report = json.loads(path.read_text())
    if report.get("schema") != 1 or not report.get("summary"):
        raise ValueError(f"{path}: incomplete or unsupported report")
    ids = [row["id"] for row in report["summary"]]
    if len(ids) != len(set(ids)):
        raise ValueError(f"{path}: duplicate summary case IDs")
    by_case = {case_id: [] for case_id in ids}
    for record in report["records"]:
        measurement = record["measurement"]
        if not measurement.get("validated"):
            raise ValueError(f"{path}: unvalidated record {record['case']['id']}")
        reference = report["references"][record["case"]["workload"]]
        for key in ("checksum", "output_hash", "logical_data_bytes"):
            if str(measurement[key]) != str(reference[key]):
                raise ValueError(f"{path}: reference mismatch for {key}")
        if record["case"]["id"] not in by_case:
            raise ValueError(f"{path}: raw record missing from summary")
        by_case[record["case"]["id"]].append(record)
    for row in report["summary"]:
        records = by_case[row["id"]]
        if sum(entry["pass"] == "memory" for entry in records) != 1:
            raise ValueError(f"{path}: expected one memory pass for {row['id']}")
        if sum(entry["pass"] == "timing" for entry in records) != report["environment"]["fresh_processes"]:
            raise ValueError(f"{path}: missing timing launches for {row['id']}")
        finite_positive(row["warm_ms"]["median"], f"{row['id']} median")
        finite_positive(row["peak_tree_pss_mib"], f"{row['id']} PSS")
        if row["warm_ms"]["q1"] > row["warm_ms"]["median"] or row["warm_ms"]["q3"] < row["warm_ms"]["median"]:
            raise ValueError(f"{path}: invalid quartiles for {row['id']}")
    report["_by_case"] = by_case
    report["_by_id"] = {row["id"]: row for row in report["summary"]}
    return report


def select(report, workload, workers, parallel_only=False):
    selected = []
    for language, backend, _label, parallel in METHODS:
        if parallel_only and not parallel:
            continue
        count = workers if parallel else 1
        case_id = f"{language}:{backend}/{workload}/{count}"
        if case_id not in report["_by_id"]:
            raise ValueError(f"missing requested chart case {case_id}")
        selected.append(report["_by_id"][case_id])
    return selected


def label(row, sparse=False):
    value = LABELS[(row["language"], row["backend"])]
    return value.replace("NumPy", "SciPy CSR") if sparse else value


def number(value):
    """Readable measured values, with no implication of excessive precision."""
    if value >= 1000:
        return f"{value:,.0f}"
    if value >= 100:
        return f"{value:.1f}"
    if value >= 10:
        return f"{value:.2f}"
    if value >= 1:
        return f"{value:.3f}"
    return f"{value:.3g}"


def init_plotting():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 12,
        "text.color": FOREGROUND, "axes.labelcolor": MUTED,
        "xtick.color": MUTED, "ytick.color": FOREGROUND,
        "axes.facecolor": BACKGROUND, "figure.facecolor": BACKGROUND,
        "savefig.facecolor": BACKGROUND, "axes.edgecolor": GRID,
        "svg.fonttype": "none", "svg.hashsalt": "feng-parallel-matvec-v1",
        "axes.unicode_minus": True, "figure.dpi": 120, "savefig.dpi": 180,
    })
    return plt


def configure_axis(ax, rows, show_labels=True, sparse=False):
    ax.set_ylim(len(rows) - 0.45, -0.85)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([label(row, sparse) for row in rows] if show_labels else [], fontsize=13)
    ax.tick_params(axis="y", length=0, pad=13)
    ax.tick_params(axis="x", length=0, pad=7, labelsize=11)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color=GRID, linewidth=0.7, alpha=0.65)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for index, row in enumerate(rows):
        if index and row["language"] != rows[index - 1]["language"]:
            ax.axhline(index - 0.5, color=GRID, linewidth=0.9, alpha=0.8)
    if show_labels:
        for tick, row in zip(ax.get_yticklabels(), rows):
            tick.set_color(COLORS[row["language"]])


def log_ticks(ax):
    from matplotlib.ticker import LogLocator, NullFormatter, FuncFormatter
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=7))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _position: f"{value:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="x", which="minor", length=0)


def point_iqr(ax, row, y, field):
    stats, color = row[field], COLORS[row["language"]]
    median = stats["median"]
    ax.hlines(y, stats["q1"], stats["q3"], color=color, linewidth=2.2)
    ax.vlines([stats["q1"], stats["q3"]], y - 0.08, y + 0.08, color=color, linewidth=1.3)
    ax.scatter([median], [y], color=color, s=34, zorder=3, edgecolors=BACKGROUND, linewidths=0.7)
    ax.annotate(number(median), (max(median, stats["q3"]), y), xytext=(7, 0),
                textcoords="offset points", va="center", fontsize=12, color=FOREGROUND)


def figure_header(fig, title, subtitle, environment):
    fig.text(0.045, 0.956, title, fontsize=19, fontweight="semibold", va="top")
    fig.text(0.045, 0.907, subtitle, fontsize=12, color=MUTED, va="top")
    cpu = next((item["data"] for item in environment["cpu"]["lscpu"] if item["field"] == "Model name:"), environment["machine"])
    cpu = cpu.replace("Intel(R)", "Intel").replace("Core(TM)", "Core")
    context = " · ".join(filter(None, [CONTEXT_NOTE, cpu, environment["measured_at"][:10]]))
    fig.text(0.045, 0.862, context, fontsize=10.5, color="#e9b47a", va="top")


def save_figure(fig, output, name, title, description):
    fig.savefig(output / f"{name}.svg", metadata={"Date": None, "Title": title, "Description": description})
    fig.savefig(output / f"{name}.png", metadata={"Title": title, "Description": description})


def two_panel(plt, rows, field, title, subtitle, xlabel, footer, output, name, environment):
    height = 7.4 if len(rows) > 13 else 6.3
    fig, axes = plt.subplots(1, 2, figsize=(12.2, height), sharey=True,
                             gridspec_kw={"width_ratios": [1.1, 1], "wspace": 0.16})
    fig.subplots_adjust(left=0.30, right=0.972, top=0.74, bottom=0.23)
    left, right = axes
    configure_axis(left, rows)
    # sharey also shares the formatter; hide labels without replacing it.
    configure_axis(right, rows)
    right.tick_params(axis="y", labelleft=False)
    log_ticks(left)
    left.set_xlim(min(row[field]["q1"] for row in rows) / 1.6,
                  max(row[field]["q3"] for row in rows) * 6)
    left.set_xlabel(xlabel, fontsize=12, labelpad=12)
    left.set_title("Warm execution", fontsize=13, loc="left", pad=18)
    right.set_title("Sampled peak process-tree PSS", fontsize=12, loc="left", pad=18)
    peak = max(row["peak_tree_pss_mib"] for row in rows)
    right.set_xlim(0, peak * 1.27)
    right.set_xlabel("Memory [MiB] · linear scale", fontsize=12, labelpad=12)
    for index, row in enumerate(rows):
        point_iqr(left, row, index, field)
        value = row["peak_tree_pss_mib"]
        right.barh(index, value, height=0.34, color=COLORS[row["language"]], alpha=0.78)
        right.annotate(f"{value:,.1f}", (value, index), xytext=(6, 0),
                       textcoords="offset points", va="center", fontsize=12)
    figure_header(fig, title, subtitle, environment)
    fig.text(0.045, 0.045, footer, fontsize=10.4, color=MUTED, va="bottom", linespacing=1.65)
    save_figure(fig, output, name, title, footer)
    plt.close(fig)


def render_dense(plt, baseline, workers, output):
    rows = select(baseline, "dense", workers)
    sample = rows[0]
    env = baseline["environment"]
    subtitle = (f"{sample['rows']:,} × {sample['cols']:,} · float64 · "
                f"parallel p={workers}, serial p=1 · {env['fresh_processes']} fresh launches")
    footer = ("Time: median and interquartile range of warm calls; logarithmic axis. Startup / first call excluded.\n"
              f"Memory: separate lifecycle pass, sampled every {env['pss_interval_ms']:g} ms; includes live children and startup.")
    two_panel(plt, rows, "warm_ms", "Dense matrix–vector product", subtitle,
              "Time [ms] · logarithmic scale", footer, output, "dense-time-memory", env)


def render_sparse(plt, report, workers, output, sustained=False):
    uniform = select(report, "sparse-uniform", workers, parallel_only=sustained)
    skew = select(report, "sparse-skew", workers, parallel_only=sustained)
    if any(a["rows"] != b["rows"] or a["cols"] != b["cols"]
           or a["measurement_details"]["nnz"] != b["measurement_details"]["nnz"]
           for a, b in zip(uniform, skew)):
        raise ValueError("sparse chart requires equal shapes and stored-entry counts")
    fig = plt.figure(figsize=(12.2, 6.3 if sustained else 7.4))
    ax = fig.add_axes([0.30, 0.23, 0.47, 0.53])
    numbers = fig.add_axes([0.79, 0.23, 0.19, 0.53])
    configure_axis(ax, uniform, sparse=True)
    log_ticks(ax)
    values = [row["warm_ms"]["median"] for row in uniform + skew]
    ax.set_xlim(min(values) / 1.8, max(values) * 1.8)
    ax.set_xlabel("Warm matvec time [ms] · logarithmic scale", fontsize=12, labelpad=13)
    numbers.set_xlim(0, 1)
    numbers.set_ylim(ax.get_ylim())
    numbers.axis("off")
    for index, (first, second) in enumerate(zip(uniform, skew)):
        color = COLORS[first["language"]]
        one, two = first["warm_ms"]["median"], second["warm_ms"]["median"]
        ax.hlines(index, min(one, two), max(one, two), linewidth=1.5, color=color, alpha=0.7)
        ax.scatter(one, index, marker="o", s=37, color=color, zorder=3)
        ax.scatter(two, index, marker="D", s=38, facecolors=BACKGROUND,
                   edgecolors=color, linewidths=1.4, zorder=3)
        numbers.text(0.42, index, number(one), ha="right", va="center", fontsize=11.5)
        numbers.text(0.98, index, number(two), ha="right", va="center", fontsize=11.5)
    numbers.text(0.42, 1.075, "Uniform", transform=numbers.transAxes,
                 ha="right", va="center", color=MUTED, fontsize=11)
    numbers.text(0.98, 1.075, "Skew", transform=numbers.transAxes,
                 ha="right", va="center", color=MUTED, fontsize=11)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=FOREGROUND, marker="o", linestyle="none", label="Uniform"),
                       Line2D([], [], color=FOREGROUND, marker="D", markerfacecolor=BACKGROUND,
                              linestyle="none", label="Skew")],
              loc="lower left", bbox_to_anchor=(-0.01, 1.02), ncol=2,
              frameon=False, fontsize=11, borderaxespad=0, handletextpad=0.35)
    sample = uniform[0]
    subtitle = (f"{sample['rows']:,} × {sample['cols']:,} CSR · "
                f"{sample['measurement_details']['nnz']:,} stored entries each · "
                + (f"p={workers} · ≥{report['environment']['stress_seconds']:g} s per launch" if sustained
                   else f"parallel p={workers}, serial p=1"))
    title = "Sustained execution: uniform vs concentrated stored entries" if sustained else "Uniform vs concentrated stored entries"
    figure_header(fig, title, subtitle, report["environment"])
    footer = ("Dots show warm medians. Uniform: 64 stored entries per row. Skew: 1/16 of rows have 784, the rest 16.\n"
              "Equal entry counts; row lengths and column access patterns differ.")
    if sustained:
        footer = ("Dots show medians during continuous warm execution; first and last outputs are checked.\n"
                  "Both CSR matrices contain the same number of stored entries. Memory and tail times for both are in the full CSV.")
    fig.text(0.045, 0.035, footer, fontsize=10.4, color=MUTED, va="bottom", linespacing=1.65)
    name = "stress-sparse-uniform-skew" if sustained else "sparse-uniform-skew"
    save_figure(fig, output, name, title, footer)
    plt.close(fig)


def render_stress(plt, stress, workers, output):
    rows = select(stress, "dense", workers, parallel_only=True)
    sample, env = rows[0], stress["environment"]
    duration = env["stress_seconds"]
    if duration <= 0:
        raise ValueError("stress figure requires a report with positive sustained duration")
    subtitle = (f"{sample['rows']:,} × {sample['cols']:,} · p={workers} · "
                f"at least {duration:g} s of continuous warm execution per launch")
    footer = (f"Throughput: calls / warm wall time; {env['fresh_processes']} timed window(s) per case. Higher means more completed calls.\n"
              f"Memory: separate lifecycle pass, {env['pss_interval_ms']:g} ms PSS sampling. Tail times and sample counts are in the CSV.")
    two_panel(plt, rows, "throughput_calls_per_second", "Sustained execution on the larger dense matrix",
              subtitle, "Throughput [calls/s] · logarithmic scale", footer, output, "stress-throughput-memory", env)


def flattened(profile, row):
    details = row["measurement_details"]
    item = {"profile": profile, "id": row["id"], "language": row["language"],
            "backend": row["backend"], "label": label(row, row["workload"] != "dense"),
            "workers": row["workers"], "workload": row["workload"], "rows": row["rows"],
            "cols": row["cols"], "nnz_per_row": row["nnz_per_row"], "nnz": details["nnz"],
            "logical_data_bytes": details["logical_data_bytes"],
            "peak_tree_pss_mib": row["peak_tree_pss_mib"], "max_processes": row["max_processes"],
            "stress_seconds_requested": details.get("stress_seconds_requested", 0),
            "checksum": details["checksum"], "output_hash": details["output_hash"]}
    for field in ("warm_ms", "first_ms", "setup_ms", "process_wall_ms", "throughput_calls_per_second"):
        for stat, value in row[field].items():
            item[f"{field}_{stat}"] = value
    return item


def publish_data(reports, paths, output, baseline_workers, stress_workers):
    flattened_rows = [flattened(profile, row) for profile, report in reports.items() for row in report["summary"]]
    with (output / "results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flattened_rows[0]))
        writer.writeheader()
        writer.writerows(flattened_rows)
    for profile, path in paths.items():
        shutil.copyfile(path, output / f"{profile}.json")
    summary = {"schema": 1, "figure_context_note": CONTEXT_NOTE, "chart_workers": {"baseline": baseline_workers, "stress": stress_workers},
               "measurement_notes": {
                   "time": "Warm calls only; first call and setup are separate fields.",
                   "memory": "Separate pass; sampled peak contemporaneous process-tree PSS across the complete lifecycle.",
                   "tails": "p95 omitted below 20 warm samples; p99 omitted below 100 samples.",
                   "scope": "Measurements on one recorded host, not general language rankings."},
               "profiles": {profile: {"signature": report["signature"],
                                      "environment": report["environment"],
                                      "references": report["references"],
                                      "cases": {row["id"]: flattened(profile, row) for row in report["summary"]}}
                            for profile, report in reports.items()}}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")


def package_sources(reports, output):
    for profile, report in reports.items():
        for relative, expected in report["environment"]["sources_sha256"].items():
            path = (SOURCE / relative).resolve()
            if not path.is_relative_to(SOURCE.resolve()) or not path.is_file() or digest(path) != expected:
                raise ValueError(f"{profile}: measured source changed or missing: {relative}")
    excluded = {"__pycache__", ".venv", ".build", "target", ".git", "build", "dist", "results", ".cache"}
    suffixes = {".py", ".cpp", ".h", ".hpp", ".rs", ".toml", ".lock", ".sh", ".txt", ".md", ".json"}
    sources = [path for path in SOURCE.rglob("*") if path.is_file()
               and not any(part in excluded for part in path.relative_to(SOURCE).parts)
               and (path.suffix in suffixes or path.name == ".gitignore")]
    sources.append(Path(__file__).resolve())
    contents = {str(path.relative_to(REPO)): path.read_bytes() for path in sorted(sources)}
    manifest = {"sha256": {name: hashlib.sha256(data).hexdigest() for name, data in contents.items()},
                "measured_report_signatures": {profile: report["signature"] for profile, report in reports.items()}}
    contents["SOURCE-MANIFEST.json"] = (json.dumps(manifest, indent=2) + "\n").encode()
    with zipfile.ZipFile(output / "source.zip", "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(contents.items()):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data)


def main():
    global CONTEXT_NOTE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=REPO / "test-results/parallel-baseline/results.json")
    parser.add_argument("--stress", type=Path, default=REPO / "test-results/parallel-stress/results.json")
    parser.add_argument("--output", type=Path, default=REPO / "content/assets/parallel-libraries")
    parser.add_argument("--context-note", default="", help="Short documented host condition shown on exported figures")
    parser.add_argument("--baseline-workers", type=int, default=4)
    parser.add_argument("--stress-workers", type=int, default=16)
    args = parser.parse_args()
    CONTEXT_NOTE = args.context_note
    if min(args.baseline_workers, args.stress_workers) < 1:
        parser.error("worker counts must be positive")
    paths = {"baseline": args.baseline, "stress": args.stress}
    reports = {profile: read_report(path) for profile, path in paths.items()}
    # Validate the requested selection and source evidence before writing.
    select(reports["baseline"], "dense", args.baseline_workers)
    select(reports["baseline"], "sparse-uniform", args.baseline_workers)
    select(reports["baseline"], "sparse-skew", args.baseline_workers)
    select(reports["stress"], "dense", args.stress_workers, parallel_only=True)
    select(reports["stress"], "sparse-uniform", args.stress_workers, parallel_only=True)
    select(reports["stress"], "sparse-skew", args.stress_workers, parallel_only=True)
    args.output.mkdir(parents=True, exist_ok=True)
    package_sources(reports, args.output)
    publish_data(reports, paths, args.output, args.baseline_workers, args.stress_workers)
    plt = init_plotting()
    render_dense(plt, reports["baseline"], args.baseline_workers, args.output)
    render_sparse(plt, reports["baseline"], args.baseline_workers, args.output)
    render_stress(plt, reports["stress"], args.stress_workers, args.output)
    render_sparse(plt, reports["stress"], args.stress_workers, args.output, sustained=True)
    print(f"Rendered 4 figures (SVG + PNG), 2 full reports, CSV, summary and source archive: {args.output}")


if __name__ == "__main__":
    main()
