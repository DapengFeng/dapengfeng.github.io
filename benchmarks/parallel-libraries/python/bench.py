#!/usr/bin/env python3
"""Dense/CSR matrix-vector benchmark; real measurements, one JSON record.

Every timed call writes y=A*x in float64. Output validation is outside timing.
The first call includes lazy worker startup or JIT, later calls reuse resources.
Measure complete process lifetime and process-tree memory externally.
With --seconds, execute continuously for at least that warm duration and the
minimum --repeats count. Verify the first/last output without intervening scans.
"""
import time
PROGRAM_START = time.perf_counter()

import argparse
from array import array
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import json
import math
import multiprocessing
from multiprocessing import shared_memory
import os
import platform
import sys

BACKENDS = ("serial", "thread", "process-copy", "process-shared", "numpy-serial",
            "numpy-thread", "numba-serial", "numba-parallel")
MASK64 = (1 << 64) - 1
_SHARED = None


def spans(rows, workers):
    return [(rows * i // workers, rows * (i + 1) // workers) for i in range(workers)]


def scalar_kernel(data, start, end, cols, sparse):
    values, x, output = data["values"], data["x"], data["output"]
    if sparse:
        indptr, indices = data["indptr"], data["indices"]
        for row in range(start, end):
            total = 0.0
            for k in range(indptr[row], indptr[row + 1]):
                total += values[k] * x[indices[k]]
            output[row] = total
    else:
        for row in range(start, end):
            total = 0.0
            offset = row * cols
            for col in range(cols):
                total += values[offset + col] * x[col]
            output[row] = total


def copy_worker(task):
    start, rows, cols, sparse, payload = task
    data = {}
    for key, (kind, raw) in payload.items():
        data[key] = array(kind)
        data[key].frombytes(raw)
    data["output"] = array("d", [0]) * rows
    scalar_kernel(data, 0, rows, cols, sparse)
    return start, data["output"].tobytes()


def attach_shared(descriptors):
    global _SHARED
    handles, data = [], {}
    for key, name, kind in descriptors:
        block = shared_memory.SharedMemory(name=name)
        handles.append(block)
        data[key] = block.buf.cast(kind)
    _SHARED = (handles, data)


def shared_worker(task):
    scalar_kernel(_SHARED[1], *task)


def make_numba_kernel(numba, parallel, sparse):
    # Identical fused scalar row loops, no fastmath. prange is serial when
    # parallel=False and statically partitioned when parallel=True.
    if sparse:
        @numba.njit(parallel=parallel, cache=False)
        def kernel(values, x, output, indptr, indices, cols):
            for row in numba.prange(output.size):
                total = 0.0
                for k in range(indptr[row], indptr[row + 1]):
                    total += values[k] * x[indices[k]]
                output[row] = total
    else:
        @numba.njit(parallel=parallel, cache=False)
        def kernel(values, x, output, indptr, indices, cols):
            for row in numba.prange(output.size):
                total = 0.0
                for col in range(cols):
                    total += values[row * cols + col] * x[col]
                output[row] = total
    return kernel


def validate_output(output):
    """The binary-fraction inputs make every y*1024 an exact integer."""
    checksum, hashed = 0, 14695981039346656037
    for value in output:
        scaled = float(value) * 1024
        if not math.isfinite(scaled):
            raise RuntimeError("non-finite matrix output")
        q = int(scaled)
        if q != scaled:
            raise RuntimeError("output does not match exact binary-fraction protocol")
        checksum += q
        hashed = ((hashed ^ (q & MASK64)) * 1099511628211) & MASK64
    return str(checksum), f"{hashed:016x}"


class Benchmark:
    def __init__(self, args):
        self.args = args
        self.pool = None
        self.shared_handles = []
        self.data = {}
        self.versions = {}
        self.np = self.numba = None
        self.csr = None
        self.shards = []
        self.sparse = args.workload != "dense"
        self.ranges = spans(args.rows, args.workers)
        self.implementation = "CPython scalar loops"
        vectorized = args.backend.startswith(("numpy", "numba"))
        if vectorized:
            import numpy as np
            from threadpoolctl import threadpool_info
            self.np = np
            self.threadpool_info = threadpool_info
            self.versions["numpy"] = np.__version__
        counts = [args.cols] * args.rows
        if self.sparse:
            counts = ([784 if row < args.rows // 16 else 16 for row in range(args.rows)]
                      if args.workload == "sparse-skew" else [args.nnz_per_row] * args.rows)
        self.nnz = sum(counts)
        specs = {"values": ("d", self.nnz), "x": ("d", args.cols),
                 "output": ("d", args.rows)}
        if self.sparse:
            specs.update(indptr=("Q", args.rows + 1), indices=("I", self.nnz))
        for key, (kind, size) in specs.items():
            if args.backend == "process-shared":
                block = shared_memory.SharedMemory(create=True, size=size * array(kind).itemsize)
                self.shared_handles.append((key, kind, block))
                self.data[key] = block.buf.cast(kind)
            elif vectorized:
                self.data[key] = self.np.zeros(size, dtype={"d": "float64", "Q": "uint64", "I": "uint32"}[kind])
            else:
                self.data[key] = array(kind, [0]) * size
        for col in range(args.cols):
            self.data["x"][col] = (col % 31 - 15) / 16.0
        offset = 0
        if vectorized:
            np = self.np
            dense_columns = np.arange(args.cols, dtype=np.int64)
            for row, count in enumerate(counts):
                columns = ((row * 13 + np.arange(count, dtype=np.int64) * 17) % args.cols
                           if self.sparse else dense_columns)
                self.data["values"][offset:offset + count] = ((row * 17 + columns * 13) % 101 - 50) / 64.0
                if self.sparse:
                    self.data["indptr"][row] = offset
                    self.data["indices"][offset:offset + count] = columns
                offset += count
        else:
            for row, count in enumerate(counts):
                if self.sparse:
                    self.data["indptr"][row] = offset
                for k in range(count):
                    col = (row * 13 + k * 17) % args.cols if self.sparse else k
                    self.data["values"][offset] = ((row * 17 + col * 13) % 101 - 50) / 64.0
                    if self.sparse:
                        self.data["indices"][offset] = col
                    offset += 1
        if self.sparse:
            self.data["indptr"][-1] = self.nnz
        self.logical_data_bytes = 8 * (self.nnz + args.cols + args.rows)
        if self.sparse:
            self.logical_data_bytes += 8 * (args.rows + 1) + 4 * self.nnz
        if args.backend.startswith("numpy"):
            self.implementation = "NumPy dot; BLAS threads=1"
            if self.sparse:
                import scipy
                from scipy.sparse import csr_matrix
                self.versions["scipy"] = scipy.__version__
                self.csr = csr_matrix((self.data["values"], self.data["indices"], self.data["indptr"]),
                                      shape=(args.rows, args.cols), copy=False)
                self.data["values"], self.data["indices"], self.data["indptr"] = (
                    self.csr.data, self.csr.indices, self.csr.indptr)
                self.implementation = "SciPy CSR matvec; contiguous row slabs"
                if args.backend == "numpy-thread":
                    for start, end in self.ranges:
                        lo, hi = self.csr.indptr[start], self.csr.indptr[end]
                        # Data and column indices are views; only the short row
                        # pointers are rebased/copied once during setup.
                        self.shards.append(csr_matrix((self.csr.data[lo:hi], self.csr.indices[lo:hi],
                                                       self.csr.indptr[start:end + 1] - lo),
                                                      shape=(end - start, args.cols), copy=False))
            else:
                self.matrix = self.data["values"].reshape(args.rows, args.cols)
        if args.backend.startswith("numba"):
            import numba
            self.numba = numba
            self.versions["numba"] = numba.__version__
            self.kernel = make_numba_kernel(numba, args.backend == "numba-parallel", self.sparse)
            self.empty_indices = self.np.empty(0, dtype=self.np.uint64)
            self.implementation = "Numba fused row loops; fastmath=False"
        if args.backend in ("thread", "numpy-thread"):
            self.pool = ThreadPoolExecutor(max_workers=args.workers)
        if args.backend == "process-copy":
            self.pool = ProcessPoolExecutor(max_workers=args.workers,
                                           mp_context=multiprocessing.get_context("spawn"))
        if args.backend == "process-shared":
            descriptors = [(key, block.name, kind) for key, kind, block in self.shared_handles]
            self.pool = ProcessPoolExecutor(max_workers=args.workers,
                                           mp_context=multiprocessing.get_context("spawn"),
                                           initializer=attach_shared, initargs=(descriptors,))

    def numpy_chunk(self, index, start, end):
        if self.sparse:
            # SciPy returns a fresh vector. Its allocation and the output write
            # are part of the measured call.
            matrix = self.shards[index] if self.shards else self.csr
            self.data["output"][start:end] = matrix @ self.data["x"]
        else:
            self.np.dot(self.matrix[start:end], self.data["x"], out=self.data["output"][start:end])

    def copy_tasks(self):
        args, data = self.args, self.data
        for start, end in self.ranges:
            lo, hi = ((data["indptr"][start], data["indptr"][end]) if self.sparse
                      else (start * args.cols, end * args.cols))
            payload = {"values": ("d", data["values"][lo:hi].tobytes()),
                       "x": ("d", data["x"].tobytes())}
            if self.sparse:
                payload["indices"] = ("I", data["indices"][lo:hi].tobytes())
                payload["indptr"] = ("Q", array("Q", (value - lo for value in data["indptr"][start:end + 1])).tobytes())
            yield start, end - start, args.cols, self.sparse, payload

    def run(self):
        args = self.args
        if args.backend == "process-copy":
            for start, payload in self.pool.map(copy_worker, self.copy_tasks()):
                values = array("d")
                values.frombytes(payload)
                self.data["output"][start:start + len(values)] = values
        elif args.backend == "process-shared":
            list(self.pool.map(shared_worker, ((start, end, args.cols, self.sparse) for start, end in self.ranges)))
        elif args.backend.startswith("numba"):
            self.kernel(self.data["values"], self.data["x"], self.data["output"],
                        self.data.get("indptr", self.empty_indices),
                        self.data.get("indices", self.empty_indices), args.cols)
        elif args.backend == "numpy-serial":
            self.numpy_chunk(0, 0, args.rows)
        elif args.backend == "numpy-thread":
            futures = [self.pool.submit(self.numpy_chunk, i, start, end) for i, (start, end) in enumerate(self.ranges)]
            for future in futures:
                future.result()
        elif args.backend == "thread":
            futures = [self.pool.submit(scalar_kernel, self.data, start, end, args.cols, self.sparse)
                       for start, end in self.ranges]
            for future in futures:
                future.result()
        else:
            scalar_kernel(self.data, 0, args.rows, args.cols, self.sparse)

    def storage_info(self):
        if self.np:
            result = {key: {"dtype": str(value.dtype), "bytes": value.nbytes} for key, value in self.data.items()}
            result["shard_row_pointer_bytes"] = sum(matrix.indptr.nbytes for matrix in self.shards)
            return result
        return {key: {"format": value.format if isinstance(value, memoryview) else value.typecode,
                      "bytes": len(value) * value.itemsize} for key, value in self.data.items()}

    def close(self):
        if self.pool is not None:
            self.pool.shutdown(wait=True, cancel_futures=True)
        for key, _kind, block in self.shared_handles:
            self.data[key].release()
            block.close()
            block.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", required=True, choices=BACKENDS)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--workload", required=True, choices=("dense", "sparse-uniform", "sparse-skew"))
    parser.add_argument("--rows", type=int, required=True)
    parser.add_argument("--cols", type=int, required=True)
    parser.add_argument("--nnz-per-row", type=int, default=64)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seconds", type=float, default=0)
    parser.add_argument("--hold-ms", type=int, default=0)
    args = parser.parse_args()
    if min(args.workers, args.rows, args.cols, args.nnz_per_row, args.repeats) < 1 or args.hold_ms < 0:
        parser.error("sizes, workers and repeats must be positive; hold-ms must be nonnegative")
    if not math.isfinite(args.seconds) or args.seconds < 0:
        parser.error("seconds must be finite and nonnegative")
    if args.workers > args.rows:
        parser.error("workers must not exceed rows")
    if args.backend in ("serial", "numpy-serial", "numba-serial") and args.workers != 1:
        parser.error("serial backends require --workers 1")
    if args.workload == "sparse-skew" and (args.rows % 16 or args.nnz_per_row != 64 or args.cols < 784):
        parser.error("sparse-skew requires rows divisible by 16, nnz-per-row=64 and cols>=784")
    if args.workload != "dense" and (args.nnz_per_row > args.cols or args.cols % 17 == 0):
        parser.error("CSR requires nnz-per-row<=cols and cols not divisible by 17 for unique columns")
    if args.cols > 0xFFFFFFFF or array("I").itemsize != 4 or array("Q").itemsize != 8:
        parser.error("requires uint32 column indices and uint64 row pointers")
    for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    os.environ["OMP_NUM_THREADS"] = str(args.workers)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["NUMBA_NUM_THREADS"] = str(args.workers)
    os.environ["NUMBA_THREADING_LAYER"] = "omp"
    bench = Benchmark(args)
    try:
        setup_ms = (time.perf_counter() - PROGRAM_START) * 1000
        start = time.perf_counter()
        bench.run()
        first_ms = (time.perf_counter() - start) * 1000
        expected = validate_output(bench.data["output"])
        samples = []
        warm_start = time.perf_counter()
        while len(samples) < args.repeats or time.perf_counter() - warm_start < args.seconds:
            start = time.perf_counter()
            bench.run()
            samples.append((time.perf_counter() - start) * 1000)
            if not args.seconds and validate_output(bench.data["output"]) != expected:
                raise RuntimeError("output checksum/hash changed between calls")
        elapsed_warm_ms = (time.perf_counter() - warm_start) * 1000
        if args.seconds and validate_output(bench.data["output"]) != expected:
            raise RuntimeError("final stress output checksum/hash differs from first call")
        result = {
            "backend": args.backend, "language": "python", "workers": args.workers,
            "workload": args.workload, "rows": args.rows, "cols": args.cols,
            "nnz_per_row": args.nnz_per_row, "nnz": bench.nnz,
            "repeats": args.repeats, "setup_ms": setup_ms, "first_ms": first_ms,
            "stress_seconds_requested": args.seconds, "elapsed_warm_ms": elapsed_warm_ms,
            "iterations": len(samples), "validation_scope": "first,last" if args.seconds else "all",
            "samples_ms": samples, "checksum": expected[0], "output_hash": expected[1],
            "validated": True, "python": platform.python_version(), "versions": bench.versions,
            "python_implementation": platform.python_implementation(),
            "python_compiler": platform.python_compiler(), "python_full": sys.version.replace("\n", " "),
            "implementation": bench.implementation, "partition": "equal contiguous rows",
            "logical_data_bytes": bench.logical_data_bytes, "storage": bench.storage_info(),
        }
        if args.backend == "numba-parallel":
            result["threading_layer"] = bench.numba.threading_layer()
            result["numba_threads"] = bench.numba.get_num_threads()
        if args.backend.startswith("process"):
            result["start_method"] = "spawn"
        if bench.np:
            # Report the libraries actually loaded, not just the environment
            # requests. Do not publish machine-specific shared-library paths.
            runtime_fields = ("user_api", "internal_api", "num_threads", "prefix",
                              "version", "threading_layer", "architecture")
            result["threadpool_runtimes"] = [
                {key: library[key] for key in runtime_fields if key in library}
                for library in bench.threadpool_info()]
        print(json.dumps(result, separators=(",", ":")), flush=True)
        if args.hold_ms:
            time.sleep(args.hold_ms / 1000.0)
    finally:
        bench.close()


if __name__ == "__main__":
    main()
