#!/usr/bin/env python3
"""Independent exact-arithmetic oracle for the matrix-vector benchmarks.

This module deliberately imports no measured kernel or numerical library. It
reconstructs each row directly from the public input formulas and accumulates
integer numerators. Since A has denominator 64 and x has denominator 16, each
row result is exactly q / 1024. The checksum is sum(q), not sum(y).

The output hash is value-wise FNV-1a: xor each signed q converted modulo 2**64
into the accumulator, then multiply by the FNV prime modulo 2**64. It is not
the byte-wise FNV hash of a serialized array. This checks output order as well
as values without depending on host byte order.

Logical storage describes float64 values/vectors, uint32 CSR column indices,
and uint64 CSR row offsets. It does not estimate runtime resident memory.
"""

import argparse
import json
import unittest


MASK64 = (1 << 64) - 1
FNV_OFFSET = 14695981039346656037
FNV_PRIME = 1099511628211
WORKLOADS = ("dense", "sparse-uniform", "sparse-skew")


def _positive_integer(name, value):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def reference(workload, rows, cols, nnz_per_row=64):
    """Return exact expected output and logical-storage metadata.

    Sparse rows contain ``nnz_per_row`` entries in the uniform case. In the
    skewed case, the first rows/16 rows contain 49*nnz_per_row/4 entries and
    the rest contain nnz_per_row/4, preserving the total stored-entry count.
    Explicit zero matrix entries remain stored, as in the measured kernels.
    Repeated column indices, if produced by a particular cols value, remain
    separate CSR entries and are each included in the dot product.
    """
    if workload not in WORKLOADS:
        raise ValueError(f"workload must be one of {', '.join(WORKLOADS)}")
    _positive_integer("rows", rows)
    _positive_integer("cols", cols)
    _positive_integer("nnz_per_row", nnz_per_row)
    if workload == "sparse-uniform" and cols < nnz_per_row:
        raise ValueError("sparse-uniform requires cols >= nnz_per_row")
    if workload == "sparse-skew":
        if rows % 16 or nnz_per_row % 4:
            raise ValueError("sparse-skew requires rows divisible by 16 and nnz_per_row by 4")
        if cols < nnz_per_row * 49 // 4:
            raise ValueError("sparse-skew requires cols >= the heavy-row entry count")

    checksum = 0
    output_hash = FNV_OFFSET
    stored_entries = 0
    nonzero_entries = 0
    max_abs_q = 0
    for row in range(rows):
        if workload == "dense":
            row_entries = cols
        elif workload == "sparse-uniform":
            row_entries = nnz_per_row
        else:
            row_entries = (nnz_per_row * 49 // 4 if row < rows // 16
                           else nnz_per_row // 4)
        q = 0
        for entry in range(row_entries):
            col = entry if workload == "dense" else (row * 13 + entry * 17) % cols
            a_numerator = (row * 17 + col * 13) % 101 - 50
            x_numerator = col % 31 - 15
            q += a_numerator * x_numerator
            nonzero_entries += a_numerator != 0
        stored_entries += row_entries
        checksum += q
        max_abs_q = max(max_abs_q, abs(q))
        output_hash = ((output_hash ^ (q & MASK64)) * FNV_PRIME) & MASK64

    if workload == "dense":
        logical_data_bytes = 8 * (stored_entries + cols + rows)
    else:
        logical_data_bytes = 12 * stored_entries + 8 * (rows + 1 + cols + rows)
    return {
        "checksum": str(checksum),
        "output_hash": f"{output_hash:016x}",
        "logical_data_bytes": logical_data_bytes,
        "stored_entries": stored_entries,
        "nonzero_entries": nonzero_entries,
        "max_abs_q": max_abs_q,
        "output_scale": 1024,
    }


class ReferenceFixtures(unittest.TestCase):
    """Small hand-derived fixtures; no production implementation is imported."""

    def test_dense_two_rows(self):
        # q = [750 + 518 + 312, 495 + 280 + 91] = [1580, 866].
        result = reference("dense", 2, 3)
        self.assertEqual(result["checksum"], "2446")
        self.assertEqual(result["max_abs_q"], 1580)
        self.assertEqual(result["logical_data_bytes"], 88)
        self.assertEqual(result["stored_entries"], 6)
        self.assertEqual(result["nonzero_entries"], 6)
        h = ((FNV_OFFSET ^ 1580) * FNV_PRIME) & MASK64
        h = ((h ^ 866) * FNV_PRIME) & MASK64
        self.assertEqual(result["output_hash"], f"{h:016x}")

    def test_sparse_uniform_two_rows(self):
        # Columns [0, 2], [3, 0] give q = [750 + 312, -72 + 495].
        result = reference("sparse-uniform", 2, 5, 2)
        self.assertEqual(result["checksum"], "1485")
        self.assertEqual(result["max_abs_q"], 1062)
        self.assertEqual(result["logical_data_bytes"], 128)
        self.assertEqual(result["stored_entries"], 4)

    def test_stored_zeros_are_counted(self):
        # Row 0, column 66 has numerator (13*66)%101 - 50 = 0.
        result = reference("dense", 1, 67)
        self.assertEqual(result["stored_entries"], 67)
        self.assertEqual(result["nonzero_entries"], 66)

    def test_full_row_fixtures(self):
        # Verified separately by materializing rational-valued dense matrices,
        # scatter-adding CSR entries and multiplying with fractions.Fraction.
        # The negative skew checksum also exercises signed-to-u64 hashing.
        cases = [
            ("dense", 17, 19, "610", "d0ef9066cad40d2d", 2872, 323, 319),
            ("sparse-uniform", 16, 64, "1056", "54e21a89de1693c5", 1544, 64, 63),
            ("sparse-skew", 16, 64, "-187", "38de260ce1edda80", 1544, 64, 63),
        ]
        for name, rows, cols, checksum, digest, size, stored, nonzero in cases:
            with self.subTest(workload=name):
                result = reference(name, rows, cols, 4)
                self.assertEqual(result["checksum"], checksum)
                self.assertEqual(result["output_hash"], digest)
                self.assertEqual(result["logical_data_bytes"], size)
                self.assertEqual(result["stored_entries"], stored)
                self.assertEqual(result["nonzero_entries"], nonzero)

    def test_invalid_inputs(self):
        invalid = [
            ("unknown", 16, 64, 4), ("dense", 0, 4, 4),
            ("dense", 1, True, 4), ("sparse-uniform", 1, 3, 4),
            ("sparse-skew", 15, 64, 4), ("sparse-skew", 16, 64, 5),
            ("sparse-skew", 16, 48, 4),
        ]
        for args in invalid:
            with self.subTest(args=args), self.assertRaises(ValueError):
                reference(*args)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--workload", choices=WORKLOADS)
    parser.add_argument("--rows", type=int)
    parser.add_argument("--cols", type=int)
    parser.add_argument("--nnz-per-row", type=int, default=64)
    args = parser.parse_args()
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(ReferenceFixtures)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        raise SystemExit(0 if result.wasSuccessful() else 1)
    if args.workload is None or args.rows is None or args.cols is None:
        parser.error("--workload, --rows and --cols are required unless --self-test is used")
    try:
        print(json.dumps(reference(args.workload, args.rows, args.cols, args.nnz_per_row)))
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
