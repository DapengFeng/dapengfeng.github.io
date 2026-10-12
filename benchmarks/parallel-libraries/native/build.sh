#!/usr/bin/env bash
set -euo pipefail

NATIVE_SOURCE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PARALLEL_BUILD_DIR="${PARALLEL_BUILD_DIR:-$NATIVE_SOURCE_DIR/.build}"
mkdir -p "$PARALLEL_BUILD_DIR"
cpp_flags=(-O3 -march=native -fno-fast-math -ffp-contract=off -std=c++20 -pthread -fopenmp -Wall -Wextra)
tbb_flags=()
if [[ "${ENABLE_TBB:-1}" == 1 ]]; then
  cpp_flags+=(-DBENCH_WITH_TBB)
  if [[ -n "${TBB_ROOT:-}" ]]; then
    cpp_flags+=(-I"$TBB_ROOT/include")
    TBB_LIB_DIR="${TBB_LIB_DIR:-$TBB_ROOT/lib/x86_64-linux-gnu}"
    tbb_flags+=(-L"$TBB_LIB_DIR" -Wl,-rpath,"$TBB_LIB_DIR")
  fi
  tbb_flags+=(-ltbb)
fi
"${CXX:-g++}" "${cpp_flags[@]}" "$NATIVE_SOURCE_DIR/parallel.cpp" "${tbb_flags[@]}" -o "$PARALLEL_BUILD_DIR/parallel-cpp"
RUSTFLAGS="${RUSTFLAGS:--C target-cpu=native}" cargo build --release --locked --manifest-path "$NATIVE_SOURCE_DIR/rust/Cargo.toml" --target-dir "$PARALLEL_BUILD_DIR/rust-target"
cp "$PARALLEL_BUILD_DIR/rust-target/release/parallel-libraries-bench" "$PARALLEL_BUILD_DIR/parallel-rust"
"${CXX:-g++}" --version | head -1
rustc --version
cargo --version
printf 'Binaries: %s/parallel-cpp and %s/parallel-rust\n' "$PARALLEL_BUILD_DIR" "$PARALLEL_BUILD_DIR"
