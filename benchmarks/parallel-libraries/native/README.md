# Native matvec runners / 原生矩阵向量乘法程序

This directory supplies the C++ and Rust runners for the parent benchmark. They implement the same row-wise dense and CSR matrix–vector product with `float64` values, a `uint32` CSR column index and a `uint64` row pointer. They are standalone programs, not website JavaScript.

这里是基准测试的 C++ 与 Rust 程序，按行执行相同的稠密矩阵／CSR 矩阵与向量乘法。数值为 `float64`，CSR 列下标为 `uint32`，行指针为 `uint64`。它们是独立程序，不是网页 JavaScript。

## Build / 构建

Requirements: a C++20 compiler with OpenMP, oneTBB headers/libraries, and Rust with Cargo. The measured versions and flags belong in the result manifest. `Cargo.lock` pins Rayon and its transitive dependencies. `build.sh` uses `-O3 -march=native -fno-fast-math -ffp-contract=off` for C++, and Rust release mode with `-C target-cpu=native`; neither enables reassociated floating-point reductions.

需要支持 OpenMP 的 C++20 编译器、oneTBB 头文件与库、Rust 和 Cargo。实际测量版本与选项记录在结果清单中。`Cargo.lock` 锁定 Rayon 及其间接依赖。`build.sh` 的 C++ 选项为 `-O3 -march=native -fno-fast-math -ffp-contract=off`，Rust 使用 release 模式和 `-C target-cpu=native`；两者均不启用浮点归约重排。

The local build used GCC 13.3.0 (Ubuntu `13.3.0-6ubuntu2~24.04.1`), oneTBB 2021.11.0 (`2021.11.0-2ubuntu2`), Rust 1.91.1 (`ed61e7d7e`, LLVM 21.1.2), Cargo 1.91.1 and Rayon 1.12.0 / rayon-core 1.13.0. Rust components were installed only in a temporary directory; mirrored `rustc` and `rust-std` archives were checked against the SHA-256 values in the official Rust 1.91.1 release manifest before installation. No global toolchain or system package was changed.

本地构建使用 GCC 13.3.0（Ubuntu `13.3.0-6ubuntu2~24.04.1`）、oneTBB 2021.11.0（`2021.11.0-2ubuntu2`）、Rust 1.91.1（`ed61e7d7e`、LLVM 21.1.2）、Cargo 1.91.1 和 Rayon 1.12.0／rayon-core 1.13.0。Rust 组件仅安装于临时目录；镜像下载的 `rustc`、`rust-std` 压缩包在安装前均按 Rust 1.91.1 官方发布清单的 SHA-256 核对。未改动全局工具链或系统软件包。

From the repository root / 从仓库根目录执行：

```bash
bash benchmarks/parallel-libraries/native/build.sh

# Optional: extracted oneTBB packages, without a system installation.
# 可选：使用解压后的 oneTBB 软件包，无需安装到系统。
TBB_ROOT=/tmp/parallel-native-deps/tbb/usr \
  PARALLEL_BUILD_DIR=/tmp/parallel-native-build \
  bash benchmarks/parallel-libraries/native/build.sh

# Without oneTBB, omit that backend explicitly / 没有 oneTBB 时显式省略该后端
ENABLE_TBB=0 bash benchmarks/parallel-libraries/native/build.sh
```

`TBB_ROOT` contains `include/` and `lib/x86_64-linux-gnu/`; set `TBB_LIB_DIR` if your library path differs. `PARALLEL_BUILD_DIR` defaults to ignored `native/.build/`. `CXX` and `RUSTFLAGS` can override the compiler/CPU target; record any override when publishing measurements. Run `cargo build --offline --locked` separately if all pinned crates are already cached and network access is unnecessary.

`TBB_ROOT` 下需有 `include/` 和 `lib/x86_64-linux-gnu/`；库目录不同时设置 `TBB_LIB_DIR`。`PARALLEL_BUILD_DIR` 默认为已忽略的 `native/.build/`。`CXX`、`RUSTFLAGS` 可以更换编译器或 CPU 目标，发表测量结果时记录这些变更。若已缓存所有锁定的 crate，可单独用 `cargo build --offline --locked` 避免联网。

## CLI and measurement boundary / 命令与测量边界

```bash
benchmarks/parallel-libraries/native/.build/parallel-cpp \
  --backend omp-dynamic --workers 4 --workload sparse-skew \
  --rows 32768 --cols 32768 --nnz-per-row 64 --repeats 5

benchmarks/parallel-libraries/native/.build/parallel-rust \
  --backend rayon --workers 4 --workload dense \
  --rows 2048 --cols 2048 --repeats 5
```

| Runtime / 运行方式 | Backend / 名称 | Scheduling / 调度 |
| --- | --- | --- |
| C++ | `serial` | Calling thread; requires one worker.<br>调用线程，要求一个 worker。 |
| C++ | `std-thread` | Balanced contiguous row ranges; spawn and join on every call.<br>连续行均分，每次调用创建并 join 线程。 |
| C++ | `omp-static` | OpenMP `schedule(static)`.<br>OpenMP 静态分配。 |
| C++ | `omp-dynamic` | OpenMP `schedule(dynamic, 32)`.<br>OpenMP 按 32 行动态领取。 |
| C++ | `tbb` | `parallel_for`, minimum range 32 rows, `auto_partitioner`.<br>最小范围 32 行，使用自动分区器。 |
| Rust | `serial` | Calling thread; requires one worker.<br>调用线程，要求一个 worker。 |
| Rust | `std-scoped` | Balanced contiguous row ranges; scoped threads spawn/join per call.<br>连续行均分，每次调用创建／回收 scoped 线程。 |
| Rust | `rayon` | A persistent custom pool, parallel chunks of 32 output rows.<br>持久的专用线程池，输出按 32 行分块并行。 |

`setup_ms` includes input construction, output allocation and explicit pool construction. Rayon constructs its pool during setup. OpenMP and oneTBB can start workers lazily during `first_ms`; no artificial warm-up precedes that first call. The repeated `samples_ms` include scheduling, matvec and waiting for completion. Standard-library backends include their per-call thread creation/join. Validation and JSON serialization are outside all kernel timings. This measures these concrete implementations and lifecycles, not an isolated universal “library overhead”.

`setup_ms` 包含输入构造、输出分配和显式线程池构造。Rayon 在线程池构造阶段创建 worker；OpenMP 和 oneTBB 可在 `first_ms` 内延迟启动 worker，首次调用前不做额外预热。后续 `samples_ms` 包含调度、矩阵向量乘法和等待完成，标准库后端还包含每次调用的线程创建与 join。验证和 JSON 序列化位于所有内核计时之外。结果针对这些实现与生命周期，不等于可以独立比较的普遍“库开销”。

The input vector is `(col % 31 - 15) / 16`; matrix values are `((row*17 + col*13) % 101 - 50) / 64`, with signed subtraction. CSR column `k` is `(row*13 + k*17) % cols`. Uniform CSR has `nnz-per-row` stored entries per row. Skewed CSR requires an average of 64: its first 1/16 rows have 784 entries each, and all remaining rows have 16. Some stored values are zero; “nnz” in the JSON is the stored-entry count. Sparse column counts must not be divisible by 17 so the generated columns are unique within a row.

输入向量为 `(col % 31 - 15) / 16`；矩阵数值为 `((row*17 + col*13) % 101 - 50) / 64`，减法按有符号数执行。CSR 每行第 `k` 列为 `(row*13 + k*17) % cols`。均匀 CSR 每行存储 `nnz-per-row` 项；倾斜 CSR 固定平均 64 项：前 1/16 行各 784 项，其余行各 16 项。一些存储数值为零，JSON 中的 `nnz` 指存储项数。稀疏列数不能被 17 整除，以保证每行生成的列下标不重复。

Every output is exactly representable in units of 1/1024 at the stated sizes. After every call, the runner checks the signed integer sum of these units and a value-wise FNV-1a-style 64-bit hash. The hash starts at 14695981039346656037 and, for each signed row result `q`, applies `(hash XOR uint64(q)) * 1099511628211` with wrapping multiplication. This is a value-wise checksum protocol, not byte-wise standard FNV-1a. `validated: true` means repeated outputs agree; the parent runner must additionally compare them to its independent integer reference.

所用尺寸下，每行输出都能精确表示为 1/1024 的整数倍。每次调用后检查这些整数的有符号总和以及逐值的 FNV-1a 风格 64 位哈希。哈希初值为 14695981039346656037，对每行有符号整数 `q` 执行 `(hash XOR uint64(q)) * 1099511628211`，乘法溢出按 64 位回绕。这是逐值协议，不是标准 FNV-1a 的逐字节协议。`validated: true` 仅表示各次结果一致；父目录的运行器还需与独立整数参考值核对。

`logical_data_bytes` counts matrix values, CSR indices/pointers where present, input vector and output vector. Linux C++ `peak_rss_kib` comes from `getrusage`; it does not include child processes and is not PSS. Use the parent runner for process-tree memory measurement. `--hold-ms 300` holds the live buffers/pools after all samples solely for a separate memory pass; leave it at zero for timing runs.

`logical_data_bytes` 统计矩阵数值、CSR 下标／行指针（若存在）、输入向量和输出向量。Linux C++ 的 `peak_rss_kib` 来自 `getrusage`，不含子进程，也不是 PSS。进程树内存由父目录运行器统计。`--hold-ms 300` 会在所有采样后暂时保留缓冲区和线程池，仅用于独立内存测量；计时轮次保持默认零值。

For sustained load, `--seconds 10` repeats the warm kernel until both the requested duration and the minimum `--repeats` count are met. It validates the first and final outputs, rather than running a full hash inside every stress iteration. `elapsed_warm_ms` covers this warm loop; `iterations` and every kernel sample remain in the JSON. With default `--seconds 0`, every iteration is validated outside its individual kernel timer. The parent runner sets CPU affinity, records the environment and keeps timing and memory passes separate.

持续负载可指定 `--seconds 10`：预热后的内核会重复执行，直到同时达到时长和 `--repeats` 的最低次数。此时只验证首次与最终输出，不在每次压测迭代之间计算完整哈希。`elapsed_warm_ms` 记录这段循环的总时长；JSON 仍保存 `iterations` 与每次内核耗时。默认 `--seconds 0` 时，每次输出都在对应内核计时之外验证。CPU 亲和性、环境记录及计时／内存分轮由父目录运行器管理。
