# Parallel dense and sparse matvec / 稠密与稀疏并行测试

This experiment supports the [benchmark article](../../content/posts/parallel-libraries-dense-sparse.html). It compares concrete CPU implementations of `y = A x`, including their scheduling and data-transfer costs. It does not rank programming languages independently of their compilers, numerical kernels or execution models. Dense matvec is also a low-arithmetic-intensity workload; these results do not represent dense matrix–matrix multiplication.

本实验对应[并行测试文章](../../content/posts/parallel-libraries-dense-sparse.html)，比较 CPU 上 `y = A x` 的具体实现，包含调度与数据传输成本。结果不能脱离编译器、数值内核与执行模型，解释成语言排名。稠密矩阵向量乘法本身也是低计算强度任务，不能代表稠密矩阵乘法。

## Reproduce / 复现

Use Linux with readable `/proc/*/smaps_rollup`, Python 3.12, a C++20 compiler with OpenMP, oneTBB, and Rust/Cargo. See [native build details](native/README.md) for flags, versions and individual backend commands. Install Python dependencies in an isolated environment. Pick CPU IDs from your own `lscpu -p=CPU,CORE,SOCKET` output; the IDs below describe the measured i7-1270P, not other machines.

需要 Linux、可读取的 `/proc/*/smaps_rollup`、Python 3.12、支持 OpenMP 的 C++20 编译器、oneTBB 和 Rust/Cargo。编译参数、版本与单个后端的命令见[原生实现说明](native/README.md)。Python 依赖安装在独立环境。根据自己机器的 `lscpu -p=CPU,CORE,SOCKET` 选择 CPU 编号；以下编号仅对应实测的 i7-1270P。

```bash
python3 -m venv benchmarks/parallel-libraries/.venv
benchmarks/parallel-libraries/.venv/bin/pip install -r benchmarks/parallel-libraries/python/requirements.txt
bash benchmarks/parallel-libraries/native/build.sh
python3 benchmarks/parallel-libraries/reference.py --self-test

# Correctness smoke; never publish these tiny inputs as performance results.
# 正确性冒烟检查；不将这些微小输入的耗时作为性能结论。
python3 benchmarks/parallel-libraries/run.py \
  --cpp benchmarks/parallel-libraries/native/.build/parallel-cpp \
  --rust benchmarks/parallel-libraries/native/.build/parallel-rust \
  --python benchmarks/parallel-libraries/.venv/bin/python \
  --smoke --launches 1 --repeats 2 --workers 2 \
  --cpus 0,2,4,6 --monitor-cpus 8,9,10,11,12,13,14,15 \
  --output test-results/parallel-smoke

# Baseline: three fresh processes, five warm samples each, separate PSS pass.
# 基础测试：三个全新进程，各五次热调用，PSS 单独测量。
python3 benchmarks/parallel-libraries/run.py \
  --cpp benchmarks/parallel-libraries/native/.build/parallel-cpp \
  --rust benchmarks/parallel-libraries/native/.build/parallel-rust \
  --python benchmarks/parallel-libraries/.venv/bin/python \
  --cpus 0,2,4,6 --monitor-cpus 8,9,10,11,12,13,14,15 \
  --workers 1,2,4 --output test-results/parallel-baseline

# Run after the baseline, never concurrently.
# 基础测试结束后再压测，两者不能同时运行。
python3 benchmarks/parallel-libraries/run.py \
  --cpp benchmarks/parallel-libraries/native/.build/parallel-cpp \
  --rust benchmarks/parallel-libraries/native/.build/parallel-rust \
  --python benchmarks/parallel-libraries/.venv/bin/python \
  --profile stress --seconds 10 --launches 1 --repeats 5 \
  --workers 16 --cpus 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 \
  --monitor-cpus 8,9,10,11,12,13,14,15 \
  --output test-results/parallel-stress
```

The baseline comprises 114 configurations and 456 fresh process launches; the stress profile comprises 33 configurations and 66 launches. Every case runs sequentially, with a deterministic shuffled order within each pass. `--only cpp:omp-static,rust:rayon` selects implementations. `--resume` accepts completed records only when the source, executable, runtime, inputs and measurement settings match. Keep unrelated builds and browser tests stopped while measuring.

基础测试共 114 种配置、456 次进程启动；压测共 33 种配置、66 次启动。每项串行执行，每轮按固定种子打乱顺序。`--only cpp:omp-static,rust:rayon` 可选择实现；`--resume` 仅在源码、可执行文件、运行时、输入与测量设置一致时使用已完成记录。测量期间停止无关构建与浏览器测试。

## Inputs and validation / 输入与校验

| Profile / 组别 | Dense / 稠密 | CSR uniform and skew / 均匀与倾斜 CSR |
| --- | --- | --- |
| Baseline / 基础 | 2048 × 2048 | 32768 × 32768; 2,097,152 stored entries / 存储项 |
| Stress / 压测 | 4096 × 4096 | 131072 × 131072; 8,388,608 stored entries / 存储项 |

Every implementation uses the same deterministic FP64 values and column-index rule. Uniform CSR has 64 stored entries per row. Skewed CSR puts 784 in each of the first 1/16 rows and 16 in each remaining row, preserving the total. Some stored values are zero. Column access patterns also differ between these matrices, so their performance ratio is not an isolated measurement of scheduling overhead.

各实现使用相同规则生成的 FP64 数值与列下标。均匀 CSR 每行 64 项；倾斜 CSR 的前 1/16 行各 784 项，其余行各 16 项，总项数不变。部分存储值为零。两种矩阵的列访问模式也有差异，因此耗时比不能单独归因于调度。

`reference.py` independently computes exact integer row results in units of 1/1024. The collector checks the signed sum, row-order-sensitive 64-bit hash and logical byte count against this oracle. Fixed-repeat runs validate every output outside its kernel timer; sustained runs validate the first and final outputs, avoiding validation work between stress iterations. A mismatch stops the run.

`reference.py` 独立计算以 1/1024 为单位的精确整数行结果。采集器对照此参考核对有符号总和、对行顺序敏感的 64 位哈希与逻辑字节数。固定次数测试在每次内核计时之外校验输出；持续压测只校验首次与最终结果，避免在压测迭代之间插入校验负载。不匹配立即停止。

## Measurement boundaries / 测量边界

- **Time:** preserve setup, first call, individual warm calls and complete process wall time separately. First-call Numba includes JIT; standard-library threads are created/joined per call. Python process pools use `spawn`. `process-copy` serializes inputs and returns outputs every call; `process-shared` prepares shared buffers once. NumPy/SciPy use compiled kernels; BLAS is fixed to one thread even inside an outer thread pool.

  **时间：**分别保存准备、首次调用、逐次热调用与完整进程耗时。Numba 首次调用包含 JIT；标准库线程每次创建并 join。Python 进程池使用 `spawn`；`process-copy` 每次序列化输入并返回结果，`process-shared` 一次准备共享缓冲区。NumPy/SciPy 使用编译内核，即使外层有线程池，BLAS 也固定为单线程。

- **Space:** sample the live parent and descendant process tree using Linux PSS in a separate run. Sum processes in each sweep, then keep the maximum. This includes imports, JIT, startup and a final 100 ms buffer hold. The target interval is 10 ms plus sampling work; snapshots are not atomic and brief peaks may be missed. This is a sampled lifetime footprint, not allocated bytes, a leak test or parent-only RSS. Only one memory pass is collected per case.

  **空间：**独立轮次用 Linux PSS 采样仍存活的父子进程树，每轮求和再取峰值。覆盖导入、JIT、启动与最后 100 ms 的缓冲区保留阶段。目标间隔为 10 ms 加采样工作本身；快照不具原子性，可能错过短促峰值。这是采样的生命周期内存占用，不是分配量、泄漏检测或仅主进程 RSS。每项仅测一轮内存。

- **Statistics:** baseline plots use the median and interquartile range of 15 warm calls across three processes. They are descriptive, correlated samples, not confidence intervals. P95 needs at least 20 samples; P99 needs at least 100. Sustained throughput is completed warm calls divided by warm-loop wall time. The baseline loop includes between-call validation, so its loop throughput is not used for performance conclusions.

  **统计：**基础图使用三个进程共 15 次热调用的中位数和四分位距，是存在相关性的描述性样本，不是置信区间。至少 20 个样本才报告 P95，至少 100 个才报告 P99。持续吞吐量为热循环完成次数除以循环墙钟时间；基础测试循环含调用间校验，不将其循环吞吐量用于性能结论。

- **Stress:** repeat a larger matrix for at least 10 seconds and at least five iterations, using all 16 logical CPUs of this hybrid 12-core processor. This changes matrix size, CPU composition and worker count together. It is a bounded sustained-load experiment, not a scaling ratio against the four-P-core baseline, a saturation proof, a thermal diagnosis or a multi-hour endurance test. The monitoring process shares E cores during this profile.

  **压测：**使用更大矩阵、混合架构 12 核处理器的全部 16 个逻辑 CPU，连续执行至少 10 秒且至少五次。这同时改变矩阵尺寸、CPU 组成与并行度，属于有边界的持续负载实验，不能直接与四个 P 核的基础测试组成扩展比，也不能证明资源饱和、诊断温度或替代数小时耐久测试。该组监控进程与负载共用 E 核。

The collector stops a process group after 120 seconds, when its conservatively summed RSS exceeds 4096 MiB, or when host available memory falls below 2048 MiB. The memory guard runs about once per second; it is a protective stop condition, not a hard OS allocation cap. Shared pages may be counted repeatedly by this guard; only PSS is reported as the space result. No incomplete result is promoted to a successful case.

采集器在单进程组超过 120 秒、保守累加 RSS 超过 4096 MiB 或主机可用内存低于 2048 MiB 时停止进程组。内存守卫约每秒检查一次，是保护性停止条件，不是操作系统强制分配上限；其中共享页可能重复计数，正式空间结果只采用 PSS。不把不完整结果记录为成功项。

## Publication and maintenance / 发布与维护

Raw runs stay under ignored `test-results/`. After completion, `scripts/render-parallel-benchmark.py` produces article assets from the reports, including raw JSON, CSV, SVG/PNG figures and a source ZIP. Install plotting dependencies separately from the measured Python environment; the renderer needs Matplotlib 3.10.7. Run `--help` for its arguments. Use `--context-note` to label documented run conditions in exported figures; the article’s published run uses `"Concurrent CPU workloads; cores not reserved"`. Its baseline was measured with other CPU-heavy jobs active, so it is not an idle-host ranking. The source ZIP is checked against the reports' measured source hashes. Never hand-edit a plotted value or reuse measurements after changing their source code.

原始运行保存在已忽略的 `test-results/`。完成后，`scripts/render-parallel-benchmark.py` 从报告生成文章资源，包括原始 JSON、CSV、SVG/PNG 图与源码 ZIP。绘图依赖与实测 Python 环境分开安装；绘图器需要 Matplotlib 3.10.7，具体参数见 `--help`。可用 `--context-note` 在导出图中标注有依据的环境条件；本次文章使用 `"Concurrent CPU workloads; cores not reserved"`。已发布的基础数据采集时存在其他高 CPU 负载，不能作为空闲主机排名。源码 ZIP 对照报告中的实测源码哈希检查。不手改图中数值，不在源码变化后沿用旧测量。

Benchmarks are manually invoked, never run during the website build or as timing thresholds on shared CI runners. Validate correctness and published artifacts in site checks; reproduce performance explicitly on an identified machine. Shared testing requirements remain in [the testing guide](../../docs/development/testing.md).

性能实验手动执行，不在网站构建时运行，也不以共享 CI runner 的耗时作为门槛。站点检查验证正确性与发布资源；性能复现应明确指定机器。共用验证要求见[测试指南](../../docs/development/testing.md)。
