// SPDX-License-Identifier: MIT
// Row-wise float64 dense and CSR matvec, matching ../python/ and rust/.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <sys/resource.h>
#include <omp.h>
#ifdef BENCH_WITH_TBB
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/partitioner.h>
#include <oneapi/tbb/version.h>
#endif

using Clock = std::chrono::steady_clock;
struct Options {
  std::string backend = "serial", workload = "dense";
  std::size_t rows = 2048, cols = 2048, nnz_per_row = 64;
  unsigned workers = 1, repeats = 5, hold_ms = 0;
  double seconds = 0;
};
static double ms(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}
static std::uint64_t number(const char* value) {
  if (!*value || *value == '-') throw std::runtime_error("expected a non-negative integer");
  std::size_t consumed = 0;
  auto result = std::stoull(value, &consumed);
  if (value[consumed] != '\0') throw std::runtime_error("invalid integer");
  return result;
}
static Options parse(int argc, char** argv) {
  Options o;
  for (int i = 1; i < argc; i += 2) {
    if (i + 1 == argc) throw std::runtime_error("each option requires a value");
    const std::string key = argv[i];
    if (key == "--backend") o.backend = argv[i + 1];
    else if (key == "--workload") o.workload = argv[i + 1];
    else if (key == "--seconds") {
      std::size_t consumed = 0;
      o.seconds = std::stod(argv[i + 1], &consumed);
      if (argv[i + 1][consumed] != '\0' || !std::isfinite(o.seconds) || o.seconds < 0)
        throw std::runtime_error("seconds must be finite and non-negative");
    }
    else {
      const auto n = number(argv[i + 1]);
      if (key == "--rows") o.rows = n;
      else if (key == "--cols") o.cols = n;
      else if (key == "--nnz-per-row") o.nnz_per_row = n;
      else {
        if (n > std::numeric_limits<unsigned>::max()) throw std::runtime_error("option too large");
        if (key == "--workers") o.workers = n;
        else if (key == "--repeats") o.repeats = n;
        else if (key == "--hold-ms") o.hold_ms = n;
        else throw std::runtime_error("unknown option: " + key);
      }
    }
  }
  if (!o.rows || !o.cols || !o.nnz_per_row || !o.workers || !o.repeats || o.workers > 4096 || o.cols > UINT32_MAX)
    throw std::runtime_error("positive dimensions/workers/repeats required; workers <= 4096; cols <= uint32 max");
  if (o.workload != "dense" && o.workload != "sparse-uniform" && o.workload != "sparse-skew") throw std::runtime_error("unknown workload");
  if (o.rows > 100000000 || o.cols > 100000000 || (o.workload == "dense" && o.rows > 1000000000 / o.cols))
    throw std::runtime_error("dimensions exceed this harness's supported bounds");
  if (o.workload == "sparse-skew" && (o.rows % 16 || o.nnz_per_row != 64 || o.cols < 784))
    throw std::runtime_error("sparse-skew requires rows divisible by 16, nnz-per-row 64, cols >= 784");
  if (o.workload != "dense" && (o.nnz_per_row > o.cols || o.cols % 17 == 0))
    throw std::runtime_error("CSR needs nnz-per-row <= cols and cols not divisible by 17 for unique columns");
  if (o.backend != "serial" && o.backend != "std-thread" && o.backend != "omp-static" && o.backend != "omp-dynamic" && o.backend != "tbb")
    throw std::runtime_error("unknown backend");
  if (o.backend == "serial" && o.workers != 1) throw std::runtime_error("serial requires --workers 1");
#ifndef BENCH_WITH_TBB
  if (o.backend == "tbb") throw std::runtime_error("build with BENCH_WITH_TBB to enable oneTBB");
#endif
  return o;
}
struct Data {
  std::vector<double> values, x, output;
  std::vector<std::uint32_t> indices;
  std::vector<std::uint64_t> indptr;
  std::uint64_t logical_bytes() const {
    return 8 * (values.size() + x.size() + output.size() + indptr.size()) + 4 * indices.size();
  }
};
static double value(std::size_t row, std::size_t col) {
  return (static_cast<std::int64_t>((row * 17 + col * 13) % 101) - 50) / 64.0;
}
static Data make_data(const Options& o) {
  Data d;
  d.x.resize(o.cols);
  d.output.resize(o.rows);
  for (std::size_t col = 0; col < o.cols; ++col) d.x[col] = (static_cast<std::int64_t>(col % 31) - 15) / 16.0;
  if (o.workload == "dense") {
    d.values.resize(o.rows * o.cols);
    for (std::size_t row = 0; row < o.rows; ++row)
      for (std::size_t col = 0; col < o.cols; ++col) d.values[row * o.cols + col] = value(row, col);
  } else {
    d.indptr.resize(o.rows + 1);
    for (std::size_t row = 0; row < o.rows; ++row)
      d.indptr[row + 1] = d.indptr[row] + (o.workload == "sparse-skew" ? (row < o.rows / 16 ? 784 : 16) : o.nnz_per_row);
    d.values.resize(d.indptr.back());
    d.indices.resize(d.indptr.back());
    for (std::size_t row = 0; row < o.rows; ++row)
      for (auto k = d.indptr[row]; k < d.indptr[row + 1]; ++k) {
        const auto col = (row * 13 + (k - d.indptr[row]) * 17) % o.cols;
        d.indices[k] = static_cast<std::uint32_t>(col);
        d.values[k] = value(row, col);
      }
  }
  return d;
}
template<bool Dense>
static void range(const Options& o, Data& d, std::size_t begin, std::size_t end) {
  for (auto row = begin; row < end; ++row) {
    double sum = 0;
    if constexpr (Dense) {
      const auto offset = row * o.cols;
      for (std::size_t col = 0; col < o.cols; ++col) sum += d.values[offset + col] * d.x[col];
    } else {
      for (auto k = d.indptr[row]; k < d.indptr[row + 1]; ++k) sum += d.values[k] * d.x[d.indices[k]];
    }
    d.output[row] = sum;
  }
}
template<bool Dense>
static void run(const Options& o, Data& d) {
  if (o.backend == "serial") { range<Dense>(o, d, 0, o.rows); return; }
  if (o.backend == "std-thread") {
    std::vector<std::thread> threads;
    threads.reserve(o.workers);
    for (unsigned w = 0; w < o.workers; ++w)
      threads.emplace_back([&, w] { range<Dense>(o, d, o.rows * w / o.workers, o.rows * (w + 1) / o.workers); });
    for (auto& thread : threads) thread.join();
  } else if (o.backend == "omp-static") {
    #pragma omp parallel for schedule(static) num_threads(o.workers)
    for (std::size_t row = 0; row < o.rows; ++row) range<Dense>(o, d, row, row + 1);
  } else if (o.backend == "omp-dynamic") {
    #pragma omp parallel for schedule(dynamic, 32) num_threads(o.workers)
    for (std::size_t row = 0; row < o.rows; ++row) range<Dense>(o, d, row, row + 1);
  } else {
#ifdef BENCH_WITH_TBB
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, o.rows, 32),
      [&](const oneapi::tbb::blocked_range<std::size_t>& r) { range<Dense>(o, d, r.begin(), r.end()); }, oneapi::tbb::auto_partitioner{});
#else
    throw std::runtime_error("backend unavailable");
#endif
  }
}
struct Signature { std::int64_t checksum; std::uint64_t hash; };
static Signature signature(const Data& d) {
  Signature s{0, UINT64_C(14695981039346656037)};
  for (auto y : d.output) {
    const double scaled = y * 1024.0;
    if (!std::isfinite(scaled) || scaled != std::trunc(scaled) || std::abs(scaled) > 9007199254740991.0)
      throw std::runtime_error("output is not an exact supported integer multiple of 1/1024");
    const auto q = static_cast<std::int64_t>(scaled);
    s.checksum += q;
    s.hash = (s.hash ^ static_cast<std::uint64_t>(q)) * UINT64_C(1099511628211);
  }
  return s;
}
int main(int argc, char** argv) {
  try {
    const auto o = parse(argc, argv);
    const auto setup_start = Clock::now();
    auto d = make_data(o);
    omp_set_dynamic(0);
#ifdef BENCH_WITH_TBB
    std::unique_ptr<oneapi::tbb::global_control> tbb_limit;
    if (o.backend == "tbb") tbb_limit = std::make_unique<oneapi::tbb::global_control>(oneapi::tbb::global_control::max_allowed_parallelism, o.workers);
#endif
    const auto setup_ms = ms(setup_start);
    const auto invoke = [&] { if (o.workload == "dense") run<true>(o, d); else run<false>(o, d); };
    const auto first_start = Clock::now();
    invoke();
    const auto first_ms = ms(first_start);
    const auto expected = signature(d);
    std::vector<double> samples;
    samples.reserve(o.repeats);
    const auto warm_start = Clock::now();
    while (samples.size() < o.repeats || ms(warm_start) < o.seconds * 1000.0) {
      const auto start = Clock::now();
      invoke();
      samples.push_back(ms(start));
      if (o.seconds == 0) {
        const auto current = signature(d);
        if (current.checksum != expected.checksum || current.hash != expected.hash) throw std::runtime_error("output differs between repetitions");
      }
    }
    const auto elapsed_warm_ms = ms(warm_start);
    if (o.seconds > 0) {
      const auto current = signature(d);
      if (current.checksum != expected.checksum || current.hash != expected.hash) throw std::runtime_error("final stress output differs from first");
    }
    if (o.hold_ms) std::this_thread::sleep_for(std::chrono::milliseconds(o.hold_ms));
    struct rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    std::cout << std::fixed << std::setprecision(6)
      << "{\"language\":\"cpp\",\"backend\":\"" << o.backend << "\",\"workers\":" << o.workers
      << ",\"workload\":\"" << o.workload << "\",\"rows\":" << o.rows << ",\"cols\":" << o.cols
      << ",\"nnz_per_row\":" << o.nnz_per_row << ",\"nnz\":" << d.values.size() << ",\"repeats\":" << o.repeats
      << ",\"setup_ms\":" << setup_ms << ",\"first_ms\":" << first_ms
      << ",\"stress_seconds_requested\":" << o.seconds << ",\"elapsed_warm_ms\":" << elapsed_warm_ms
      << ",\"iterations\":" << samples.size() << ",\"validations\":\"" << (o.seconds > 0 ? "first,last" : "all") << "\",\"samples_ms\":[";
    for (std::size_t i = 0; i < samples.size(); ++i) std::cout << (i ? "," : "") << samples[i];
    std::cout << "],\"checksum\":\"" << expected.checksum << "\",\"output_hash\":\""
      << std::hex << std::setw(16) << std::setfill('0') << expected.hash << std::dec << "\",\"validated\":true"
      << ",\"logical_data_bytes\":" << d.logical_bytes() << ",\"peak_rss_kib\":" << usage.ru_maxrss
      << ",\"compiler\":\"" << __VERSION__ << "\",\"openmp\":" << _OPENMP;
#ifdef BENCH_WITH_TBB
    std::cout << ",\"tbb_version\":\"" << TBB_VERSION_STRING << "\"";
#endif
    std::cout << "}\n";
  } catch (const std::exception& error) { std::cerr << "benchmark error: " << error.what() << '\n'; return 1; }
}
