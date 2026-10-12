// SPDX-License-Identifier: MIT
use rayon::prelude::*;
use std::{env, process, thread, time::{Duration, Instant}};

struct Options {
    backend: String, workload: String,
    rows: usize, cols: usize, nnz_per_row: usize, workers: usize, repeats: usize, hold_ms: u64,
    seconds: f64,
}
fn parse() -> Result<Options, String> {
    let mut o = Options { backend: "serial".into(), workload: "dense".into(), rows: 2048, cols: 2048, nnz_per_row: 64, workers: 1, repeats: 5, hold_ms: 0, seconds: 0.0 };
    let mut args = env::args().skip(1);
    while let Some(key) = args.next() {
        let value = args.next().ok_or("each option requires a value")?;
        let integer = || value.parse::<usize>().map_err(|_| "expected a non-negative integer".to_string());
        match key.as_str() {
            "--backend" => o.backend = value,
            "--workload" => o.workload = value,
            "--rows" => o.rows = integer()?,
            "--cols" => o.cols = integer()?,
            "--nnz-per-row" => o.nnz_per_row = integer()?,
            "--workers" => o.workers = integer()?,
            "--repeats" => o.repeats = integer()?,
            "--hold-ms" => o.hold_ms = integer()? as u64,
            "--seconds" => o.seconds = value.parse::<f64>().map_err(|_| "seconds must be finite and non-negative")?,
            _ => return Err(format!("unknown option: {key}")),
        }
    }
    if !o.seconds.is_finite() || o.seconds < 0.0 { return Err("seconds must be finite and non-negative".into()); }
    if o.rows == 0 || o.cols == 0 || o.nnz_per_row == 0 || o.workers == 0 || o.workers > 4096 || o.repeats == 0 || o.cols > u32::MAX as usize {
        return Err("positive dimensions/workers/repeats required; workers <= 4096; cols <= uint32 max".into());
    }
    if !["dense", "sparse-uniform", "sparse-skew"].contains(&o.workload.as_str()) { return Err("unknown workload".into()); }
    if o.rows > 100000000 || o.cols > 100000000 || (o.workload == "dense" && o.rows > 1000000000 / o.cols) { return Err("dimensions exceed this harness's supported bounds".into()); }
    if o.workload == "sparse-skew" && (o.rows % 16 != 0 || o.nnz_per_row != 64 || o.cols < 784) { return Err("sparse-skew requires rows divisible by 16, nnz-per-row 64, cols >= 784".into()); }
    if o.workload != "dense" && (o.nnz_per_row > o.cols || o.cols % 17 == 0) { return Err("CSR needs nnz-per-row <= cols and cols not divisible by 17 for unique columns".into()); }
    if !["serial", "std-scoped", "rayon"].contains(&o.backend.as_str()) { return Err("unknown backend".into()); }
    if o.backend == "serial" && o.workers != 1 { return Err("serial requires --workers 1".into()); }
    Ok(o)
}
struct Data { values: Vec<f64>, x: Vec<f64>, indices: Vec<u32>, indptr: Vec<u64> }
fn value(row: usize, col: usize) -> f64 { (((row * 17 + col * 13) % 101) as i64 - 50) as f64 / 64.0 }
fn make_data(o: &Options) -> Data {
    let x = (0..o.cols).map(|col| ((col % 31) as i64 - 15) as f64 / 16.0).collect();
    let mut d = Data { values: Vec::new(), x, indices: Vec::new(), indptr: Vec::new() };
    if o.workload == "dense" {
        d.values = vec![0.0; o.rows * o.cols];
        for row in 0..o.rows { for col in 0..o.cols { d.values[row * o.cols + col] = value(row, col); } }
    } else {
        d.indptr = vec![0; o.rows + 1];
        for row in 0..o.rows {
            let count = if o.workload == "sparse-skew" { if row < o.rows / 16 { 784 } else { 16 } } else { o.nnz_per_row };
            d.indptr[row + 1] = d.indptr[row] + count as u64;
        }
        let nnz = d.indptr[o.rows] as usize;
        d.values = vec![0.0; nnz]; d.indices = vec![0; nnz];
        for row in 0..o.rows {
            for k in d.indptr[row] as usize..d.indptr[row + 1] as usize {
                let col = (row * 13 + (k - d.indptr[row] as usize) * 17) % o.cols;
                d.indices[k] = col as u32; d.values[k] = value(row, col);
            }
        }
    }
    d
}
fn range<const DENSE: bool>(o: &Options, d: &Data, row_start: usize, output: &mut [f64]) {
    for (local, out) in output.iter_mut().enumerate() {
        let row = row_start + local;
        let mut sum = 0.0;
        if DENSE {
            let values = &d.values[row * o.cols..(row + 1) * o.cols];
            for (&v, &x) in values.iter().zip(d.x.iter()) { sum += v * x; }
        } else {
            for k in d.indptr[row] as usize..d.indptr[row + 1] as usize { sum += d.values[k] * d.x[d.indices[k] as usize]; }
        }
        *out = sum;
    }
}
fn run<const DENSE: bool>(o: &Options, d: &Data, output: &mut [f64], pool: Option<&rayon::ThreadPool>) {
    match o.backend.as_str() {
        "serial" => range::<DENSE>(o, d, 0, output),
        "std-scoped" => thread::scope(|scope| {
            let mut remaining = output;
            let mut start = 0;
            for worker in 0..o.workers {
                let end = o.rows * (worker + 1) / o.workers;
                let (chunk, rest) = remaining.split_at_mut(end - start); remaining = rest;
                let row_start = start;
                scope.spawn(move || range::<DENSE>(o, d, row_start, chunk));
                start = end;
            }
        }),
        "rayon" => pool.expect("Rayon pool").install(|| {
            output.par_chunks_mut(32).enumerate().for_each(|(index, out)| range::<DENSE>(o, d, index * 32, out));
        }),
        _ => unreachable!(),
    }
}
fn signature(output: &[f64]) -> Result<(i64, u64), String> {
    let (mut sum, mut hash) = (0i64, 14695981039346656037u64);
    for &y in output {
        let scaled = y * 1024.0;
        if !scaled.is_finite() || scaled != scaled.trunc() || scaled.abs() > 9007199254740991.0 { return Err("output is not an exact supported integer multiple of 1/1024".into()); }
        let q = scaled as i64;
        sum = sum.checked_add(q).ok_or("checksum overflow")?;
        hash = (hash ^ q as u64).wrapping_mul(1099511628211);
    }
    Ok((sum, hash))
}
fn main_result() -> Result<(), String> {
    let o = parse()?;
    let setup = Instant::now();
    let d = make_data(&o);
    let mut output = vec![0.0; o.rows];
    let pool = if o.backend == "rayon" { Some(rayon::ThreadPoolBuilder::new().num_threads(o.workers).build().map_err(|e| e.to_string())?) } else { None };
    let setup_ms = setup.elapsed().as_secs_f64() * 1000.0;
    let invoke = |output: &mut [f64]| if o.workload == "dense" { run::<true>(&o, &d, output, pool.as_ref()) } else { run::<false>(&o, &d, output, pool.as_ref()) };
    let first = Instant::now(); invoke(&mut output);
    let first_ms = first.elapsed().as_secs_f64() * 1000.0;
    let expected = signature(&output)?;
    let mut samples = Vec::with_capacity(o.repeats);
    let warm = Instant::now();
    while samples.len() < o.repeats || warm.elapsed().as_secs_f64() < o.seconds {
        let start = Instant::now(); invoke(&mut output);
        samples.push(start.elapsed().as_secs_f64() * 1000.0);
        if o.seconds == 0.0 && signature(&output)? != expected { return Err("output differs between repetitions".into()); }
    }
    let elapsed_warm_ms = warm.elapsed().as_secs_f64() * 1000.0;
    if o.seconds > 0.0 && signature(&output)? != expected { return Err("final stress output differs from first".into()); }
    if o.hold_ms > 0 { thread::sleep(Duration::from_millis(o.hold_ms)); }
    let logical = 8 * (d.values.len() + d.x.len() + output.len() + d.indptr.len()) + 4 * d.indices.len();
    let samples_json = samples.iter().map(|x| format!("{x:.6}")).collect::<Vec<_>>().join(",");
    println!("{{\"language\":\"rust\",\"backend\":\"{}\",\"workers\":{},\"workload\":\"{}\",\"rows\":{},\"cols\":{},\"nnz_per_row\":{},\"nnz\":{},\"repeats\":{},\"setup_ms\":{:.6},\"first_ms\":{:.6},\"stress_seconds_requested\":{},\"elapsed_warm_ms\":{:.6},\"iterations\":{},\"validations\":\"{}\",\"samples_ms\":[{}],\"checksum\":\"{}\",\"output_hash\":\"{:016x}\",\"validated\":true,\"logical_data_bytes\":{},\"rayon_version\":\"1.12.0\"}}", o.backend, o.workers, o.workload, o.rows, o.cols, o.nnz_per_row, d.values.len(), o.repeats, setup_ms, first_ms, o.seconds, elapsed_warm_ms, samples.len(), if o.seconds > 0.0 { "first,last" } else { "all" }, samples_json, expected.0, expected.1, logical);
    Ok(())
}
fn main() { if let Err(error) = main_result() { eprintln!("benchmark error: {error}"); process::exit(1); } }
