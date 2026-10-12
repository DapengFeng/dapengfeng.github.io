import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {fileURLToPath} from 'node:url';
import {load} from 'cheerio';

// These are the published experiment's fixtures, not performance requirements.
// CI reads evidence only: it never compiles or runs a measured implementation.
const root = fileURLToPath(new URL('../', import.meta.url));
const sourceRoot = path.join(root, 'benchmarks/parallel-libraries');
const assets = path.join(root, 'content/assets/parallel-libraries');
const implementations = [
  ['cpp', 'serial', false], ['cpp', 'std-thread', true],
  ['cpp', 'omp-static', true], ['cpp', 'omp-dynamic', true], ['cpp', 'tbb', true],
  ['rust', 'serial', false], ['rust', 'std-scoped', true], ['rust', 'rayon', true],
  ['python', 'serial', false], ['python', 'thread', true],
  ['python', 'process-copy', true], ['python', 'process-shared', true],
  ['python', 'numpy-serial', false], ['python', 'numpy-thread', true],
  ['python', 'numba-serial', false], ['python', 'numba-parallel', true],
];
const workloadNames = ['dense', 'sparse-uniform', 'sparse-skew'];
const reportCache = new Map();

test('compact article tables preserve the published measurements and their units', async () => {
  const $ = load(await fs.readFile(path.join(root, 'content/posts/parallel-libraries-dense-sparse.html'), 'utf8'));
  const tables = $('.pb-mobile-results').toArray();
  assert.equal(tables.length, 4);
  for (const table of tables) {
    const profile = $(table).attr('data-measurement-profile');
    const workload = $(table).attr('data-measurement-workload');
    const report = await reportFor(profile);
    const byId = new Map(report.summary.map(row => [row.id, row]));
    const rows = $(table).find('tbody tr').toArray();
    assert.equal(rows.length, profile === 'baseline' ? 16 : 11);
    for (const row of rows) {
      const id = $(row).attr('data-case-id'), record = byId.get(id);
      assert.ok(record, `Unknown table record ${id}`);
      const expected = workload === 'dense'
        ? [profile === 'stress' ? record.throughput_calls_per_second.median : record.warm_ms.median, record.peak_tree_pss_mib]
        : [record.warm_ms.median, byId.get(id.replace('/sparse-uniform/', '/sparse-skew/')).warm_ms.median];
      const cells = $(row).find('td data').toArray();
      assert.equal(cells.length, 2);
      cells.forEach((cell, index) => {
        close(Number($(cell).attr('value')), expected[index], `${profile}/${id}/${index}`);
        const shown = Number($(cell).text().replaceAll(',', ''));
        assert.ok(Math.abs(shown - expected[index]) <= Math.max(0.051, expected[index] * 0.00051), `${id}: visible number differs from its measurement`);
      });
    }
  }
});

async function reportFor(profile) {
  if (!reportCache.has(profile)) {
    reportCache.set(profile, fs.readFile(path.join(assets, `${profile}.json`), 'utf8').then(JSON.parse));
  }
  return reportCache.get(profile);
}

function expectedCases(profile) {
  const cases = new Map();
  for (const [language, backend, parallel] of implementations) {
    if (profile === 'stress' && !parallel) continue;
    const workers = profile === 'stress' ? [16] : parallel ? [1, 2, 4] : [1];
    for (const count of workers) for (const workload of workloadNames) {
      const size = workload === 'dense'
        ? (profile === 'baseline' ? 2048 : 4096)
        : (profile === 'baseline' ? 32768 : 131072);
      const id = `${language}:${backend}/${workload}/${count}`;
      cases.set(id, {id, language, backend, workers: count, workload, rows: size, cols: size, nnz_per_row: 64});
    }
  }
  return cases;
}

function positive(value, label) {
  assert.ok(Number.isFinite(value) && value > 0, `${label}: expected a positive finite measurement`);
}

function positiveInteger(value, label) {
  assert.ok(Number.isSafeInteger(value) && value > 0, `${label}: expected a positive integer`);
}

function close(actual, expected, label) {
  assert.ok(Number.isFinite(actual), `${label}: non-finite summary`);
  const tolerance = Math.max(1e-9, Math.abs(expected) * 1e-11);
  assert.ok(Math.abs(actual - expected) <= tolerance,
    `${label}: summary ${actual} differs from raw-derived ${expected}`);
}

function checkCase(actual, expected, label) {
  assert.ok(expected, `${label}: unexpected experiment configuration`);
  for (const [key, value] of Object.entries(expected)) assert.equal(actual[key], value, `${label}.${key}`);
}

function distribution(samples) {
  assert.ok(samples.length, 'distribution needs raw samples');
  const values = [...samples].sort((a, b) => a - b);
  const percentile = p => {
    const position = (values.length - 1) * p;
    const lower = Math.floor(position), upper = Math.ceil(position);
    return values[lower] + (values[upper] - values[lower]) * (position - lower);
  };
  return {
    count: values.length, min: values[0], max: values.at(-1),
    median: percentile(0.5), q1: percentile(0.25), q3: percentile(0.75),
    p95: values.length >= 20 ? percentile(0.95) : null,
    p99: values.length >= 100 ? percentile(0.99) : null,
  };
}

function checkDistribution(actual, samples, label) {
  assert.ok(actual && typeof actual === 'object', `${label}: missing distribution`);
  for (const [key, value] of Object.entries(distribution(samples))) {
    if (value === null || key === 'count') assert.equal(actual[key], value, `${label}.${key}`);
    else close(actual[key], value, `${label}.${key}`);
  }
}

for (const profile of ['baseline', 'stress']) {
  test(`parallel benchmark ${profile} contains every promised case and raw pass exactly once`, async () => {
    const report = await reportFor(profile), expected = expectedCases(profile), env = report.environment;
    assert.equal(report.schema, 1);
    assert.match(report.signature, /^[a-f0-9]{64}$/);
    assert.equal(env.profile, profile);
    positiveInteger(env.fresh_processes, 'fresh_processes');
    positiveInteger(env.warm_samples_per_process, 'warm_samples_per_process');
    assert.equal(env.memory_passes_per_case, 1);
    assert.equal(env.timing_memory_separate, true);
    assert.equal(env.cold_cache, false);
    positive(env.pss_interval_ms, 'PSS interval');
    if (profile === 'baseline') {
      assert.equal(env.fresh_processes, 3);
      assert.equal(env.warm_samples_per_process, 5);
      assert.equal(env.stress_seconds, 0);
    } else {
      assert.ok(env.stress_seconds >= 10, 'published sustained experiment must last at least 10 seconds per warm loop');
    }
    assert.equal(expected.size, profile === 'baseline' ? 114 : 33);
    assert.deepEqual(Object.keys(report.references).sort(), [...workloadNames].sort());

    const expectedPasses = new Set();
    for (const id of expected.keys()) {
      for (let launch = 0; launch < env.fresh_processes; launch++) expectedPasses.add(`${id}|timing|${launch}`);
      expectedPasses.add(`${id}|memory|0`);
    }
    const seenPasses = new Set();
    for (const record of report.records) {
      checkCase(record.case, expected.get(record.case.id), record.case.id);
      const key = `${record.case.id}|${record.pass}|${record.launch}`;
      assert.ok(expectedPasses.has(key), `unexpected case/pass/launch: ${key}`);
      assert.ok(!seenPasses.has(key), `duplicate case/pass/launch: ${key}`);
      seenPasses.add(key);
    }
    assert.deepEqual([...seenPasses].sort(), [...expectedPasses].sort(), 'missing raw experiment passes');
    const seenSummaries = new Set();
    for (const summary of report.summary) {
      checkCase(summary, expected.get(summary.id), summary.id);
      assert.ok(!seenSummaries.has(summary.id), `duplicate summary ${summary.id}`);
      seenSummaries.add(summary.id);
    }
    assert.deepEqual([...seenSummaries].sort(), [...expected.keys()].sort(), 'missing case summaries');
  });

  test(`parallel benchmark ${profile} validates every raw output and separates timing from sampled memory`, async () => {
    const report = await reportFor(profile);
    for (const record of report.records) {
      const label = `${record.case.id}/${record.pass}/${record.launch}`;
      const measurement = record.measurement, reference = report.references[record.case.workload];
      for (const key of ['language', 'backend', 'workers', 'workload', 'rows', 'cols', 'nnz_per_row']) {
        assert.equal(measurement[key], record.case[key], `${label}.${key}`);
      }
      assert.match(reference.checksum, /^-?\d+$/);
      assert.match(reference.output_hash, /^[a-f0-9]{16}$/);
      for (const key of ['checksum', 'output_hash', 'logical_data_bytes']) {
        assert.equal(String(measurement[key]), String(reference[key]), `${label}.${key} differs from integer oracle`);
      }
      assert.equal(measurement.validated, true, `${label}: missing validation`);
      const {rows, cols, workload, nnz_per_row} = record.case;
      const storedEntries = rows * (workload === 'dense' ? cols : nnz_per_row);
      const logicalBytes = 8 * (storedEntries + cols + rows)
        + (workload === 'dense' ? 0 : 4 * storedEntries + 8 * (rows + 1));
      assert.equal(reference.stored_entries, storedEntries, `${label}: stored matrix size`);
      assert.equal(measurement.nnz, storedEntries, `${label}: stored entries`);
      assert.equal(reference.logical_data_bytes, logicalBytes, `${label}: protocol storage size`);

      const minimum = record.pass === 'timing' ? report.environment.warm_samples_per_process : 2;
      assert.equal(measurement.repeats, minimum, `${label}: requested minimum repeats`);
      assert.equal(measurement.iterations, measurement.samples_ms.length, `${label}: iteration count`);
      assert.ok(measurement.samples_ms.length >= minimum, `${label}: missing iterations`);
      measurement.samples_ms.forEach(value => positive(value, `${label}: raw kernel sample`));
      positive(measurement.first_ms, `${label}: first call`);
      positive(measurement.setup_ms, `${label}: setup`);
      positive(measurement.elapsed_warm_ms, `${label}: warm-loop elapsed`);
      positive(record.process_wall_ms, `${label}: process wall`);
      assert.equal(measurement.stress_seconds_requested, report.environment.stress_seconds);
      const validationScope = measurement.validation_scope ?? measurement.validations;
      if (profile === 'stress') {
        assert.ok(measurement.elapsed_warm_ms >= report.environment.stress_seconds * 1000,
          `${label}: sustained warm loop ended before its requested duration`);
        assert.equal(validationScope, 'first,last', `${label}: stress validation boundary`);
      } else {
        assert.equal(measurement.samples_ms.length, minimum, `${label}: unexpected baseline iterations`);
        assert.equal(validationScope, 'all', `${label}: baseline validation boundary`);
      }
      if (record.pass === 'memory') {
        positive(record.sampled_peak_tree_pss_bytes, `${label}: sampled tree PSS`);
        positiveInteger(record.pss_samples, `${label}: PSS samples`);
        positiveInteger(record.max_observed_processes, `${label}: observed processes`);
      } else {
        for (const field of ['sampled_peak_tree_pss_bytes', 'pss_samples', 'max_observed_processes']) {
          assert.equal(record[field], null, `${label}: timing pass must not include PSS polling`);
        }
      }
    }
  });

  test(`parallel benchmark ${profile} summaries are derived from the published raw measurements`, async () => {
    const report = await reportFor(profile);
    for (const summary of report.summary) {
      const records = report.records.filter(record => record.case.id === summary.id);
      const timing = records.filter(record => record.pass === 'timing');
      const memory = records.filter(record => record.pass === 'memory');
      assert.equal(memory.length, 1, `${summary.id}: expected one separate memory pass`);
      checkDistribution(summary.warm_ms, timing.flatMap(record => record.measurement.samples_ms), `${summary.id}.warm_ms`);
      for (const field of ['first_ms', 'setup_ms']) {
        checkDistribution(summary[field], timing.map(record => record.measurement[field]), `${summary.id}.${field}`);
      }
      checkDistribution(summary.process_wall_ms, timing.map(record => record.process_wall_ms), `${summary.id}.process_wall_ms`);
      // This is loop throughput. Baseline elapsed includes validation between
      // samples; only the sustained profile excludes those intervening scans.
      checkDistribution(summary.throughput_calls_per_second, timing.map(record =>
        record.measurement.samples_ms.length / (record.measurement.elapsed_warm_ms / 1000)), `${summary.id}.loop_throughput`);
      close(summary.peak_tree_pss_mib, memory[0].sampled_peak_tree_pss_bytes / 2**20, `${summary.id}.sampled_PSS`);
      assert.equal(summary.max_processes, memory[0].max_observed_processes, `${summary.id}: process count`);
      assert.deepEqual(summary.measurement_details, timing[0].measurement, `${summary.id}: detail provenance`);
    }
  });
}

test('parallel benchmark evidence identifies the measured sources currently shipped in the repository', async () => {
  const baseline = await reportFor('baseline'), stress = await reportFor('stress');
  assert.deepEqual(stress.environment.sources_sha256, baseline.environment.sources_sha256,
    'baseline and sustained measurements must use the same frozen source files');
  const hashes = baseline.environment.sources_sha256;
  for (const required of ['run.py', 'reference.py', 'python/bench.py', 'python/requirements.txt',
    'native/parallel.cpp', 'native/build.sh', 'native/rust/Cargo.toml',
    'native/rust/Cargo.lock', 'native/rust/src/main.rs']) {
    assert.ok(Object.hasOwn(hashes, required), `missing measured source hash: ${required}`);
  }
  for (const [relative, expected] of Object.entries(hashes)) {
    assert.match(expected, /^[a-f0-9]{64}$/, relative);
    const absolute = path.resolve(sourceRoot, relative);
    assert.ok(absolute.startsWith(`${sourceRoot}${path.sep}`), `unsafe measured source path: ${relative}`);
    const actual = createHash('sha256').update(await fs.readFile(absolute)).digest('hex');
    assert.equal(actual, expected, `${relative}: source changed after measurement; republish matching evidence`);
  }
  for (const report of [baseline, stress]) {
    for (const name of ['cpp', 'rust', 'python', 'python_script']) {
      assert.match(report.environment.executables[name].sha256, /^[a-f0-9]{64}$/, `${name} binary fingerprint`);
    }
    assert.equal(report.environment.executables.python_script.sha256, hashes['python/bench.py']);
  }
});

test('parallel benchmark article exposes both raw reports and the reproducible source archive', async () => {
  const article = await fs.readFile(path.join(root, 'content/posts/parallel-libraries-dense-sparse.html'), 'utf8');
  for (const name of ['baseline.json', 'stress.json', 'source.zip']) {
    assert.match(article, new RegExp(`href=["'][^"']*(?:/|\\b)parallel-libraries/${name.replace('.', '\\.')}["']`),
      `article must link the ${name} download`);
    const info = await fs.stat(path.join(assets, name));
    assert.ok(info.isFile() && info.size > 0, `empty or missing evidence asset: ${name}`);
  }
  const zip = await fs.readFile(path.join(assets, 'source.zip'));
  assert.deepEqual([...zip.subarray(0, 4)], [0x50, 0x4b, 0x03, 0x04], 'source archive must be a ZIP file');
});
