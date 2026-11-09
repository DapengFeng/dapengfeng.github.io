import {test} from 'node:test';
import assert from 'node:assert/strict';
import {bernoulliRunningMeans, bernoulliMeanVariance, chebyshevBound, standardizedUniformTail, conditionalProbability, covarianceStatistics, gaussianDensity, gaussianMixtureDensity, uniformMaximumCDF, uniformMaximumDensity, uniformMaximumMean, triangularCDF, triangularQuantile, bernoulliKL, probabilityRenderers, drawProbabilityConstruction} from '../src/scripts/daily-math-probability.js';
const near = (a, b, tolerance = 1e-10) => assert.ok(Math.abs(a - b) <= tolerance, `${a} differs from ${b}`);
const integrate = (f, a, b, n = 4000) => { const h = (b - a) / n; let sum = f(a) + f(b); for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * f(a + h * i); return sum * h / 3; };

test('running means retain a reproducible Bernoulli sequence and mean variance matches exact enumeration', () => {
  const path = bernoulliRunningMeans(256, .5, 11);
  assert.deepEqual(path, bernoulliRunningMeans(256, .5, 11));
  assert.notDeepEqual(path, bernoulliRunningMeans(256, .5, 12));
  path.forEach((mean, i) => {assert.ok(mean >= 0 && mean <= 1); const bit = mean * (i + 1) - (i ? path[i - 1] * i : 0); near(bit, Math.round(bit)); assert.ok(Math.round(bit) === 0 || Math.round(bit) === 1);});
  assert.ok(bernoulliRunningMeans(20, 0).every(x => x === 0)); assert.ok(bernoulliRunningMeans(20, 1).every(x => x === 1));
  const n = 8, means = Array.from({length: 1 << n}, (_, state) => state.toString(2).replaceAll('0', '').length / n);
  const variance = means.reduce((sum, mean) => sum + (mean - .5) ** 2, 0) / means.length;
  near(variance, bernoulliMeanVariance(n, .5)); near(bernoulliMeanVariance(400, .3), bernoulliMeanVariance(100, .3) / 4);
  assert.throws(() => bernoulliRunningMeans(0), RangeError); assert.throws(() => bernoulliRunningMeans(5, .5, -1), RangeError);
});

test('Chebyshev bounds the exact standardized-uniform tail and is capped at probability one', () => {
  for (let i = 1; i <= 100; i++) {const k = i / 20; assert.ok(standardizedUniformTail(k) <= chebyshevBound(1, k));}
  near(chebyshevBound(1, 2), .25); near(standardizedUniformTail(2), 0); near(chebyshevBound(1, .5), 1); near(chebyshevBound(0, 2), 0);
  near(standardizedUniformTail(Math.sqrt(3) / 2), .5);
  assert.throws(() => chebyshevBound(-1, 1), RangeError); assert.throws(() => chebyshevBound(1, 0), RangeError);
});

test('conditional probability renormalizes intersections and their complements', () => {
  for (const condition of [.2, .5, 1]) for (const proportion of [0, .3, 1]) {
    const joint = condition * proportion; near(conditionalProbability(joint, condition), proportion);
    near(conditionalProbability(joint, condition) + conditionalProbability(condition - joint, condition), 1);
  }
  near(conditionalProbability(.2, .6), 1 / 3);
  assert.throws(() => conditionalProbability(.3, .2), RangeError); assert.throws(() => conditionalProbability(0, 0), RangeError);
});

test('covariance carries units while standardized correlation is invariant to positive affine rescaling', () => {
  const points = [[-3, 2], [-1, -2], [0, 1], [2, 4], [5, 3]], a = covarianceStatistics(points);
  const b = covarianceStatistics(points.map(([x, y]) => [3 * x + 7, 2 * y - 4]));
  near(b.covariance, 6 * a.covariance); near(b.correlation, a.correlation);
  const standardized = covarianceStatistics(a.standardized);
  near(standardized.mean[0], 0); near(standardized.mean[1], 0); near(standardized.variance[0], 1); near(standardized.variance[1], 1); near(standardized.covariance, a.correlation);
  near(covarianceStatistics(points.map(([x, y]) => [-x, y])).correlation, -a.correlation);
  assert.throws(() => covarianceStatistics([[0, 1], [0, 2]]), RangeError);
});

test('Gaussian mixture is normalized, resolves its components, and can change its number of peaks', () => {
  for (const weight of [0, .2, .5, 1]) for (const d of [0, .5, 1.75]) near(integrate(x => gaussianMixtureDensity(x, weight, d), -10, 10), 1, 1e-8);
  near(gaussianMixtureDensity(.4, 1, 1), gaussianDensity(.4, -1, .65)); near(gaussianMixtureDensity(.4, 0, 1), gaussianDensity(.4, 1, .65));
  assert.ok(gaussianMixtureDensity(0, .5, .2) > gaussianMixtureDensity(.5, .5, .2));
  assert.ok(gaussianMixtureDensity(1.5, .5, 1.5) > gaussianMixtureDensity(0, .5, 1.5));
  near(gaussianMixtureDensity(1.5, .5, 1.5), gaussianMixtureDensity(-1.5, .5, 1.5));
  assert.throws(() => gaussianMixtureDensity(0, 1.1), RangeError); assert.throws(() => gaussianDensity(0, 0, 0), RangeError);
});

test('uniform maximum follows the independent-event product, derivative, normalization and expectation', () => {
  for (const n of [1, 2, 4, 8]) {
    near(uniformMaximumCDF(.6, n), Array(n).fill(.6).reduce((a, b) => a * b, 1));
    near(integrate(x => uniformMaximumDensity(x, n), 0, 1), 1, 1e-9);
    near(integrate(x => x * uniformMaximumDensity(x, n), 0, 1), uniformMaximumMean(n), 1e-9);
    const h = 1e-5; near((uniformMaximumCDF(.4 + h, n) - uniformMaximumCDF(.4 - h, n)) / (2 * h), uniformMaximumDensity(.4, n), 1e-8);
    near(uniformMaximumCDF(-1, n), 0); near(uniformMaximumCDF(2, n), 1);
  }
  assert.throws(() => uniformMaximumCDF(.4, 1.5), RangeError);
});

test('inverse CDF maps uniform probability mass to the triangular density', () => {
  for (let i = 0; i <= 100; i++) {const u = i / 100; near(triangularCDF(triangularQuantile(u)), u); near(triangularQuantile(triangularCDF(u)), u);}
  const a = .2, b = .7; near(triangularCDF(b) - triangularCDF(a), integrate(x => 2 * x, a, b));
  near(triangularCDF(-1), 0); near(triangularCDF(2), 1);
  assert.throws(() => triangularQuantile(-.01), RangeError);
});

test('Bernoulli KL is nonnegative, directed, and handles zero-probability support exactly', () => {
  for (let i = 0; i <= 30; i++) for (let j = 0; j <= 30; j++) { const p = i / 30, q = j / 30; assert.ok(bernoulliKL(p, q) >= 0); if (i === j) near(bernoulliKL(p, q), 0); }
  assert.ok(Math.abs(bernoulliKL(.3, .8) - bernoulliKL(.8, .3)) > .01);
  near(bernoulliKL(0, .3), -Math.log(.7)); near(bernoulliKL(1, .3), -Math.log(.3));
  assert.equal(bernoulliKL(.2, 0), Infinity); assert.equal(bernoulliKL(.8, 1), Infinity);
  near(bernoulliKL(.5, .25), .5 * Math.log(4 / 3));
  assert.throws(() => bernoulliKL(.3, NaN), RangeError);
});

test('probability renderers remain finite, in frame, and moving throughout repeated cycles', () => {
  let points = [], labels = [];
  const point = p => {assert.equal(p.length, 2); assert.ok(p.every(Number.isFinite)); assert.ok(Math.abs(p[0]) <= 5.3 && Math.abs(p[1]) <= 2.8, `out of frame ${p}`); points.push(p);};
  const helpers = {line: path => path.forEach(point), dot: point, label: (p, text) => {point(p); labels.push(text);}, curve: (fn, a, b) => { for (let i = 0; i <= 50; i++) point(fn(a + (b - a) * i / 50)); }};
  for (const id of probabilityRenderers) {
    const signatures = [];
    for (const time of [0, .1, 2.8, 8, 17, 101]) {points = []; labels = []; assert.equal(drawProbabilityConstruction({id, time, phase: (time * .22) % 1, progress: (1 - Math.cos(time * .35)) / 2, ...helpers}), true); assert.ok(points.length > 15); assert.ok(labels.length > 1); signatures.push(JSON.stringify(points));}
    assert.equal(new Set(signatures).size, signatures.length, id);
  }
  assert.equal(drawProbabilityConstruction({id: 'unknown'}), false);
});
