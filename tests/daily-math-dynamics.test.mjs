import {test} from 'node:test';
import assert from 'node:assert/strict';
import {logisticGrowth, twoStateProbability, binaryEntropy, exponentialSurvival, jensenQuadratic, eulerDecay, dampedState, varianceDecomposition, dynamicsRenderers, drawDynamicsConstruction} from '../src/scripts/daily-math-dynamics.js';
const near = (a, b, tolerance = 1e-10) => assert.ok(Math.abs(a - b) < tolerance, `${a} differs from ${b}`);

test('logistic solution satisfies its ODE, initial value, equilibria and capacity bound', () => {
  near(logisticGrowth(0, .6, 20, 3), 3);
  for (const t of [.2, 2, 8]) {
    const h = 1e-5, n = logisticGrowth(t, .6, 20, 3);
    near((logisticGrowth(t + h, .6, 20, 3) - logisticGrowth(t - h, .6, 20, 3)) / (2 * h), .6 * n * (1 - n / 20), 1e-8);
    assert.ok(n >= 3 && n < 20);
  }
  near(logisticGrowth(100, 1, 20, 0), 0); near(logisticGrowth(100, 1, 20, 20), 20);
  near(logisticGrowth(100, 0, 20, 3), 3); near(logisticGrowth(1000), 1);
  assert.throws(() => logisticGrowth(-1), RangeError);
});

test('two-state probabilities match the transition law and handle nonmixing boundaries', () => {
  for (const [a, b] of [[.2, .3], [.8, .7], [0, .4], [.4, 0]]) for (const initial of [0, .3, 1]) {
    near(twoStateProbability(0, initial, a, b), initial);
    for (let n = 0; n < 20; n++) near(twoStateProbability(n + 1, initial, a, b), b + (1 - a - b) * twoStateProbability(n, initial, a, b));
    near(twoStateProbability(200, initial, a, b), b / (a + b));
  }
  near(twoStateProbability(10, .4, 0, 0), .4);
  near(twoStateProbability(1, 1, 1, 1), 0); near(twoStateProbability(2, 1, 1, 1), 1);
  assert.throws(() => twoStateProbability(.5, .5), RangeError);
});

test('binary entropy is symmetric with maximum one bit and exact endpoint convention', () => {
  near(binaryEntropy(0), 0); near(binaryEntropy(1), 0); near(binaryEntropy(.5), 1);
  for (let i = 0; i <= 100; i++) {
    const p = i / 100; near(binaryEntropy(p), binaryEntropy(1 - p));
    assert.ok(binaryEntropy(p) >= 0 && binaryEntropy(p) <= 1);
  }
  assert.throws(() => binaryEntropy(-.1), RangeError); assert.throws(() => binaryEntropy(NaN), RangeError);
});

test('exponential conditional survival is unchanged by elapsed waiting time', () => {
  for (const rate of [.1, .7, 3]) for (const s of [0, 1, 5]) for (const t of [0, .3, 4]) {
    near(exponentialSurvival(s + t, rate) / exponentialSurvival(s, rate), exponentialSurvival(t, rate));
  }
  near(exponentialSurvival(0), 1); assert.throws(() => exponentialSurvival(2, 0), RangeError);
});

test('quadratic Jensen gap equals weighted squared separation with equality at boundaries', () => {
  for (const weight of [0, .1, .5, .9, 1]) for (const [a, b] of [[-2, 1.5], [2, 2]]) {
    const result = jensenQuadratic(a, b, weight);
    near(result.meanTransformed - result.transformedMean, result.gap); assert.ok(result.gap >= 0);
  }
  assert.throws(() => jensenQuadratic(1, 2, 1.1), RangeError);
});

test('Euler endpoint error is first order while one exact-start step has second-order error', () => {
  const coarse = eulerDecay(256), fine = eulerDecay(512);
  near(coarse.error / fine.error, 2, .01);
  near(coarse.points.at(-1)[0], 3);
  for (const steps of [4, 8, 16, 32]) {
    const result = eulerDecay(steps);
    near(result.points.at(-1)[1], (1 - 3 / steps) ** steps);
    assert.ok(result.points.every(([, value]) => value >= 0));
  }
  const error = h => Math.exp(-h) - (1 - h);
  near(error(.002) / error(.001), 4, .002);
  assert.throws(() => eulerDecay(0), RangeError);
});

test('damped motion obeys the ODE and dissipates exact mechanical energy', () => {
  const options = {gamma: .16, omega0: 1.7, x0: 1, v0: .2, mass: 2};
  near(dampedState(0, options).x, 1); near(dampedState(0, options).v, .2);
  for (const t of [.2, 1.7, 5, 12]) {
    const h = 1e-5, before = dampedState(t - h, options), now = dampedState(t, options), after = dampedState(t + h, options);
    near((after.x - before.x) / (2 * h), now.v, 1e-8);
    near((after.v - before.v) / (2 * h), -2 * options.gamma * now.v - options.omega0 ** 2 * now.x, 1e-8);
    near((after.energy - before.energy) / (2 * h), now.energyRate, 1e-8);
    assert.ok(now.energyRate <= 0 && Math.abs(now.x) <= now.amplitudeBound + 1e-12);
  }
  near(dampedState(10, {gamma: 0}).energy, dampedState(0, {gamma: 0}).energy);
  assert.throws(() => dampedState(1, {gamma: 2, omega0: 1}), RangeError);
});

test('total variance decomposes with observation-weighted unequal group sizes', () => {
  for (const groups of [[[1, 2], [4, 5, 6]], [[2], [2]], [[1, 2, 3]]]) {
    const result = varianceDecomposition(groups);
    near(result.within + result.between, result.total);
    near(result.groupStats.reduce((sum, group) => sum + group.weight, 0), 1);
  }
  const offsets = [-.7, -.35, 0, .35, .7];
  for (const d of [0, .2, 1.6]) {
    const result = varianceDecomposition([offsets.map(x => x - d), offsets.map(x => x + d)]);
    near(result.within, .245); near(result.between, d * d); near(result.total, .245 + d * d);
  }
  assert.throws(() => varianceDecomposition([[]]), RangeError);
});

test('all dynamics renderers produce finite in-frame geometry and reject unknown IDs', () => {
  let marks = 0;
  const point = p => {assert.equal(p.length, 2); assert.ok(p.every(Number.isFinite)); assert.ok(Math.abs(p[0]) <= 5.001 && Math.abs(p[1]) <= 2.701, `out of frame: ${p}`); marks++;};
  const helpers = {line: points => points.forEach(point), dot: point, label: point, curve: (fn, a, b) => {for (let i = 0; i <= 40; i++) point(fn(a + (b - a) * i / 40));}};
  for (const id of dynamicsRenderers) for (const time of [0, .1, 2.8, 8, 17, 101]) {
    marks = 0;
    assert.equal(drawDynamicsConstruction({id, time, phase: (time * .22) % 1, progress: (1 - Math.cos(time * .35)) / 2, ...helpers}), true);
    assert.ok(marks > 5, id);
  }
  assert.equal(drawDynamicsConstruction({id: 'unrecognized'}), false);
});
