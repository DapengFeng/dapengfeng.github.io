import {test} from 'node:test';
import assert from 'node:assert/strict';
import {accumulatedArea, compositeTrapezoid, simpsonPanel, secantStep, secantSequence, lagrangeValue, interpolationNodes, rungeInterpolation, parabolaCurvature, parabolaArcLength, analysisRenderers, drawAnalysisConstruction} from '../src/scripts/daily-math-analysis.js';
const near = (a, b, tolerance = 1e-10) => assert.ok(Math.abs(a - b) < tolerance, `${a} differs from ${b}`);

test('accumulated area starts at zero and its numerical derivative equals the integrand', () => {
  near(accumulatedArea(0).value, 0);
  for (const x of [.1, .5, 1, 1.9]) {
    const h = 1e-5, result = accumulatedArea(x);
    near((accumulatedArea(x + h).value - accumulatedArea(x - h).value) / (2 * h), result.derivative, 1e-9);
    near(simpsonPanel(t => .5 + t * t / 4, 0, x), result.value);
  }
  assert.throws(() => accumulatedArea(Infinity), RangeError);
});

test('trapezoidal integration is exact for affine functions and quadratic error quarters on refinement', () => {
  for (const n of [1, 2, 7, 16]) {
    near(compositeTrapezoid(x => 2 * x + 3, -1, 2, n), 12);
    const error = compositeTrapezoid(x => x * x, 0, 1, n) - 1 / 3;
    near(error, 1 / (6 * n * n));
    near(error / (compositeTrapezoid(x => x * x, 0, 1, 2 * n) - 1 / 3), 4, 1e-10);
  }
  assert.throws(() => compositeTrapezoid(x => x, 0, 1, 0), RangeError);
  assert.throws(() => compositeTrapezoid(() => NaN, 0, 1, 4), RangeError);
});

test('Simpson integrates every monomial up to cubic exactly, including asymmetric intervals', () => {
  for (const [a, b] of [[0, 1], [-1.3, 2.1], [2, 5]]) for (let power = 0; power <= 3; power++) {
    near(simpsonPanel(x => x ** power, a, b), (b ** (power + 1) - a ** (power + 1)) / (power + 1));
  }
  near(simpsonPanel(x => x ** 4, 0, 1) - .2, 1 / 120);
  const residualPrimitive = x => x ** 4 / 4 - x ** 3 / 2 + x ** 2 / 4;
  near(residualPrimitive(.5), 1 / 64); near(residualPrimitive(1), 0);
  assert.throws(() => simpsonPanel(x => x, 1, 0), RangeError);
});

test('secant iterates use two evaluations, find the line root, and approach sqrt(2)', () => {
  const f = x => x * x - 2, points = secantSequence(f, 1, 2, 5);
  near(points[2], 4 / 3); near(points[3], 7 / 5); near(points.at(-1), Math.sqrt(2), 1e-8);
  for (let i = 2; i < points.length; i++) {
    const a = points[i - 2], b = points[i - 1], next = points[i];
    near(f(b) + (f(b) - f(a)) * (next - b) / (b - a), 0, 1e-10);
  }
  near(secantStep(x => 3 * x - 7, -8, 10), 7 / 3);
  assert.throws(() => secantStep(() => 2, 1, 2), RangeError);
});

test('Lagrange interpolation reproduces nodes, constants and independent polynomial coefficients', () => {
  const nodes = [-2, -.3, .8, 2.1], f = x => .7 * x ** 3 - 2 * x * x + .3 * x - 1;
  const values = nodes.map(f);
  for (const x of [...nodes, -1.5, 0, 1.25]) near(lagrangeValue(nodes, values, x), f(x));
  for (const x of [-3, 0, 2]) near(lagrangeValue(nodes, [1, 1, 1, 1], x), 1);
  // A difference of degree <= 3 with these four roots must be identically zero.
  for (const x of [-1.7, -.2, .4, 2.4]) near(lagrangeValue(nodes, values.map(v => 2 * v), x), 2 * f(x));
  assert.throws(() => lagrangeValue([0, 0], [1, 2], 1), RangeError);
});

test('Chebyshev endpoints cluster nodes and reduce Runge error at the displayed degrees', () => {
  for (const n of [6, 8, 10]) {
    const nodes = interpolationNodes(n, true), equal = rungeInterpolation(n), cheb = rungeInterpolation(n, true);
    near(nodes[0], -1); near(nodes.at(-1), 1); near(nodes[n / 2], 0);
    assert.ok(nodes[1] - nodes[0] < nodes[n / 2] - nodes[n / 2 - 1]);
    let equalError = 0, chebError = 0;
    for (let i = 0; i <= 2000; i++) {
      const x = -1 + i / 1000, exact = 1 / (1 + 25 * x * x);
      equalError = Math.max(equalError, Math.abs(equal.evaluate(x) - exact));
      chebError = Math.max(chebError, Math.abs(cheb.evaluate(x) - exact));
      // No hidden plot clipping: this bounds every plotted degree's full range.
      assert.ok(equal.evaluate(x) > -1.01 && equal.evaluate(x) < 1.97);
    }
    assert.ok(chebError < equalError, `${n}: ${chebError} >= ${equalError}`);
    nodes.forEach((x, k) => near(cheb.evaluate(x), cheb.values[k]));
  }
});

test('osculating circle touches the parabola normally and shares its curvature', () => {
  for (const x of [-1, -.55, 0, .27, .55, 1]) {
    const {point, center, radius, curvature} = parabolaCurvature(x), normal = [center[0] - point[0], center[1] - point[1]];
    near(Math.hypot(...normal), radius); near(normal[0] + 2 * x * normal[1], 0); near(radius * curvature, 1);
    assert.ok(normal[1] > 0);
    // Circle's local graph y=c_y-sqrt(R²-(u-c_x)²) shares first and second derivatives.
    const circleSlope = (x - center[0]) / (center[1] - x * x);
    const circleSecond = radius ** 2 / (center[1] - x * x) ** 3;
    near(circleSlope, 2 * x); near(circleSecond, 2);
  }
  near(parabolaCurvature(0).radius, .5);
});

test('nested polygonal lengths increase to the integral and agree with analytic endpoint limits', () => {
  let previous = 0;
  for (const n of [1, 2, 4, 8, 16, 32, 64, 1024]) {
    const {length, exact} = parabolaArcLength(-1.5, 1.5, n);
    assert.ok(length > previous && length < exact); previous = length;
    near(exact, 1.5 * Math.sqrt(3.25) + Math.log(1.5 + Math.sqrt(3.25)));
    if (n === 1024) near(length, exact, 1e-6);
  }
  near(parabolaArcLength(0, 2, 1).length, Math.sqrt(8));
  assert.throws(() => parabolaArcLength(0, 1, 1.5), RangeError);
});

test('every analysis drawing stays finite in frame and has changing mathematical marks', () => {
  const snapshots = new Map();
  const inspect = (id, time) => {
    const marks = [], point = p => {assert.equal(p.length, 2); assert.ok(p.every(Number.isFinite)); assert.ok(Math.abs(p[0]) <= 5.3 && Math.abs(p[1]) <= 2.8, `${id} at ${time}: ${p}`); marks.push(p);};
    const helpers = {line: points => points.forEach(point), dot: point, label: point, curve: (fn, a, b) => {for (let i = 0; i <= 100; i++) point(fn(a + (b - a) * i / 100));}, circle: (p, r) => {assert.ok(r > 0); for (let i = 0; i < 80; i++) point([p[0] + r * Math.cos(i * Math.PI / 40), p[1] + r * Math.sin(i * Math.PI / 40)]);}};
    assert.equal(drawAnalysisConstruction({id, time, phase: time * .22 % 1, progress: (1 - Math.cos(time * .35)) / 2, ...helpers}), true);
    assert.ok(marks.length > 15); return JSON.stringify(marks);
  };
  for (const id of analysisRenderers) {
    for (const time of [0, .01, .4, 2.8, 5.99, 6.01, 8, 10, 17, 21, 31, 63, 120]) inspect(id, time);
    snapshots.set(id, inspect(id, 1.2)); assert.notEqual(snapshots.get(id), inspect(id, 1.25), `${id} must animate within a stage`);
  }
  assert.equal(drawAnalysisConstruction({id: 'unknown'}), false);
});
