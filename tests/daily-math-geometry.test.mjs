import test from 'node:test';
import assert from 'node:assert/strict';
import {rightTriangleAreas, unityRoots, orthogonalProjection, cauchySchwarzValues, geometricRemainder, invertPoint, polygonPiBounds, barycentricPoint, triangleArea, geometryRenderers, drawGeometryConstruction} from '../src/scripts/daily-math-geometry.js';

const close = (actual, expected, tolerance = 1e-10) => assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} differs from ${expected}`);
const dot = (a, b) => a.reduce((sum, value, i) => sum + value * b[i], 0);
const difference = (a, b) => a.map((value, i) => value - b[i]);

test('right-triangle side squares agree with an independently measured hypotenuse', () => {
  for (const [a, b] of [[3, 4], [5, 12], [.1, .7], [Math.sqrt(2), Math.sqrt(3)]]) {
    const areas = rightTriangleAreas(a, b);
    close(areas.c ** 2, areas.a2 + areas.b2);
    close(areas.c2, Math.hypot(a, b) ** 2);
  }
  assert.equal(rightTriangleAreas(3, 4).c, 5);
  assert.throws(() => rightTriangleAreas(0, 2), RangeError);
});

test('unity roots satisfy z^n = 1 and nontrivial root sums cancel', () => {
  for (const n of [1, 2, 3, 7, 16]) {
    const roots = unityRoots(n);
    assert.equal(roots.length, n);
    for (const [x, y] of roots) {
      let real = 1, imaginary = 0;
      for (let k = 0; k < n; k++) [real, imaginary] = [real * x - imaginary * y, real * y + imaginary * x];
      close(real, 1); close(imaginary, 0);
      close(Math.hypot(x, y), 1);
    }
    if (n > 1) {close(roots.reduce((sum, p) => sum + p[0], 0), 0); close(roots.reduce((sum, p) => sum + p[1], 0), 0);}
  }
  assert.throws(() => unityRoots(0), RangeError);
  assert.throws(() => unityRoots(2.5), RangeError);
});

test('projection is idempotent, orthogonal, and nearest on the line', () => {
  const v = [2, -3, 5], u = [1, 2, -1], p = orthogonalProjection(v, u), residual = difference(v, p);
  close(dot(residual, u), 0);
  orthogonalProjection(p, u).forEach((value, i) => close(value, p[i]));
  orthogonalProjection(v, u.map(x => -7 * x)).forEach((value, i) => close(value, p[i]));
  const error = dot(residual, residual);
  for (const coefficient of [-10, -1, 0, 1, 3]) {
    const alternative = difference(v, u.map(x => coefficient * x));
    assert.ok(dot(alternative, alternative) >= error - 1e-10);
  }
  assert.deepEqual(orthogonalProjection([0, 0], [1, 0]), [0, 0]);
  assert.throws(() => orthogonalProjection([1, 2], [0, 0]), RangeError);
});

test('Cauchy–Schwarz includes equality, perpendicular and zero-vector cases', () => {
  for (const [u, v] of [[[3, 4], [-6, -8]], [[1, 0], [0, 2]], [[0, 0], [2, 3]], [[2, 1, -4], [7, -2, 3]]]) {
    const {absoluteDot, normProduct} = cauchySchwarzValues(u, v);
    assert.ok(absoluteDot <= normProduct + 1e-10);
  }
  close(cauchySchwarzValues([3, 4], [-6, -8]).absoluteDot, 50);
  close(cauchySchwarzValues([3, 4], [-6, -8]).normProduct, 50);
  assert.equal(cauchySchwarzValues([1, 0], [0, 2]).absoluteDot, 0);
  assert.throws(() => cauchySchwarzValues([1, 2], [1]), RangeError);
});

test('geometric remainders match independent summation and keep the sign of omitted tails', () => {
  for (const r of [0, .5, -.5, .85, -.8]) for (const n of [0, 1, 2, 7, 30]) {
    const result = geometricRemainder(r, n);
    let partial = 0, tail = 0;
    for (let k = 0; k < n; k++) partial += r ** k;
    for (let k = n; k < n + 500; k++) tail += r ** k;
    close(result.partial, partial); close(result.remainder, tail);
    close(result.partial + result.remainder, result.total);
    close(geometricRemainder(r, n + 1).remainder, result.remainder * r);
  }
  assert.ok(geometricRemainder(-.5, 1).remainder < 0);
  assert.throws(() => geometricRemainder(1, 3), RangeError);
  assert.throws(() => geometricRemainder(.5, -1), RangeError);
});

test('circle inversion is an involution and fixes its defining circle', () => {
  for (const radius of [.5, 1, 2]) for (const p of [[1, 0], [.3, -.6], [-2, 4], [3, 5]]) {
    const q = invertPoint(p, radius), twice = invertPoint(q, radius);
    twice.forEach((value, i) => close(value, p[i]));
    close(Math.hypot(...p) * Math.hypot(...q), radius ** 2);
    close(p[0] * q[1] - p[1] * q[0], 0);
    assert.ok(dot(p, q) > 0);
  }
  const p = [Math.sqrt(2), Math.sqrt(2)];
  invertPoint(p, 2).forEach((value, i) => close(value, p[i]));
  assert.throws(() => invertPoint([0, 0]), RangeError);
});

test('polygon perimeters enclose pi and doubling tightens both bounds', () => {
  let previous = {lower: 0, upper: Infinity};
  for (const n of [6, 12, 24, 48, 96, 192, 384]) {
    const bounds = polygonPiBounds(n);
    assert.ok(bounds.lower < Math.PI && bounds.upper > Math.PI);
    assert.ok(bounds.lower > previous.lower && bounds.upper < previous.upper);
    const roots = unityRoots(n);
    const measuredPerimeter = roots.reduce((sum, p, i) => sum + Math.hypot(...difference(p, roots[(i + 1) % n])), 0);
    close(bounds.lower, measuredPerimeter / 2);
    previous = bounds;
  }
  close(polygonPiBounds(6).lower, 3);
  close(polygonPiBounds(6).upper, 2 * Math.sqrt(3));
  assert.throws(() => polygonPiBounds(2), RangeError);
});

test('barycentric coordinates reproduce vertices and actual subtriangle area fractions', () => {
  const vertices = [[-3, -1.5], [2.5, -1.5], [-.5, 1.8]], [A, B, C] = vertices;
  const total = triangleArea(A, B, C);
  for (const weights of [[1, 0, 0], [0, 1, 0], [0, 0, 1], [.2, .3, .5], [1 / 3, 1 / 3, 1 / 3]]) {
    const P = barycentricPoint(vertices, weights);
    const fractions = [triangleArea(P, B, C), triangleArea(P, C, A), triangleArea(P, A, B)].map(area => area / total);
    fractions.forEach((value, i) => close(value, weights[i]));
  }
  assert.deepEqual(barycentricPoint(vertices, [1, 0, 0]), A);
  const weights = [.2, .3, .5], shifted = vertices.map(p => [p[0] + 3, p[1] - 7]);
  const before = barycentricPoint(vertices, weights), after = barycentricPoint(shifted, weights);
  close(after[0] - before[0], 3); close(after[1] - before[1], -7);
  assert.throws(() => barycentricPoint(vertices, [1, 1, 1]), RangeError);
});

test('each renderer keeps its construction finite and in the shared drawing window throughout animation', () => {
  const point = p => {assert.equal(p.length, 2); assert.ok(p.every(Number.isFinite)); assert.ok(Math.abs(p[0]) <= 5.01 && Math.abs(p[1]) <= 2.71, JSON.stringify(p));};
  const line = points => points.forEach(point), dot = point, label = point;
  const curve = (fn, a, b) => {for (let n = 0; n <= 24; n++) point(fn(a + (b - a) * n / 24));};
  const circle = (center, radius) => curve(angle => [center[0] + radius * Math.cos(angle), center[1] + radius * Math.sin(angle)], 0, Math.PI * 2);
  for (const id of geometryRenderers) for (let step = 0; step <= 160; step++) {
    assert.equal(drawGeometryConstruction({id, time: step * .25, line, dot, label, curve, circle, axes() {}}), true);
  }
  assert.equal(drawGeometryConstruction({id: 'not-a-geometry-topic'}), false);
});
