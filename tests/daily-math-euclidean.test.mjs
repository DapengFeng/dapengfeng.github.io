import test from 'node:test';
import assert from 'node:assert/strict';
import {semicircleTriangle, angleBisectorPoint, circleChordThrough, perpendicularFoot, vivianiDistances, triangleCenters, parabolaFocusPoint, hyperbolaFocusPoint, sineRuleValues, euclideanRenderers, drawEuclideanConstruction} from '../src/scripts/daily-math-euclidean.js';

const sub = (a, b) => a.map((v, i) => v - b[i]);
const dot = (a, b) => a.reduce((sum, v, i) => sum + v * b[i], 0);
const cross = (a, b) => a[0] * b[1] - a[1] * b[0];
const distance = (a, b) => Math.hypot(...sub(a, b));
const close = (a, b, tolerance = 1e-10) => assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);

test('every nondegenerate point on a semicircle sees its diameter at a right angle', () => {
  for (const radius of [.3, 1, 4]) for (let i = 1; i < 40; i++) {
    const {A, B, P} = semicircleTriangle(radius, Math.PI * i / 40);
    close(dot(sub(A, P), sub(B, P)), 0);
    close(Math.hypot(...P), radius);
    close(distance(A, P) ** 2 + distance(B, P) ** 2, distance(A, B) ** 2);
  }
  assert.throws(() => semicircleTriangle(1, 0), RangeError);
  assert.throws(() => semicircleTriangle(-1, 1), RangeError);
});

test('angle bisector intersects the opposite side with the side-length ratio and equal angles', () => {
  for (const A of [[-1, 2], [2, 4], [-5, 1]]) {
    const B = [-3, -1], C = [3, -1], D = angleBisectorPoint(A, B, C);
    close(cross(sub(D, B), sub(C, B)), 0);
    assert.ok(D[0] > B[0] && D[0] < C[0]);
    close(distance(B, D) / distance(D, C), distance(A, B) / distance(A, C));
    const direction = sub(D, A), u = sub(B, A), v = sub(C, A);
    close(dot(direction, u) / (Math.hypot(...direction) * Math.hypot(...u)), dot(direction, v) / (Math.hypot(...direction) * Math.hypot(...v)));
  }
  assert.throws(() => angleBisectorPoint([0, 0], [1, 1], [2, 2]), RangeError);
});

test('intersecting chord products agree with the interior point power for every direction', () => {
  for (const radius of [.5, 2, 5]) for (const fraction of [0, .3, .8, .99]) {
    const P = [radius * fraction * Math.cos(.4), radius * fraction * Math.sin(.4)];
    for (let k = 0; k < 30; k++) {
      const chord = circleChordThrough(P, k * .21, radius);
      close(Math.hypot(...chord.A), radius); close(Math.hypot(...chord.B), radius);
      close(distance(P, chord.A), chord.PA); close(distance(P, chord.B), chord.PB);
      close(cross(sub(chord.A, P), sub(chord.B, P)), 0);
      assert.ok(dot(sub(chord.A, P), sub(chord.B, P)) <= 0);
      close(chord.PA * chord.PB, radius ** 2 - dot(P, P));
    }
  }
  assert.throws(() => circleChordThrough([1, 0], 0, 1), RangeError);
});

test('Viviani distances are perpendicular and sum to the altitude, including boundary points', () => {
  for (const side of [.5, 2, 4]) for (const weights of [[1, 0, 0], [.5, .5, 0], [.1, .7, .2], [1 / 3, 1 / 3, 1 / 3]]) {
    const {vertices, P, feet, distances, altitude} = vivianiDistances(side, weights);
    close(distances.reduce((sum, value) => sum + value, 0), altitude);
    close(altitude, side * Math.sqrt(3) / 2);
    feet.forEach((F, i) => {
      const A = vertices[i], B = vertices[(i + 1) % 3];
      close(dot(sub(P, F), sub(B, A)), 0); close(cross(sub(F, A), sub(B, A)), 0);
    });
  }
  assert.throws(() => vivianiDistances(2, [-.1, .5, .6]), RangeError);
  assert.throws(() => perpendicularFoot([0, 0], [1, 1], [1, 1]), RangeError);
});

test('Euler centers meet their independent definitions and retain the 1:2 division', () => {
  for (const vertices of [[[0, 0], [4, 0], [1, 3]], [[1, 1], [5, 1], [-1, 2]], [[0, 0], [3, 0], [0, 4]], [[-2, -1], [2, -1], [0, 2]]]) {
    const {O, G, H, radius} = triangleCenters(...vertices);
    vertices.forEach((A, i) => {
      close(distance(A, O), radius);
      const B = vertices[(i + 1) % 3], C = vertices[(i + 2) % 3];
      close(dot(sub(A, H), sub(C, B)), 0);
      const midpoint = B.map((value, axis) => (value + C[axis]) / 2);
      close(cross(sub(G, A), sub(midpoint, A)), 0);
    });
    close(cross(sub(G, O), sub(H, O)), 0);
    close(distance(G, H), 2 * distance(O, G));
  }
  const {O, G, H} = triangleCenters([0, 0], [2, 0], [1, Math.sqrt(3)]);
  close(distance(O, G), 0); close(distance(O, H), 0);
});

test('parabola points have equal focus and perpendicular directrix distances', () => {
  for (const p of [.1, .55, 2]) for (let y = -4; y <= 4; y += .25) {
    const {P, F, D, focusDistance, directrixDistance} = parabolaFocusPoint(p, y);
    close(P[1] ** 2, 4 * p * P[0]); close(D[0], -p); close(D[1], P[1]);
    close(distance(P, F), distance(P, D)); close(focusDistance, directrixDistance);
    close(focusDistance, P[0] + p);
  }
  assert.throws(() => parabolaFocusPoint(0, 1), RangeError);
});

test('both hyperbola branches preserve the absolute focal difference and quadratic equation', () => {
  for (const [a, b] of [[1, 1], [1.25, .85], [.3, 2]]) for (const branch of [-1, 1]) for (let t = -2; t <= 2; t += .125) {
    const {P, foci, distances, difference} = hyperbolaFocusPoint(a, b, t, branch);
    close(P[0] ** 2 / a ** 2 - P[1] ** 2 / b ** 2, 1);
    close(difference, 2 * a); close(distances[0] - distances[1], branch * 2 * a);
    foci.forEach((F, i) => close(distance(P, F), distances[i]));
  }
  assert.throws(() => hyperbolaFocusPoint(1, 1, 0, 0), RangeError);
  assert.throws(() => hyperbolaFocusPoint(1, 1, 1000), RangeError);
});

test('extended sine law matches the circumcircle for acute, obtuse and right triangles', () => {
  for (const vertices of [[[0, 0], [4, 0], [1, 3]], [[0, 0], [4, 0], [-1, 1]], [[0, 0], [3, 0], [0, 4]]]) {
    const {sides, sines, ratios, radius} = sineRuleValues(...vertices);
    ratios.forEach(value => close(value, 2 * radius));
    for (let i = 0; i < 3; i++) {
      const a = sides[i], b = sides[(i + 1) % 3], c = sides[(i + 2) % 3];
      const cosine = (b * b + c * c - a * a) / (2 * b * c);
      close(sines[i] ** 2 + cosine ** 2, 1);
    }
  }
  assert.throws(() => sineRuleValues([0, 0], [1, 1], [2, 2]), RangeError);
});

function frame(id, time) {
  const points = [];
  const point = p => {
    assert.ok(Array.isArray(p) && p.length === 2 && p.every(Number.isFinite));
    assert.ok(Math.abs(p[0]) <= 5.3 && Math.abs(p[1]) <= 2.8, `${id}: ${JSON.stringify(p)}`);
    points.push(p);
  };
  const curve = (fn, a, b) => {for (let i = 0; i <= 24; i++) point(fn(a + (b - a) * i / 24));};
  assert.equal(drawEuclideanConstruction({id, time, line: points => points.forEach(point), dot: point, label: point, curve, circle: (p, r) => curve(t => [p[0] + r * Math.cos(t), p[1] + r * Math.sin(t)], 0, 2 * Math.PI)}), true);
  return points.flat();
}

test('all eight constructions stay finite and visible, animate continuously, and do not stop', () => {
  for (const id of euclideanRenderers) {
    for (let step = 0; step <= 200; step++) frame(id, step * .31);
    for (const t of [0, 5, 17, 63, 180, 3600]) {
      const a = frame(id, t), b = frame(id, t + 1e-5);
      assert.equal(a.length, b.length);
      assert.ok(a.some((x, i) => Math.abs(x - b[i]) > 1e-9), `${id} has stopped at ${t}`);
      assert.ok(a.every((x, i) => Math.abs(x - b[i]) < 1e-3), `${id} jumps at ${t}`);
    }
  }
  assert.equal(drawEuclideanConstruction({id: 'unknown'}), false);
});
