const TAU = 2 * Math.PI;
const GREEN = '#c5ec96', GOLD = '#e5bf82', CYAN = '#83bec8', FAINT = '#9ab79d45';
export const geometryRenderers = ['pythagorean-area', 'roots-of-unity', 'orthogonal-projection', 'cauchy-schwarz', 'geometric-remainder', 'circle-inversion', 'polygon-pi', 'barycentric-coordinates'];

const inner = (a, b) => a.reduce((sum, value, i) => sum + value * b[i], 0);
export function rightTriangleAreas(a, b) {
  if (!(a > 0 && b > 0) || !Number.isFinite(a + b)) throw new RangeError('Positive finite legs required');
  return {a2: a * a, b2: b * b, c: Math.hypot(a, b), c2: a * a + b * b};
}
export function unityRoots(n) {
  if (!Number.isInteger(n) || n < 1) throw new RangeError('Positive integer order required');
  return Array.from({length: n}, (_, k) => [Math.cos(TAU * k / n), Math.sin(TAU * k / n)]);
}
export function orthogonalProjection(v, u) {
  if (!v.length || v.length !== u.length || ![...v, ...u].every(Number.isFinite)) throw new RangeError('Equal finite vector dimensions required');
  const norm2 = inner(u, u);
  if (norm2 === 0) throw new RangeError('Projection direction must be nonzero');
  const coefficient = inner(v, u) / norm2;
  return u.map(value => coefficient * value);
}
export function cauchySchwarzValues(u, v) {
  if (!u.length || u.length !== v.length || ![...u, ...v].every(Number.isFinite)) throw new RangeError('Equal finite vector dimensions required');
  return {absoluteDot: Math.abs(inner(u, v)), normProduct: Math.hypot(...u) * Math.hypot(...v)};
}
// n is the number of retained terms: 1 + r + ... + r^(n-1).
export function geometricRemainder(r, n) {
  if (!Number.isFinite(r) || Math.abs(r) >= 1 || !Number.isInteger(n) || n < 0) throw new RangeError('|r| < 1 and a nonnegative term count required');
  return {partial: (1 - r ** n) / (1 - r), remainder: r ** n / (1 - r), total: 1 / (1 - r)};
}
export function invertPoint(point, radius = 1) {
  if (point.length !== 2 || !point.every(Number.isFinite) || !(radius > 0) || !Number.isFinite(radius)) throw new RangeError('Finite plane point and positive radius required');
  const length2 = inner(point, point);
  if (length2 === 0) throw new RangeError('The inversion center has no finite image');
  return point.map(value => radius * radius * value / length2);
}
export function polygonPiBounds(n) {
  if (!Number.isInteger(n) || n < 3) throw new RangeError('A polygon needs at least three sides');
  return {lower: n * Math.sin(Math.PI / n), upper: n * Math.tan(Math.PI / n)};
}
export function barycentricPoint(vertices, weights) {
  if (vertices.length !== 3 || weights.length !== 3 || vertices.some(p => p.length !== 2 || !p.every(Number.isFinite)) || !weights.every(Number.isFinite)) throw new RangeError('Three finite plane vertices and weights required');
  if (Math.abs(weights.reduce((a, b) => a + b, 0) - 1) > 1e-10) throw new RangeError('Normalized weights must sum to one');
  return [0, 1].map(axis => weights.reduce((sum, weight, i) => sum + weight * vertices[i][axis], 0));
}
export function triangleArea(a, b, c) {
  return Math.abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2;
}

export function drawGeometryConstruction({id, time, line, dot, label, curve, circle, axes}) {
  if (!geometryRenderers.includes(id)) return false;
  const polygon = (points, color, fill) => line([...points, points[0]], color, 1.8, fill);
  const add = (a, b) => [a[0] + b[0], a[1] + b[1]];
  if (id === 'pythagorean-area') {
    const a = 1.1 + .3 * Math.sin(time * .45), b = 1.2 - .2 * Math.sin(time * .45);
    const A = [-1.3, -.4], B = [A[0] + a, A[1]], C = [A[0], A[1] + b];
    polygon([A, B, [B[0], B[1] - a], [A[0], A[1] - a]], GREEN, '#c5ec9620');
    polygon([A, C, [C[0] - b, C[1]], [A[0] - b, A[1]]], CYAN, '#83bec820');
    polygon([B, C, add(C, [b, a]), add(B, [b, a])], GOLD, '#e5bf8220');
    polygon([A, B, C], FAINT);
    line([add(A, [.16, 0]), add(A, [.16, .16]), add(A, [0, .16])], FAINT, 1);
    label([A[0] + a / 2 - .18, A[1] - a / 2], 'a²', GREEN);
    label([A[0] - b / 2 - .18, A[1] + b / 2], 'b²', CYAN);
    label([(B[0] + C[0] + b) / 2 - .18, (B[1] + C[1] + a) / 2], 'c²', GOLD);
    const areas = rightTriangleAreas(a, b), scale = 5.6 / areas.c2;
    line([[-2.8, -2.2], [-2.8 + areas.a2 * scale, -2.2]], GREEN, 7);
    line([[-2.8 + areas.a2 * scale, -2.2], [2.8, -2.2]], CYAN, 7);
    line([[-2.8, -2.5], [2.8, -2.5]], GOLD, 7);
    label([-2.6, 2.55], 'a² + b² = c²');
  } else if (id === 'roots-of-unity') {
    const roots = unityRoots(7), center = [-1.5, 0], radius = 1.85;
    const points = roots.map(p => add(center, p.map(v => v * radius))), step = Math.floor(time * .65) % 7;
    circle(center, radius, FAINT); polygon(points, CYAN);
    points.forEach((point, k) => {dot(point, CYAN, 3); if (k === step) {line([center, point], GREEN, 2.5); dot(point, GREEN, 5 + Math.sin(time * 4));}});
    line([center, points[0]], GOLD, 1.2);
    curve(angle => add(center, [.5 * Math.cos(angle), .5 * Math.sin(angle)]), 0, TAU / 7, GOLD, 35);
    label([1.2, 1.3], 'n = 7', CYAN); label([1.2, .55], `k = ${step}`, GREEN);
    label([1.2, -.2], 'z⁷ = 1', GOLD); label([-2.3, -2.5], 'Δθ = 2π/7', GOLD);
  } else if (id === 'orthogonal-projection') {
    const u = [1, .35], v = [2.3 * Math.cos(time * .38), 1.8 * Math.sin(time * .38)], p = orthogonalProjection(v, u);
    const anchor = [-.6, 0], at = point => add(point, anchor);
    line([at([-3.8, -1.33]), at([3.8, 1.33])], FAINT);
    line([anchor, at(v)], CYAN, 2); line([anchor, at(p)], GREEN, 3); line([at(p), at(v)], GOLD, 2);
    const norm = Math.hypot(...u), du = u.map(x => .17 * x / norm), sign = v[1] - p[1] < 0 ? -1 : 1, dn = [-du[1] * sign, du[0] * sign];
    line([at(add(p, du)), at(add(add(p, du), dn)), at(add(p, dn))], GOLD, 1);
    dot(at(v), CYAN, 5); dot(at(p), GREEN, 4);
    label([-4.3, 2.5], 'v', CYAN); label([-2.8, 2.5], 'p', GREEN); label([-1.3, 2.5], 'v − p', GOLD);
    label([-3.4, -2.5], '(v − p) · u = 0', GOLD);
  } else if (id === 'cauchy-schwarz') {
    const theta = time * .5, u = [1.8, 0], v = [1.6 * Math.cos(theta), 1.6 * Math.sin(theta)], center = [-2.1, .35];
    circle(center, 1.6, FAINT); line([center, add(center, u)], GOLD, 2.5); line([center, add(center, v)], GREEN, 2.5);
    dot(add(center, v), GREEN, 5); line([add(center, v), add(center, [v[0], 0])], CYAN, 1.5);
    const {absoluteDot, normProduct} = cauchySchwarzValues(u, v), x = 1.2, length = 2.6;
    line([[x, .8], [x + length, .8]], GOLD, 7);
    line([[x, -.1], [x + length * absoluteDot / normProduct, -.1]], GREEN, 7);
    label([x, 1.25], '‖u‖ ‖v‖', GOLD); label([x, -.75], '|u · v|', GREEN);
    label([-3.4, -2.25], `|cos θ| = ${Math.abs(Math.cos(theta)).toFixed(2)}`, CYAN);
  } else if (id === 'geometric-remainder') {
    const n = 1 + Math.floor(time * .6) % 7, r = .5, start = -3.8, width = 7.6;
    const {partial, remainder, total} = geometricRemainder(r, n);
    let sum = 0;
    for (let k = 0; k < n; k++) {
      const next = sum + r ** k, left = start + width * sum / total, right = start + width * next / total;
      polygon([[left, -.45], [right, -.45], [right, .8], [left, .8]], GREEN, '#c5ec9622'); sum = next;
    }
    const boundary = start + width * partial / total;
    polygon([[boundary, -.45], [start + width, -.45], [start + width, .8], [boundary, .8]], GOLD, '#e5bf8235');
    const pulse = (Math.sin(time * 3) + 1) / 2;
    dot([boundary, -1.1], GOLD, 3 + pulse * 2);
    label([-3.8, 2], `n = ${n}`, GREEN); label([-.1, 2], 'r = 1/2', CYAN);
    label([-3.8, -2.1], `S = ${partial.toFixed(4)}`, GREEN); label([.2, -2.1], `R = ${remainder.toFixed(4)}`, GOLD);
  } else if (id === 'circle-inversion') {
    const radius = 1.35, angle = time * .36, distance = 1.3 + .55 * Math.sin(time * .7);
    const p = [distance * Math.cos(angle), distance * Math.sin(angle)], q = invertPoint(p, radius);
    circle([0, 0], radius, CYAN); line([[0, 0], p], GREEN, 2); line([[0, 0], q], GOLD, 1.1);
    dot(p, GREEN, 6); dot(q, GOLD, 4); dot([0, 0], CYAN, 2.5);
    label([-4, 2.5], 'P', GREEN); label([-2.8, 2.5], 'P′', GOLD); label([-1.5, 2.5], 'R = 1.35', CYAN);
    label([-3.2, -2.6], `OP·OP′ = ${(radius * radius).toFixed(4)}`, GOLD);
  } else if (id === 'polygon-pi') {
    const n = [6, 12, 24, 48, 96][Math.floor(time * .42) % 5], center = [-2.1, 0], radius = 1.65;
    const inner = unityRoots(n).map(p => add(center, p.map(v => v * radius)));
    const outer = Array.from({length: n}, (_, k) => {const angle = TAU * (k + .5) / n; return add(center, [Math.cos(angle) * radius / Math.cos(Math.PI / n), Math.sin(angle) * radius / Math.cos(Math.PI / n)]);});
    polygon(outer, GOLD); circle(center, radius, CYAN); polygon(inner, GREEN);
    const angle = time * .5; dot(add(center, [radius * Math.cos(angle), radius * Math.sin(angle)]), CYAN, 4);
    const {lower, upper} = polygonPiBounds(n);
    label([.6, 1.3], `n = ${n}`, CYAN); label([.6, .4], lower.toFixed(5), GREEN);
    label([.6, -.35], '< π <', CYAN); label([.6, -1.1], upper.toFixed(5), GOLD);
  } else if (id === 'barycentric-coordinates') {
    const A = [-3, -1.5], B = [2.5, -1.5], C = [-.5, 1.8], vertices = [A, B, C];
    const raw = [0, TAU / 3, 2 * TAU / 3].map(offset => 1 + .8 * Math.sin(time * .45 + offset)), weights = raw.map(v => v / 3), P = barycentricPoint(vertices, weights);
    polygon([P, B, C], GREEN, '#c5ec9620'); polygon([P, C, A], GOLD, '#e5bf8220'); polygon([P, A, B], CYAN, '#83bec820');
    polygon(vertices, FAINT); dot(P, '#ffffff', 4);
    label([-3.4, -1.9], 'A', GREEN); label([2.55, -1.9], 'B', GOLD); label([-.6, 2.15], 'C', CYAN);
    weights.forEach((weight, index) => label([-3.7 + index * 2.65, -2.6], `λ${['A', 'B', 'C'][index]}=${weight.toFixed(2)}`, [GREEN, GOLD, CYAN][index]));
  }
  return true;
}
