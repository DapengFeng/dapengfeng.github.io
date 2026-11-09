const TAU = 2 * Math.PI;
const GREEN = '#c5ec96', GOLD = '#e5bf82', CYAN = '#83bec8', FAINT = '#9ab79d45';
export const euclideanRenderers = ['thales-angle', 'angle-bisector', 'intersecting-chords', 'viviani-distance', 'triangle-euler-line', 'parabola-focus', 'hyperbola-foci', 'sine-rule'];

const add = (a, b) => a.map((x, i) => x + b[i]);
const sub = (a, b) => a.map((x, i) => x - b[i]);
const mul = (a, s) => a.map(x => x * s);
const inner = (a, b) => a[0] * b[0] + a[1] * b[1];
const norm = a => Math.hypot(...a);
const distance = (a, b) => norm(sub(a, b));
const cross = (a, b) => a[0] * b[1] - a[1] * b[0];
const finitePoint = p => Array.isArray(p) && p.length === 2 && p.every(Number.isFinite);
function triangle(A, B, C) {
  if (![A, B, C].every(finitePoint) || Math.abs(cross(sub(B, A), sub(C, A))) < 1e-12) throw new RangeError('Three finite non-collinear vertices required');
}

export function semicircleTriangle(radius, theta) {
  if (!(radius > 0 && theta > 0 && theta < Math.PI) || !Number.isFinite(radius)) throw new RangeError('Positive finite radius and angle strictly between 0 and pi required');
  return {A: [-radius, 0], B: [radius, 0], P: [radius * Math.cos(theta), radius * Math.sin(theta)]};
}

export function angleBisectorPoint(A, B, C) {
  triangle(A, B, C);
  const AB = distance(A, B), AC = distance(A, C);
  // The ratio BD/DC = AB/AC places D this fraction along BC.
  return add(B, mul(sub(C, B), AB / (AB + AC)));
}

export function circleChordThrough(point, angle, radius = 1) {
  if (!finitePoint(point) || !Number.isFinite(angle) || !(radius > 0) || !Number.isFinite(radius) || norm(point) >= radius) throw new RangeError('Finite angle and point strictly inside a positive circle required');
  const direction = [Math.cos(angle), Math.sin(angle)], projection = inner(point, direction);
  const discriminant = Math.sqrt(projection ** 2 + radius ** 2 - inner(point, point));
  const negative = -projection - discriminant, positive = -projection + discriminant;
  return {A: add(point, mul(direction, negative)), B: add(point, mul(direction, positive)), PA: -negative, PB: positive};
}

export function perpendicularFoot(P, A, B) {
  if (![P, A, B].every(finitePoint) || distance(A, B) === 0) throw new RangeError('Finite point and a nonzero line required');
  const direction = sub(B, A);
  return add(A, mul(direction, inner(sub(P, A), direction) / inner(direction, direction)));
}

export function vivianiDistances(side, weights) {
  if (!(side > 0) || !Number.isFinite(side) || !Array.isArray(weights) || weights.length !== 3 || !weights.every(x => Number.isFinite(x) && x >= 0) || Math.abs(weights.reduce((a, b) => a + b, 0) - 1) > 1e-10) throw new RangeError('Positive side and three nonnegative normalized weights required');
  const altitude = Math.sqrt(3) * side / 2, vertices = [[-side / 2, 0], [side / 2, 0], [0, altitude]];
  const P = [0, 1].map(axis => vertices.reduce((sum, v, i) => sum + weights[i] * v[axis], 0));
  const feet = vertices.map((v, i) => perpendicularFoot(P, v, vertices[(i + 1) % 3]));
  return {vertices, P, feet, distances: feet.map(p => distance(P, p)), altitude};
}

export function triangleCenters(A, B, C) {
  triangle(A, B, C);
  // Solve perpendicular-bisector equations in coordinates relative to A.
  const b = sub(B, A), c = sub(C, A), determinant = 2 * cross(b, c);
  const b2 = inner(b, b), c2 = inner(c, c);
  const O = add(A, [(b2 * c[1] - c2 * b[1]) / determinant, (b[0] * c2 - c[0] * b2) / determinant]);
  const G = mul(add(add(A, B), C), 1 / 3);
  const H = sub(add(add(A, B), C), mul(O, 2));
  return {O, G, H, radius: distance(O, A)};
}

export function parabolaFocusPoint(p, y) {
  if (!(p > 0) || !Number.isFinite(p + y)) throw new RangeError('Positive finite focal parameter and finite ordinate required');
  const P = [y * y / (4 * p), y], F = [p, 0], D = [-p, y];
  return {P, F, D, focusDistance: distance(P, F), directrixDistance: distance(P, D)};
}

export function hyperbolaFocusPoint(a, b, parameter, branch = 1) {
  if (!(a > 0 && b > 0) || !Number.isFinite(a + b + parameter) || ![-1, 1].includes(branch)) throw new RangeError('Positive finite semi-axes, finite parameter and branch ±1 required');
  const c = Math.hypot(a, b), P = [branch * a * Math.cosh(parameter), b * Math.sinh(parameter)];
  if (!P.every(Number.isFinite)) throw new RangeError('Hyperbola point exceeds the numerical range');
  const foci = [[-c, 0], [c, 0]], distances = foci.map(F => distance(P, F));
  return {P, foci, distances, difference: Math.abs(distances[0] - distances[1])};
}

export function sineRuleValues(A, B, C) {
  triangle(A, B, C);
  const vertices = [A, B, C], sides = [distance(B, C), distance(C, A), distance(A, B)];
  const sines = vertices.map((v, i) => {
    const u = sub(vertices[(i + 1) % 3], v), w = sub(vertices[(i + 2) % 3], v);
    return Math.abs(cross(u, w)) / (norm(u) * norm(w));
  });
  return {sides, sines, ratios: sides.map((side, i) => side / sines[i]), radius: triangleCenters(A, B, C).radius};
}

export function drawEuclideanConstruction({id, time, line, dot, label, curve, circle}) {
  if (!euclideanRenderers.includes(id)) return false;
  const polygon = (points, color = FAINT) => line([...points, points[0]], color, 1.5);
  const angleArc = (vertex, A, B, radius, color) => {
    const u = sub(A, vertex), v = sub(B, vertex), start = Math.atan2(u[1], u[0]);
    const sweep = Math.atan2(cross(u, v), inner(u, v));
    curve(t => add(vertex, [radius * Math.cos(t), radius * Math.sin(t)]), start, start + sweep, color, 30);
  };
  const rightMark = (vertex, A, B, color) => {
    const u = mul(sub(A, vertex), .16 / distance(A, vertex)), v = mul(sub(B, vertex), .16 / distance(B, vertex));
    line([add(vertex, u), add(add(vertex, u), v), add(vertex, v)], color, 1);
  };
  if (id === 'thales-angle') {
    const radius = 2.4, theta = Math.PI / 2 + 1.1 * Math.sin(time * .35), center = [-.7, -.7];
    const t = semicircleTriangle(radius, theta), [A, B, P] = [t.A, t.B, t.P].map(v => add(center, v));
    curve(t => add(center, [radius * Math.cos(t), radius * Math.sin(t)]), 0, Math.PI, FAINT);
    line([A, B], CYAN); line([A, P, B], GREEN, 2.2); rightMark(P, A, B, GOLD);
    [A, B].forEach(v => dot(v, CYAN)); dot(P, GOLD, 5);
    label(add(A, [-.15, -.75]), 'A', CYAN); label(add(B, [-.15, -.75]), 'B', CYAN);
    label(add(P, [-.08, .4]), 'P', GOLD); label([-2, -2.25], '∠APB = 90°', GOLD);
  } else if (id === 'angle-bisector') {
    const A = [-.8 + .65 * Math.sin(time * .4), 1.8], B = [-3.2, -1.4], C = [1.4, -1.4], D = angleBisectorPoint(A, B, C);
    polygon([A, B, C]); line([A, B], GREEN, 2); line([A, C], CYAN, 2); line([A, D], GOLD, 2);
    line([B, D], GREEN, 3); line([D, C], CYAN, 3);
    angleArc(A, B, D, .4, GOLD); angleArc(A, D, C, .52, GOLD);
    dot(D, GOLD, 4); label(add(A, [-.1, .4]), 'A'); label(add(B, [-.15, -.75]), 'B', GREEN);
    label(add(C, [-.15, -.75]), 'C', CYAN); label(add(D, [-.15, -.75]), 'D');
    const ratio = distance(A, B) / distance(A, C);
    label([2.3, 1], 'AB/AC', GREEN); label([2.3, .45], ratio.toFixed(3), GREEN);
    label([2.3, -.45], 'BD/DC', CYAN); label([2.3, -1], ratio.toFixed(3), CYAN);
  } else if (id === 'intersecting-chords') {
    const center = [-1.75, 0], radius = 1.9, P = [.65 * Math.cos(time * .3), .45 * Math.sin(time * .3)];
    const first = circleChordThrough(P, time * .2 + .2, radius), second = circleChordThrough(P, time * .2 + 1.5, radius);
    circle(center, radius, FAINT);
    line([center, add(center, P)], FAINT, 1); dot(center, CYAN, 2);
    line([center, add(center, [0, .7])], FAINT, 1); label(add(center, [-.15, 1]), 'O', CYAN);
    [first, second].forEach((chord, i) => {
      const color = [GREEN, GOLD][i]; line([add(center, chord.A), add(center, chord.B)], color, 2);
      [chord.A, chord.B].forEach((v, j) => {
        dot(add(center, v), color);
        label(add(add(center, v), [.6 * v[0] / radius - .16, .6 * v[1] / radius - .15]), ['A', 'B', 'C', 'D'][2 * i + j], color);
      });
    });
    dot(add(center, P), CYAN, 5); label(add(add(center, P), [.32, .42]), 'P', CYAN);
    label([1.1, 1.2], 'PA · PB', GREEN); label([1.1, .65], (first.PA * first.PB).toFixed(3), GREEN);
    label([1.1, -.35], 'PC · PD', GOLD); label([1.1, -.9], (second.PA * second.PB).toFixed(3), GOLD);
    label([1.1, -2.3], 'R²−OP²', CYAN);
  } else if (id === 'viviani-distance') {
    const raw = [0, TAU / 3, 2 * TAU / 3].map(offset => 1 + .6 * Math.sin(time * .4 + offset));
    const result = vivianiDistances(4, raw.map(v => v / 3)), offset = [-1.35, -1.35];
    polygon(result.vertices.map(v => add(offset, v)), FAINT);
    let base = offset[1];
    result.feet.forEach((foot, i) => {
      const color = [GREEN, GOLD, CYAN][i], P = add(offset, result.P), F = add(offset, foot);
      line([P, F], color, 2.5); dot(F, color, 3);
      const side = result.vertices[i];
      if (distance(foot, side) > .02 && distance(foot, result.P) > .02) rightMark(F, add(offset, side), P, color);
      line([[1.7, base], [1.7, base + result.distances[i]]], color, 7);
      label([1.95, base + result.distances[i] / 2 - .06], `d${['₁', '₂', '₃'][i]}`, color);
      base += result.distances[i];
    });
    dot(add(offset, result.P), '#ffffff', 4);
    line([[3.3, offset[1]], [3.3, offset[1] + result.altitude]], FAINT, 7); label([3.55, .3], 'h', GOLD);
    label([-2.5, -2.3], 'd₁ + d₂ + d₃ = h', GOLD);
  } else if (id === 'triangle-euler-line') {
    const center = [-1.8, .5], r = 1.7, at = angle => add(center, [r * Math.cos(angle), r * Math.sin(angle)]);
    // Keep an obtuse triangle well away from the equilateral collapse O = G = H.
    const vertices = [-.2 + .08 * Math.sin(time * .35), 3.2, 4.6].map(at), {O, G, H} = triangleCenters(...vertices);
    circle(O, r, FAINT); polygon(vertices, CYAN);
    vertices.forEach((v, i) => {
      const mid = mul(add(vertices[(i + 1) % 3], vertices[(i + 2) % 3]), .5);
      line([v, mid], FAINT, 1);
    });
    line([O, H], GOLD, 2.6);
    [[O, 'O', GREEN, [.4, .2]], [G, 'G', GOLD, [.4, 0]], [H, 'H', CYAN, [-.15, -.75]]].forEach(([v, text, color, offset]) => {dot(v, color, 4); label(add(v, offset), text, color);});
    // The separate ratio ruler enlarges the center separations; the triangle retains its true scale.
    line([[1.1, .2], [3.8, .2]], GOLD, 2);
    [[1.1, 'O', GREEN], [2, 'G', GOLD], [3.8, 'H', CYAN]].forEach(([x, text, color]) => {dot([x, .2], color, 5); label([x - .12, -.65], text, color);});
    label([1.15, 1.7], 'OG : GH', GOLD); label([1.9, 1.05], '1 : 2', GOLD);
  } else if (id === 'parabola-focus') {
    const p = .55, offset = [-.75, 0], y = 1.9 * Math.sin(time * .38), result = parabolaFocusPoint(p, y);
    const [P, F, D] = [result.P, result.F, result.D].map(v => add(offset, v));
    curve(y => add(offset, parabolaFocusPoint(p, y).P), -2.15, 2.15, CYAN);
    line([[-p + offset[0], -2.25], [-p + offset[0], 2.25]], FAINT);
    line([P, F], GREEN, 2); line([P, D], GOLD, 2); dot(P, '#ffffff', 4); dot(F, GREEN, 4); dot(D, GOLD, 3);
    label(add(F, [.12, -.75]), 'F', GREEN); label(add(P, [.22, .4]), 'P', CYAN); label(add(D, [-.65, -.15]), 'D', GOLD);
    label([2.15, 1.25], 'PF', GREEN); label([2.15, -.15], 'PD', GOLD);
    line([[2.15, .7], [2.15 + result.focusDistance, .7]], GREEN, 6);
    line([[2.15, -.7], [2.15 + result.directrixDistance, -.7]], GOLD, 6);
    label([-2.5, -2.65], 'PF = PD', GOLD);
  } else if (id === 'hyperbola-foci') {
    const a = 1.25, b = .85, parameter = 1.2 * Math.sin(time * .35), result = hyperbolaFocusPoint(a, b, parameter);
    for (const branch of [-1, 1]) curve(t => hyperbolaFocusPoint(a, b, t, branch).P, -1.55, 1.55, FAINT);
    result.foci.forEach((F, i) => {const color = [GREEN, GOLD][i]; line([F, result.P], color, 2); dot(F, color, 4); label(add(F, [-.25, -.8]), `F${i + 1}`, color);});
    dot(result.P, CYAN, 5);
    const [d1, d2] = result.distances, scale = .9, x = -2.6;
    line([[x, -2.2], [x + d1 * scale, -2.2]], GREEN, 6);
    line([[x, -2.55], [x + d2 * scale, -2.55]], GOLD, 6);
    line([[x + d2 * scale, -2.55], [x + d1 * scale, -2.55]], CYAN, 3);
    label([-1.8, 2.55], '|PF₁ − PF₂| = 2a', CYAN);
  } else if (id === 'sine-rule') {
    const center = [-1.8, 0], radius = 1.8, angles = [.9 + .35 * Math.sin(time * .35), 2.8, 4.8];
    const vertices = angles.map(a => add(center, [radius * Math.cos(a), radius * Math.sin(a)]));
    const values = sineRuleValues(...vertices), colors = [GREEN, GOLD, CYAN];
    circle(center, radius, FAINT);
    vertices.forEach((v, i) => {
      const next = vertices[(i + 1) % 3], last = vertices[(i + 2) % 3];
      line([next, last], colors[i], 2); angleArc(v, next, last, .33, colors[i]);
      label(add(v, [[-.1, .45], [-.65, .35], [-.15, -.75]][i]), ['A', 'B', 'C'][i], colors[i]);
      label([.45, 1.1 - i * .85], ['a/sin A', 'b/sin B', 'c/sin C'][i], colors[i]);
      label([3.6, 1.1 - i * .85], values.ratios[i].toFixed(2), colors[i]);
    });
    dot(vertices[0], GREEN, 5); dot(center, FAINT, 3);
    label([1.65, -1.65], '= 2R', GOLD);
  }
  return true;
}
