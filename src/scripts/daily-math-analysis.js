// Independently testable models and constructions for calculus and approximation.
// Coordinates belong to the shared mathematical canvas; no DOM or canvas state here.
export const analysisRenderers = ['fundamental-calculus', 'trapezoid-error', 'simpson-exactness', 'secant-root', 'lagrange-interpolation', 'chebyshev-nodes', 'osculating-circle', 'arc-length'];
const GREEN = '#c5ec96', GOLD = '#e5bf82', CYAN = '#83bec8', FAINT = '#9ab79d45';
function finite(...values) { if (!values.every(Number.isFinite)) throw new RangeError('Finite arguments required'); }
function count(n, minimum = 1) { if (!Number.isInteger(n) || n < minimum) throw new RangeError('Invalid subdivision count'); }

export function accumulatedArea(x) {
  finite(x);
  return {value: x / 2 + x ** 3 / 12, derivative: .5 + x * x / 4};
}
export function compositeTrapezoid(fn, a, b, n) {
  finite(a, b); count(n);
  if (b <= a) throw new RangeError('Increasing integration interval required');
  const h = (b - a) / n;
  const values = Array.from({length: n + 1}, (_, i) => fn(a + i * h));
  finite(...values);
  return h * (values.reduce((sum, y) => sum + y, 0) - (values[0] + values[n]) / 2);
}
export function simpsonPanel(fn, a, b) {
  finite(a, b);
  if (b <= a) throw new RangeError('Increasing integration interval required');
  const values = [fn(a), fn((a + b) / 2), fn(b)]; finite(...values);
  return (b - a) * (values[0] + 4 * values[1] + values[2]) / 6;
}
export function secantStep(fn, previous, current) {
  finite(previous, current);
  const f0 = fn(previous), f1 = fn(current); finite(f0, f1);
  if (f0 === f1) throw new RangeError('Secant slope must be nonzero');
  const next = current - f1 * (current - previous) / (f1 - f0); finite(next);
  return next;
}
export function secantSequence(fn, a, b, steps = 4) {
  count(steps, 0); const points = [a, b]; finite(a, b);
  for (let i = 0; i < steps; i++) points.push(secantStep(fn, points.at(-2), points.at(-1)));
  return points;
}
export function lagrangeValue(nodes, values, x) {
  finite(x);
  if (!Array.isArray(nodes) || !nodes.length || nodes.length !== values.length) throw new RangeError('Matching nonempty nodes and values required');
  finite(...nodes, ...values);
  if (new Set(nodes).size !== nodes.length) throw new RangeError('Distinct nodes required');
  return values.reduce((sum, value, k) => sum + value * nodes.reduce((basis, node, j) => j === k ? basis : basis * (x - node) / (nodes[k] - node), 1), 0);
}
// n is polynomial degree; second-kind Chebyshev nodes include both endpoints.
export function interpolationNodes(n, chebyshev = false) {
  count(n);
  return Array.from({length: n + 1}, (_, k) => chebyshev ? -Math.cos(k * Math.PI / n) : -1 + 2 * k / n);
}
export function rungeInterpolation(n, chebyshev = false) {
  const nodes = interpolationNodes(n, chebyshev), values = nodes.map(x => 1 / (1 + 25 * x * x));
  return {nodes, values, evaluate: x => lagrangeValue(nodes, values, x)};
}
// For y=x²: center = point + normal/curvature, with normal pointing upward.
export function parabolaCurvature(x) {
  finite(x);
  const curvature = 2 / (1 + 4 * x * x) ** 1.5;
  return {point: [x, x * x], center: [-4 * x ** 3, .5 + 3 * x * x], radius: 1 / curvature, curvature};
}
// Polygonal length of y=x²/2. The exact antiderivative integrates sqrt(1+x²).
export function parabolaArcLength(a, b, n) {
  finite(a, b); count(n);
  if (b <= a) throw new RangeError('Increasing interval required');
  const points = Array.from({length: n + 1}, (_, i) => {const x = a + (b - a) * i / n; return [x, x * x / 2];});
  const length = points.slice(1).reduce((sum, p, i) => sum + Math.hypot(p[0] - points[i][0], p[1] - points[i][1]), 0);
  const primitive = x => (x * Math.hypot(1, x) + Math.asinh(x)) / 2;
  return {points, length, exact: primitive(b) - primitive(a)};
}

const samples = (fn, a, b, n = 60) => Array.from({length: n + 1}, (_, i) => fn(a + (b - a) * i / n));
export function drawAnalysisConstruction({id, time, line, dot, label, curve, circle, progress}) {
  if (!analysisRenderers.includes(id)) return false;
  if (id === 'fundamental-calculus') {
    const x = .15 + 1.7 * progress, f = s => accumulatedArea(s).derivative;
    const left = (s, y) => [-4.5 + 1.9 * s, -1.65 + 1.8 * y], right = (s, y) => [.65 + 1.9 * s, -1.65 + 1.8 * y];
    for (const at of [left, right]) line([at(0, 0), at(2, 0)], FAINT, 1);
    line([left(0, 0), ...samples(s => left(s, f(s)), 0, x), left(x, 0)], FAINT, 1, '#c5ec9625');
    curve(s => left(s, f(s)), 0, 2, GREEN); line([left(x, 0), left(x, f(x))], GOLD, 1.5); dot(left(x, f(x)), GOLD, 4);
    curve(s => right(s, accumulatedArea(s).value), 0, 2, CYAN);
    const area = accumulatedArea(x);
    line([right(x - .14, area.value - .14 * f(x)), right(x + .14, area.value + .14 * f(x))], GOLD, 3);
    dot(right(x, area.value), GOLD, 4);
    label([-4.5, 2.4], 'f(x)', GREEN); label([.65, 2.4], 'A(x)', CYAN);
    label([-4.5, -2.35], `x = ${x.toFixed(2)}`, GOLD); label([.65, -2.35], `A′ = ${f(x).toFixed(2)}`, GOLD);
  } else if (id === 'trapezoid-error') {
    const n = [2, 4, 8, 16][Math.floor(time / 8) % 4], at = (x, y) => [-4.3 + 8.6 * x, -1.5 + 2.9 * y];
    line([at(0, 0), at(1, 0)], FAINT, 1);
    for (let i = 0; i < n; i++) {
      const a = i / n, b = (i + 1) / n;
      line([at(a, a * a), at(b, b * b), ...samples(x => at(x, x * x), b, a, 15)], GOLD, 1, '#e5bf8235');
      line([at(a, 0), at(a, a * a)], FAINT, 1);
    }
    curve(x => at(x, x * x), 0, 1, GREEN);
    const x = progress, a = Math.min(n - 1, Math.floor(x * n)) / n, b = a + 1 / n, chord = a * a + (a + b) * (x - a);
    line([at(x, x * x), at(x, chord)], CYAN, 3); dot(at(x, chord), GOLD, 4); dot(at(x, x * x), GREEN, 3);
    const error = compositeTrapezoid(x => x * x, 0, 1, n) - 1 / 3;
    label([-4.3, 2.4], `n = ${n}`, GOLD); label([.65, 2.4], 'I = 1/3', GREEN);
    label([-4.3, -2.35], `h = ${(1 / n).toFixed(4)}`, GOLD); label([.65, -2.35], `E = ${error.toFixed(5)}`, CYAN);
  } else if (id === 'simpson-exactness') {
    const f = x => x ** 3, q = x => 1.5 * x * x - .5 * x, at = (x, y) => [-4.3 + 8.6 * x, -1.25 + 2.7 * y];
    line([at(0, 0), at(1, 0)], FAINT, 1);
    for (const [a, b, color] of [[0, .5, '#c5ec9630'], [.5, 1, '#e5bf8230']]) line([...samples(x => at(x, f(x)), a, b), ...samples(x => at(x, q(x)), b, a)], FAINT, 1, color);
    curve(x => at(x, f(x)), 0, 1, GREEN); curve(x => at(x, q(x)), 0, 1, GOLD);
    for (const x of [0, .5, 1]) dot(at(x, f(x)), CYAN, 4);
    const x = progress; line([at(x, f(x)), at(x, q(x))], CYAN, 2); dot(at(x, f(x)), GREEN, 4);
    label([-4.3, 2.4], 'f = x³', GREEN); label([.3, 2.4], 'q₂(x)', GOLD);
    label([-4.3, -2.35], '∫f = ∫q = 1/4', CYAN);
  } else if (id === 'secant-root') {
    const f = x => x * x - 2, sequence = secantSequence(f, 1, 2, 4), step = Math.floor(time / 6) % 4;
    const a = sequence[step], b = sequence[step + 1], next = sequence[step + 2], lo0 = Math.min(a, b, next), hi0 = Math.max(a, b, next);
    const span = hi0 - lo0, lo = lo0 - .18 * span, hi = hi0 + .18 * span, magnitude = Math.max(Math.abs(f(lo)), Math.abs(f(hi)));
    const at = (x, y) => [-4.2 + 8.4 * (x - lo) / (hi - lo), 1.7 * y / magnitude];
    line([at(lo, 0), at(hi, 0)], FAINT, 1); curve(x => at(x, f(x)), lo, hi, GREEN);
    line([at(a, f(a)), at(b, f(b)), at(next, 0)], GOLD, 1.8);
    dot(at(a, f(a)), GOLD, 4); dot(at(b, f(b)), GOLD, 4); dot(at(next, 0), CYAN, 4);
    const fraction = (time % 6) / 6, x = b + fraction * (next - b), y = f(b) * (1 - fraction); dot(at(x, y), CYAN, 4);
    label([-4.3, 2.4], `k = ${step + 1}`, GOLD); label([-.8, 2.4], `x₊ = ${next.toFixed(6)}`, CYAN);
    label([-4.3, -2.35], 'x² − 2', GREEN); label([.65, -2.35], `W = ${(hi - lo).toFixed(4)}`, GOLD);
  } else if (id === 'lagrange-interpolation') {
    const nodes = [-1, 0, 1], values = [.2, .4 + .5 * Math.sin(time * .35), .4], f = x => lagrangeValue(nodes, values, x), at = (x, y) => [3.4 * x, -.5 + 1.9 * y];
    line([at(-1.25, 0), at(1.25, 0)], FAINT, 1);
    curve(x => at(x, f(x)), -1.25, 1.25, GREEN);
    line([at(0, -.1), at(0, .9)], FAINT, 1);
    nodes.forEach((x, i) => {dot(at(x, values[i]), i === 1 ? GOLD : CYAN, 5); label([3.4 * x - .15, -1.65], `${x}`, i === 1 ? GOLD : CYAN);});
    label([-4.3, 2.4], 'p(xᵢ) = yᵢ', GREEN); label([.65, 2.4], 'deg p ≤ 2', CYAN);
    label([-4.3, -2.35], `y₁ = ${values[1].toFixed(2)}`, GOLD);
  } else if (id === 'chebyshev-nodes') {
    const n = [6, 8, 10][Math.floor(time / 10) % 3], equal = rungeInterpolation(n), cheb = rungeInterpolation(n, true), f = x => 1 / (1 + 25 * x * x);
    const at = (x, y) => [4.3 * x, -.1 + .95 * y];
    curve(x => at(x, equal.evaluate(x)), -1, 1, GOLD, 360); curve(x => at(x, cheb.evaluate(x)), -1, 1, CYAN, 360); curve(x => at(x, f(x)), -1, 1, GREEN, 360);
    const x = -.98 + 1.96 * progress;
    for (const [model, color, y] of [[equal, GOLD, -1.95], [cheb, CYAN, -2.3]]) {
      line([at(x, f(x)), at(x, model.evaluate(x))], color, 1.5); dot(at(x, model.evaluate(x)), color, 4);
      line([[-4.3, y], [4.3, y]], FAINT, 1); model.nodes.forEach(node => dot([4.3 * node, y], color, 2.8));
    }
    dot(at(x, f(x)), GREEN, 3);
    label([-4.3, 2.4], `n = ${n}`, GREEN); label([-1.8, 2.4], 'equal', GOLD); label([1.5, 2.4], 'Chebyshev', CYAN);
  } else if (id === 'osculating-circle') {
    const x = .55 * Math.sin(time * .35), model = parabolaCurvature(x), scale = .95, at = ([x, y]) => [scale * x, scale * y - 1];
    curve(x => at([x, x * x]), -1.75, 1.75, GREEN);
    circle(at(model.center), scale * model.radius, GOLD); line([at(model.point), at(model.center)], CYAN, 1.5);
    const length = .6, tangent = [1 / Math.hypot(1, 2 * x), 2 * x / Math.hypot(1, 2 * x)];
    line([at(model.point.map((v, i) => v - length * tangent[i])), at(model.point.map((v, i) => v + length * tangent[i]))], FAINT, 1.5);
    dot(at(model.point), GREEN, 5); dot(at(model.center), CYAN, 3);
    label([-4.3, 2.4], 'y = x²', GREEN); label([1.8, 2.4], 'Rκ = 1', GOLD);
    label([-4.3, -2.35], `κ = ${model.curvature.toFixed(3)}`, GREEN); label([.65, -2.35], `R = ${model.radius.toFixed(3)}`, GOLD);
  } else if (id === 'arc-length') {
    const n = [2, 4, 8, 16, 32][Math.floor(time / 8) % 5], model = parabolaArcLength(-1.5, 1.5, n), at = ([x, y]) => [2.3 * x, 2.3 * y - 1.2];
    curve(x => at([x, x * x / 2]), -1.5, 1.5, GREEN); line(model.points.map(at), GOLD, 1.6); model.points.forEach(p => dot(at(p), GOLD, 2.6));
    const x = -1.5 + 3 * progress, index = Math.min(n - 1, Math.floor(progress * n)), a = model.points[index], b = model.points[index + 1], q = (x - a[0]) / (b[0] - a[0]);
    const chord = [x, a[1] * (1 - q) + b[1] * q], point = [x, x * x / 2];
    line([at(a), at(b)], CYAN, 3); line([at(chord), at(point)], FAINT, 1); dot(at(point), GREEN, 4);
    label([-4.3, 2.4], `n = ${n}`, GOLD); label([.65, 2.4], 'y = x²/2', GREEN);
    label([-4.3, -2.35], `Lₙ = ${model.length.toFixed(4)}`, GOLD); label([.65, -2.35], `L = ${model.exact.toFixed(4)}`, GREEN);
  }
  return true;
}
