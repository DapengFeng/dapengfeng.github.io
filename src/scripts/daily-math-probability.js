// Probability constructions: every displayed quantity comes from these pure models.
export const probabilityRenderers = ['law-large-numbers', 'chebyshev-bound', 'conditional-probability', 'covariance-correlation', 'mixture-density', 'uniform-maximum', 'inverse-cdf', 'kl-divergence'];
const GREEN = '#c5ec96', GOLD = '#e5bf82', CYAN = '#83bec8', FAINT = '#9ab79d45';
function finite(...values) { if (!values.every(Number.isFinite)) throw new RangeError('Finite values required'); }
function probability(p) { finite(p); if (p < 0 || p > 1) throw new RangeError('Probability must be in [0, 1]'); }
function count(n) { if (!Number.isSafeInteger(n) || n < 1) throw new RangeError('Positive integer count required'); }

// Seeded, finite Bernoulli illustrations; no claim that a finite path proves a limit.
export function bernoulliRunningMeans(n, p = .5, seed = 1) {
  count(n); probability(p);
  if (!Number.isInteger(seed) || seed < 0 || seed > 0xffffffff) throw new RangeError('Unsigned 32-bit seed required');
  let state = seed >>> 0, sum = 0;
  return Array.from({length: n}, (_, i) => {
    state = (Math.imul(1664525, state) + 1013904223) >>> 0;
    sum += state / 0x100000000 < p ? 1 : 0;
    return sum / (i + 1);
  });
}
export function bernoulliMeanVariance(n, p) { count(n); probability(p); return p * (1 - p) / n; }
export function chebyshevBound(variance, deviation) {
  finite(variance, deviation);
  if (variance < 0 || deviation <= 0) throw new RangeError('Nonnegative variance and positive deviation required');
  return Math.min(1, variance / deviation ** 2);
}
// Exact tail for the zero-mean, unit-variance uniform distribution.
export function standardizedUniformTail(k) {
  finite(k); if (k < 0) throw new RangeError('Nonnegative threshold required');
  return Math.max(0, 1 - k / Math.sqrt(3));
}
export function conditionalProbability(intersection, condition) {
  probability(intersection); probability(condition);
  if (condition === 0 || intersection > condition) throw new RangeError('Positive condition probability containing the intersection required');
  return intersection / condition;
}
// Population moments of equally weighted points (division by n, not n - 1).
export function covarianceStatistics(points) {
  if (!Array.isArray(points) || !points.length || points.some(p => !Array.isArray(p) || p.length !== 2 || !p.every(Number.isFinite))) throw new RangeError('Nonempty finite plane points required');
  const n = points.length, mean = [0, 1].map(axis => points.reduce((s, p) => s + p[axis], 0) / n);
  const centered = points.map(p => p.map((x, axis) => x - mean[axis]));
  const variance = [0, 1].map(axis => centered.reduce((s, p) => s + p[axis] ** 2, 0) / n);
  if (variance.some(v => !(v > 0))) throw new RangeError('Both variables must have positive variance');
  const covariance = centered.reduce((s, p) => s + p[0] * p[1], 0) / n;
  const standardDeviation = variance.map(Math.sqrt);
  return {mean, variance, covariance, correlation: covariance / (standardDeviation[0] * standardDeviation[1]), standardized: centered.map(p => p.map((x, axis) => x / standardDeviation[axis]))};
}
export function gaussianDensity(x, mean = 0, sigma = 1) {
  finite(x, mean, sigma); if (sigma <= 0) throw new RangeError('Positive standard deviation required');
  return Math.exp(-.5 * ((x - mean) / sigma) ** 2) / (sigma * Math.sqrt(2 * Math.PI));
}
export function gaussianMixtureDensity(x, weight = .5, separation = 1, sigma = .65) {
  probability(weight); finite(separation); if (separation < 0) throw new RangeError('Nonnegative half-separation required');
  return weight * gaussianDensity(x, -separation, sigma) + (1 - weight) * gaussianDensity(x, separation, sigma);
}
export function uniformMaximumCDF(x, n) {
  finite(x); count(n); return x <= 0 ? 0 : x >= 1 ? 1 : x ** n;
}
export function uniformMaximumDensity(x, n) {
  finite(x); count(n); return x < 0 || x > 1 ? 0 : n * x ** (n - 1);
}
export function uniformMaximumMean(n) { count(n); return n / (n + 1); }
export function triangularCDF(x) { finite(x); return x <= 0 ? 0 : x >= 1 ? 1 : x * x; }
export function triangularQuantile(u) { probability(u); return Math.sqrt(u); }
export function bernoulliKL(p, q) {
  probability(p); probability(q);
  if (p === q) return 0;
  if (q === 0 || q === 1) return Infinity;
  const first = p === 0 ? 0 : p * Math.log(p / q);
  const second = p === 1 ? 0 : (1 - p) * Math.log((1 - p) / (1 - q));
  return Math.max(0, first + second); // Avoid a negative rounding residue at equality.
}

const MEAN_PATHS = [11, 102, 3033].map(seed => bernoulliRunningMeans(256, .5, seed));
const CLOUD = Array.from({length: 32}, (_, i) => {
  const theta = 2 * Math.PI * i / 32;
  return [Math.SQRT2 * Math.cos(theta), Math.SQRT2 * (.65 * Math.cos(theta) + Math.sqrt(1 - .65 ** 2) * Math.sin(theta))];
});

export function drawProbabilityConstruction({id, time, phase, progress, line, dot, label, curve}) {
  if (!probabilityRenderers.includes(id)) return false;
  const rect = (left, right, bottom, top, color, fill) => line([[left, bottom], [right, bottom], [right, top], [left, top], [left, bottom]], color, 1.4, fill);
  if (id === 'law-large-numbers') {
    const x = n => -4.3 + (n - 1) * 8.6 / 255, y = p => -1.7 + 3.4 * p;
    const position = 1 + phase * 255, n = Math.min(256, Math.floor(position)), fraction = position - n;
    line([[x(1), y(.5)], [x(256), y(.5)]], GOLD, 1.5);
    MEAN_PATHS.forEach((path, i) => {
      const color = [GREEN, CYAN, '#b6a7d1'][i];
      line(path.map((mean, index) => [x(index + 1), y(mean)]), FAINT, 1);
      const value = path[n - 1] * (1 - fraction) + path[Math.min(n, 255)] * fraction;
      line([...path.slice(0, n).map((mean, index) => [x(index + 1), y(mean)]), [x(position), y(value)]], color, 1.5);
      dot([x(position), y(value)], color, 4);
    });
    label([-4.3, 2.45], 'p = 0.5', GOLD); label([1.35, 2.45], `n = ${n}`, GREEN);
    label([-4.3, -2.45], '1'); label([3.55, -2.45], '256');
  } else if (id === 'chebyshev-bound') {
    const k = .65 + 1.75 * progress, a = Math.sqrt(3), x = v => 1.7 * v;
    const actual = standardizedUniformTail(k), bound = chebyshevBound(1, k);
    rect(x(-a), x(a), .15, 1.65, CYAN, '#83bec814');
    if (k < a) { rect(x(-a), x(-k), .15, 1.65, GREEN, '#c5ec9650'); rect(x(k), x(a), .15, 1.65, GREEN, '#c5ec9650'); }
    for (const sign of [-1, 1]) line([[x(sign * k), -.05], [x(sign * k), 1.85]], GOLD, 1.6);
    label([-4.3, 2.5], 'X ~ U[−√3, √3]', CYAN); label([2.4, 2.5], `k=${k.toFixed(2)}`, GOLD);
    line([[-4.3, -.9], [4.3, -.9]], FAINT, 6); line([[-4.3, -.9], [-4.3 + 8.6 * bound, -.9]], GOLD, 6);
    line([[-4.3, -1.8], [4.3, -1.8]], FAINT, 6); line([[-4.3, -1.8], [-4.3 + 8.6 * actual, -1.8]], GREEN, 6);
    label([-4.3, -.48], `bound = ${bound.toFixed(3)}`, GOLD); label([-4.3, -2.45], `tail = ${actual.toFixed(3)}`, GREEN);
  } else if (id === 'conditional-probability') {
    const b = .05 + .75 * progress, x = u => -4.3 + 8.6 * u;
    const intersection = Math.max(0, .6 - b), p = conditionalProbability(intersection, 1 - b);
    rect(x(0), x(1), -.15, 1.65, FAINT, '#83bec80c');
    rect(x(0), x(.6), -.15, 1.65, CYAN, '#83bec828');
    if (b < .6) rect(x(b), x(.6), -.15, 1.65, GREEN, '#c5ec9650');
    rect(x(b), x(1), -.15, 1.65, GOLD);
    label([-4.3, 2.5], 'A: U≤0.6', CYAN); label([1.35, 2.5], `B: U≥${b.toFixed(2)}`, GOLD);
    line([[-4.3, -1.4], [4.3, -1.4]], GOLD, 6);
    line([[-4.3, -1.4], [-4.3 + 8.6 * p, -1.4]], GREEN, 8);
    label([-4.3, -2.45], `P(A|B) = ${p.toFixed(3)}`, GREEN); label([3.3, -2.45], '1', GOLD);
  } else if (id === 'covariance-correlation') {
    const scale = .55 + 1.1 * progress, points = CLOUD.map(([x, y]) => [scale * x, y]);
    const stats = covarianceStatistics(points), centers = [-2.35, 2.35];
    for (const center of centers) {line([[center - 1.75, 0], [center + 1.75, 0]], FAINT, 1); line([[center, -1.5], [center, 1.5]], FAINT, 1);}
    points.forEach((p, i) => { dot([centers[0] + .7 * p[0], .7 * p[1]], CYAN, 2.7); dot([centers[1] + .7 * stats.standardized[i][0], .7 * stats.standardized[i][1]], GREEN, 2.7); });
    const theta = time * .5, u = Math.SQRT2 * Math.cos(theta), v = Math.SQRT2 * (.65 * Math.cos(theta) + Math.sqrt(1 - .65 ** 2) * Math.sin(theta));
    dot([centers[0] + .7 * scale * u, .7 * v], GOLD, 4); dot([centers[1] + .7 * u, .7 * v], GOLD, 4);
    label([-4.3, 2.45], `a=${scale.toFixed(2)}`, CYAN); label([.7, 2.45], 'Zₓ, Zᵧ', GREEN);
    label([-4.3, -2.15], `Cov=${stats.covariance.toFixed(3)}`, CYAN); label([.7, -2.15], `ρ=${stats.correlation.toFixed(2)}`, GREEN);
  } else if (id === 'mixture-density') {
    const d = .2 + 1.55 * progress, sigma = .65, at = (x, density) => [x, -1.7 + 5.3 * density];
    line([[-4.3, -1.7], [4.3, -1.7]], FAINT, 1);
    curve(x => at(x, .5 * gaussianDensity(x, -d, sigma)), -4.3, 4.3, CYAN);
    curve(x => at(x, .5 * gaussianDensity(x, d, sigma)), -4.3, 4.3, GOLD);
    curve(x => at(x, gaussianMixtureDensity(x, .5, d, sigma)), -4.3, 4.3, GREEN);
    const x = -3.5 + 7 * phase; dot(at(x, gaussianMixtureDensity(x, .5, d, sigma)), GREEN, 4);
    label([-4.3, 2.5], '½N(−d, σ²)', CYAN); label([.7, 2.5], '½N(d, σ²)', GOLD);
    label([-4.3, -2.45], `d = ${d.toFixed(2)}`, GREEN); label([1.1, -2.45], 'σ = 0.65', GREEN);
  } else if (id === 'uniform-maximum') {
    const x = .04 + .92 * progress, at = (u, v) => [-4.2 + 8.4 * u, -1.65 + 3.4 * v];
    line([at(0, 0), at(1, 0)], FAINT, 1); line([at(0, 0), at(0, 1)], FAINT, 1);
    for (const n of [1, 2, 4, 8]) { const color = n === 4 ? GREEN : FAINT; curve(u => at(u, uniformMaximumCDF(u, n)), 0, 1, color); if (n !== 4) label(at(.8, uniformMaximumCDF(.8, n) + .035), `n=${n}`, CYAN); }
    line([at(x, 0), at(x, uniformMaximumCDF(x, 4))], GOLD, 1.5); dot(at(x, uniformMaximumCDF(x, 4)), GOLD, 5);
    label([-4.3, 2.5], 'P(Mₙ ≤ x) = xⁿ', GREEN);
    label([-4.3, -2.4], `x=${x.toFixed(2)}`, GOLD); label([.1, -2.4], `F₄=${(x ** 4).toFixed(3)}`, GREEN);
  } else if (id === 'inverse-cdf') {
    const u = .03 + .94 * progress, q = triangularQuantile(u), at = (x, y) => [-3.9 + 7.8 * x, -1.8 + 3.5 * y];
    line([at(0, 1), at(0, 0), at(1, 0)], FAINT, 1);
    for (let i = 1; i < 10; i++) {const v = i / 10, x = triangularQuantile(v); line([at(0, v), at(x, v), at(x, 0)], FAINT, 1);}
    curve(x => at(x, triangularCDF(x)), 0, 1, GREEN);
    line([at(0, u), at(q, u), at(q, 0)], GOLD, 2); dot(at(q, u), GOLD, 5); dot(at(q, 0), CYAN, 4);
    label([-4.3, 2.5], `u = ${u.toFixed(2)}`, GOLD); label([.8, 2.5], 'F(x) = x²', GREEN);
    label([-4.3, -2.5], '0'); label([.3, -2.5], `√u=${q.toFixed(3)}`, CYAN);
  } else if (id === 'kl-divergence') {
    const p = .3, q = .04 + .92 * progress, at = (q, value) => [-4.3 + 8.6 * q, -1.7 + 1.8 * value];
    line([at(0, 0), at(1, 0)], FAINT, 1);
    curve(q => at(q, bernoulliKL(p, q)), .04, .96, GREEN);
    curve(q => at(q, bernoulliKL(q, p)), .04, .96, CYAN);
    line([at(p, 0), at(p, 1.75)], FAINT, 1); dot(at(p, 0), GOLD, 3);
    for (const [value, color] of [[bernoulliKL(p, q), GREEN], [bernoulliKL(q, p), CYAN]]) dot(at(q, value), color, 4.5);
    label([-4.3, 2.5], 'p = 0.30', GOLD); label([1.2, 2.5], `q = ${q.toFixed(2)}`, GOLD);
    label([-4.3, -2.45], `p∥q:${bernoulliKL(p, q).toFixed(3)}`, GREEN); label([.3, -2.45], `q∥p:${bernoulliKL(q, p).toFixed(3)}`, CYAN);
  }
  return true;
}
