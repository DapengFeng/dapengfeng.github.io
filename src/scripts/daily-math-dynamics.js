// Numerical models and drawings for probability, evolution and numerical error.
// The parent renderer owns the canvas transform and save/restore lifecycle.
export const dynamicsRenderers = ['logistic-growth', 'markov-equilibrium', 'binary-entropy', 'exponential-memoryless', 'jensen-convexity', 'euler-global-error', 'damped-energy', 'total-variance'];
const GREEN = '#c5ec96', GOLD = '#e5bf82', CYAN = '#83bec8', FAINT = '#9ab79d45';

function finite(...values) {
  if (values.some(value => !Number.isFinite(value))) throw new RangeError('Arguments must be finite');
}
function probability(value) {
  finite(value);
  if (value < 0 || value > 1) throw new RangeError('Probability must be in [0, 1]');
}

export function logisticGrowth(t, rate = 1, capacity = 1, initial = .08) {
  finite(t, rate, capacity, initial);
  if (t < 0 || rate < 0 || capacity <= 0 || initial < 0 || initial > capacity) throw new RangeError('Invalid logistic growth parameters');
  if (initial === 0) return 0;
  return capacity / (1 + (capacity / initial - 1) * Math.exp(-rate * t));
}

// a = P(0 to 1), b = P(1 to 0); return P(X_n = 0).
export function twoStateProbability(n, initial, a = .2, b = .3) {
  probability(initial); probability(a); probability(b);
  if (!Number.isInteger(n) || n < 0) throw new RangeError('Step must be a nonnegative integer');
  if (a + b === 0) return initial;
  const equilibrium = b / (a + b);
  return equilibrium + (initial - equilibrium) * (1 - a - b) ** n;
}

export function binaryEntropy(p) {
  probability(p);
  if (p === 0 || p === 1) return 0;
  return -p * Math.log2(p) - (1 - p) * Math.log2(1 - p);
}

export function exponentialSurvival(t, rate = .7) {
  finite(t, rate);
  if (t < 0 || rate <= 0) throw new RangeError('Time must be nonnegative and rate positive');
  return Math.exp(-rate * t);
}

export function jensenQuadratic(a, b, weight) {
  finite(a, b); probability(weight);
  const mean = (1 - weight) * a + weight * b;
  const transformedMean = mean * mean;
  const meanTransformed = (1 - weight) * a * a + weight * b * b;
  return {mean, transformedMean, meanTransformed, gap: weight * (1 - weight) * (a - b) ** 2};
}

// Explicit Euler applied to y'=-y, y(0)=1, over a fixed duration.
export function eulerDecay(steps, duration = 3) {
  finite(duration);
  if (!Number.isInteger(steps) || steps < 1 || duration <= 0) throw new RangeError('Positive step count and duration required');
  const h = duration / steps, points = [[0, 1]];
  for (let n = 1; n <= steps; n++) points.push([n * h, points[n - 1][1] * (1 - h)]);
  return {h, points, error: Math.abs(points.at(-1)[1] - Math.exp(-duration))};
}

// Underdamped oscillator: x''+2 gamma x'+omega0^2 x=0.
// Energy is exact; no averaging or exponential-energy approximation is used.
export function dampedState(t, {gamma = .16, omega0 = 1.7, x0 = 1, v0 = 0, mass = 1} = {}) {
  finite(t, gamma, omega0, x0, v0, mass);
  if (t < 0 || gamma < 0 || omega0 <= gamma || mass <= 0) throw new RangeError('An underdamped oscillator with nonnegative time is required');
  const omega = Math.sqrt(omega0 * omega0 - gamma * gamma);
  const b = (v0 + gamma * x0) / omega, envelope = Math.exp(-gamma * t);
  const shape = x0 * Math.cos(omega * t) + b * Math.sin(omega * t);
  const x = envelope * shape;
  const v = envelope * (-gamma * shape - x0 * omega * Math.sin(omega * t) + b * omega * Math.cos(omega * t));
  const energy = mass * (v * v + omega0 * omega0 * x * x) / 2;
  return {x, v, energy, energyRate: -2 * mass * gamma * v * v, amplitudeBound: Math.hypot(x0, b) * envelope};
}

// Population variances: each observation has equal weight, so group weights
// are proportional to their sizes. This is not the unbiased sample estimator.
export function varianceDecomposition(groups) {
  if (!Array.isArray(groups) || !groups.length || groups.some(group => !Array.isArray(group) || !group.length)) throw new RangeError('Nonempty groups required');
  const values = groups.flat(); finite(...values);
  const mean = values.reduce((sum, value) => sum + value, 0) / values.length;
  const groupStats = groups.map(group => {
    const groupMean = group.reduce((sum, value) => sum + value, 0) / group.length;
    return {mean: groupMean, variance: group.reduce((sum, value) => sum + (value - groupMean) ** 2, 0) / group.length, weight: group.length / values.length};
  });
  const within = groupStats.reduce((sum, group) => sum + group.weight * group.variance, 0);
  const between = groupStats.reduce((sum, group) => sum + group.weight * (group.mean - mean) ** 2, 0);
  const total = values.reduce((sum, value) => sum + (value - mean) ** 2, 0) / values.length;
  return {mean, within, between, total, groupStats};
}

export function drawDynamicsConstruction({id, time, line, dot, label, curve, phase, progress}) {
  if (!dynamicsRenderers.includes(id)) return false;
  if (id === 'logistic-growth') {
    const x = t => -4.3 + t * 8.6 / 8, y = n => -1.8 + n * 3.6, t = phase * 8;
    const n = logisticGrowth(t), slope = n * (1 - n) * 3.6 / (8.6 / 8);
    line([[-4.3, y(0)], [4.3, y(0)]], FAINT, 1);
    line([[-4.3, y(.5)], [4.3, y(.5)]], FAINT, 1);
    line([[-4.3, y(1)], [4.3, y(1)]], GOLD, 1);
    curve(s => [x(s), y(logisticGrowth(s))], 0, 8, FAINT);
    curve(s => [x(s), y(logisticGrowth(s))], 0, t, GREEN);
    const left = Math.max(-4.3, x(t) - .5), right = Math.min(4.3, x(t) + .5);
    line([[left, y(n) + slope * (left - x(t))], [right, y(n) + slope * (right - x(t))]], CYAN, 2);
    dot([x(t), y(n)], GREEN, 5);
    label([-4.3, 2.5], 'N/K'); label([3.25, 2.1], 'K', GOLD);
    label([-4.3, -2.45], `t = ${t.toFixed(1)}`, CYAN);
    label([1.3, -2.45], `N = ${n.toFixed(2)}`, GREEN);
  } else if (id === 'markov-equilibrium') {
    const n = Math.floor(time * 1.4) % 19, x = k => -4.3 + k * 8.6 / 18, y = p => -1.7 + p * 3.4;
    line([[-4.3, y(.6)], [4.3, y(.6)]], GOLD, 1.5);
    for (const [initial, color] of [[1, GREEN], [0, CYAN]]) {
      const all = Array.from({length: 19}, (_, k) => [x(k), y(twoStateProbability(k, initial))]);
      line(all, FAINT, 1); line(all.slice(0, n + 1), color, 1.6);
      for (const point of all.slice(0, n + 1)) dot(point, color, 2.8);
      dot(all[n], color, 4.3);
    }
    label([-4.3, 2.5], 'P(Xₙ = 0)'); label([2.7, 2.5], 'π₀ = 0.6', GOLD);
    label([-4.3, -2.45], `n = ${n}`, CYAN);
  } else if (id === 'binary-entropy') {
    const p = progress, x = v => -4 + 8 * v, y = v => -1.3 + 3 * v, h = binaryEntropy(p);
    line([[x(0), y(0)], [x(1), y(0)]], FAINT, 1);
    line([[x(.5), y(0)], [x(.5), y(1)]], FAINT, 1);
    curve(v => [x(v), y(binaryEntropy(v))], 0, 1, GREEN);
    line([[x(p), y(0)], [x(p), y(h)]], GOLD, 1); dot([x(p), y(h)], GOLD, 5);
    line([[-4, -2], [x(p), -2]], GREEN, 6); line([[x(p), -2], [4, -2]], CYAN, 6);
    label([-4, 2.45], `p = ${p.toFixed(2)}`, GREEN); label([1.1, 2.45], `H = ${h.toFixed(2)}`, GOLD);
    label([-4, -2.55], '0'); label([3.65, -2.55], '1'); label([-.45, 1.95], '1 bit', GOLD);
  } else if (id === 'exponential-memoryless') {
    const s = .4 + 2.1 * progress, x = t => -4.3 + t * 8.6 / 8, y = v => .2 + 1.55 * v;
    curve(t => [x(t), y(exponentialSurvival(t))], 0, 8, FAINT);
    curve(t => [x(t), y(exponentialSurvival(t))], s, 8, GREEN);
    for (const [t, color] of [[s, GOLD], [s + 2, CYAN]]) {
      line([[x(t), .2], [x(t), y(exponentialSurvival(t))]], color, 1);
      dot([x(t), y(exponentialSurvival(t))], color, 4);
    }
    const lower = (u, v) => [-4.3 + u * 8.6 / 5, -2.1 + 1.55 * v];
    const reference = Array.from({length: 161}, (_, i) => lower(i / 32, exponentialSurvival(i / 32)));
    line(reference, GOLD, 4);
    curve(u => lower(u, exponentialSurvival(s + u) / exponentialSurvival(s)), 0, 5, GREEN);
    const u = phase * 5; dot(lower(u, exponentialSurvival(u)), CYAN, 4);
    label([-4.3, 2.45], `s = ${s.toFixed(2)}`, GOLD); label([1.4, 2.45], 'S(s+2)/S(s)', CYAN);
    label([2.75, 1.95], exponentialSurvival(2).toFixed(3), CYAN);
    label([-4.3, -2.6], 'u = t − s', GREEN);
  } else if (id === 'jensen-convexity') {
    const a = -2, b = 1.5, weight = progress, result = jensenQuadratic(a, b, weight);
    const at = (x, value) => [1.8 * x, -1.7 + .7 * value];
    curve(x => at(x, x * x), -2.2, 2.2, GREEN);
    line([at(a, a * a), at(b, b * b)], GOLD, 1.8); dot(at(a, a * a), GOLD); dot(at(b, b * b), GOLD);
    line([at(result.mean, result.transformedMean), at(result.mean, result.meanTransformed)], CYAN, 3);
    dot(at(result.mean, result.transformedMean), GREEN, 4); dot(at(result.mean, result.meanTransformed), GOLD, 4);
    label([-4.3, 2.45], `w = ${weight.toFixed(2)}`, GOLD); label([1.25, 2.45], `gap = ${result.gap.toFixed(2)}`, CYAN);
    label([-4.3, -2.45], 'f(x) = x²', GREEN);
  } else if (id === 'euler-global-error') {
    const steps = [4, 8, 16, 32][Math.floor(time / 10) % 4], model = eulerDecay(steps), t = (time * .4) % 3;
    const at = (t, y) => [-4.3 + t * 8.6 / 3, -1.7 + 3.3 * y];
    curve(t => at(t, Math.exp(-t)), 0, 3, GOLD);
    line(model.points.map(([t, y]) => at(t, y)), GREEN, 1.8);
    model.points.forEach(([t, y]) => dot(at(t, y), GREEN, 2.3));
    const step = Math.min(steps - 1, Math.floor(t / model.h)), fraction = (t - step * model.h) / model.h;
    const approximate = model.points[step][1] * (1 - fraction) + model.points[step + 1][1] * fraction;
    line([at(t, approximate), at(t, Math.exp(-t))], CYAN, 2); dot(at(t, approximate), GREEN, 4.5); dot(at(t, Math.exp(-t)), GOLD, 3);
    label([-4.3, 2.45], `h = ${model.h.toFixed(4)}`, GREEN); label([1.4, 2.45], `E₃ = ${model.error.toFixed(4)}`, CYAN);
    label([-4.3, -2.45], 'y′ = −y', GOLD); label([3.1, -2.45], 't = 3', GOLD);
  } else if (id === 'damped-energy') {
    const t = phase * 18, x = t => -4.3 + t * 8.6 / 18, initial = dampedState(0).energy;
    line([[-4.3, 1.05], [4.3, 1.05]], FAINT, 1);
    curve(s => [x(s), 1.05 + .9 * dampedState(s).x], 0, 18, FAINT);
    for (const sign of [-1, 1]) curve(s => [x(s), 1.05 + sign * .9 * dampedState(s).amplitudeBound], 0, 18, GOLD);
    curve(s => [x(s), 1.05 + .9 * dampedState(s).x], 0, t, GREEN);
    curve(s => [x(s), -2.05 + 1.4 * dampedState(s).energy / initial], 0, 18, FAINT);
    curve(s => [x(s), -2.05 + 1.4 * dampedState(s).energy / initial], 0, t, CYAN);
    const state = dampedState(t);
    dot([x(t), 1.05 + .9 * state.x], GREEN, 4); dot([x(t), -2.05 + 1.4 * state.energy / initial], CYAN, 4);
    label([-4.3, 2.45], 'x(t)', GREEN); label([-4.3, -.35], 'E/E₀', CYAN); label([2.7, -2.6], `t = ${t.toFixed(1)}`, GOLD);
  } else if (id === 'total-variance') {
    const d = .2 + 1.4 * progress, offsets = [-.7, -.35, 0, .35, .7];
    const groups = [offsets.map(value => value - d), offsets.map(value => value + d)], result = varianceDecomposition(groups);
    line([[0, -.35], [0, 1.8]], FAINT, 1);
    groups.forEach((group, i) => {
      const y = 1.3 - i * 1.1, color = i ? GREEN : CYAN, mean = result.groupStats[i].mean;
      line([[group[0] * 1.8, y], [group.at(-1) * 1.8, y]], FAINT, 1);
      group.forEach(value => dot([value * 1.8, y], color, 4));
      line([[mean * 1.8, y - .18], [mean * 1.8, y + .18]], GOLD, 2.5);
      label([-4.8, y + .1], i ? 'B' : 'A', color);
    });
    const start = -4.3, unit = 8.6 / 2.805, middle = start + result.within * unit;
    line([[start, -1.45], [middle, -1.45]], CYAN, 8);
    line([[middle, -1.45], [middle + result.between * unit, -1.45]], GREEN, 8);
    label([-4.3, 2.45], `Var = ${result.total.toFixed(3)}`, GOLD);
    label([-4.3, -2.25], `W = ${result.within.toFixed(3)}`, CYAN); label([1.1, -2.25], `B = ${result.between.toFixed(3)}`, GREEN);
  }
  return true;
}
