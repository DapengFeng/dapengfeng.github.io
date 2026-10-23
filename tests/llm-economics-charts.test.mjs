import test from 'node:test';
import assert from 'node:assert/strict';
import {evaluateTokenScenario} from '../content/assets/llm-economics/token-model.js';
import {prepareChartData} from '../content/assets/llm-economics/charts.js';

function parameters(overrides = {}) {
  return {
    monthlyInputTokens: 100e6, monthlyOutputTokens: 25e6,
    inputTokensPerSecond: 200, outputTokensPerSecond: 50, capacityMultiplier: 1,
    inputPeakShare: .5, outputPeakShare: .5, inputAnnualGrowth: 0, outputAnnualGrowth: 0,
    opexAnnualGrowth: 0, perMillionLocalInputCost: .05, perMillionLocalOutputCost: .1,
    apiCacheHit: .5, fallbackCacheHit: .25, peakHours: 360, offHours: 360,
    availability: .95, headroom: .8, capex: 1000, monthlyOpex: 20,
    horizon: 12, annualDiscount: .1, terminalResidual: 100, apiAnnualChange: 0,
    prices: {peak: {hit: .1, miss: 1, output: 4}, off: {hit: .05, miss: .5, output: 2}},
    ...overrides,
  };
}
const close = (actual, expected) => assert.ok(Math.abs(actual - expected) <= Math.max(1, Math.abs(expected)) * 1e-10, `${actual} != ${expected}`);

test('zero demand retains every entered price and selected model, with no invented sensitivity', () => {
  const r = evaluateTokenScenario(parameters({monthlyInputTokens: 0, monthlyOutputTokens: 0}));
  const d = prepareChartData(r, {name: 'Model A', version: 'rev-7', currency: 'USD'});
  assert.equal(d.hasDemand, false); assert.deepEqual(d.sensitivity, []);
  assert.equal(d.model, 'Model A'); assert.equal(d.version, 'rev-7'); assert.equal(d.currency, 'USD');
  assert.deepEqual(d.prices, [{key: 'hit', peak: .1, off: .05}, {key: 'miss', peak: 1, off: .5}, {key: 'output', peak: 4, off: 2}]);
  assert.equal(d.flatTariff, false);
});

test('model changes replace displayed tariff data, including flat and zero prices', () => {
  const a = prepareChartData(evaluateTokenScenario(parameters()), {name: 'A', currency: 'CNY'});
  const flat = {hit: .2, miss: .2, output: 1};
  const b = prepareChartData(evaluateTokenScenario(parameters({prices: {peak: flat, off: flat}})), {name: '<img src=x>', currency: 'USD'});
  assert.notDeepEqual(a.prices, b.prices); assert.equal(b.flatTariff, true);
  assert.equal(b.model, '<img src=x>'); // Renderer must keep identity as text, never markup.
  assert.equal(b.currency, 'USD');
  const zero = {hit: 0, miss: 0, output: 0};
  const z = prepareChartData(evaluateTokenScenario(parameters({prices: {peak: zero, off: zero}})));
  assert.ok(z.prices.every(p => p.peak === 0 && p.off === 0));
});

test('the same model keeps distinct provider identities and flags edited quote inputs', () => {
  const result = evaluateTokenScenario(parameters());
  const a = prepareChartData(result, {name: 'Shared model', pricing: {provider: 'Provider A', providerId: 'provider-a'}});
  const b = prepareChartData(result, {name: 'Shared model', pricing: {provider: 'Provider B', providerId: 'provider-b'}, modified: true});
  assert.equal(a.model, b.model);
  assert.equal(a.provider, 'Provider A'); assert.equal(a.providerId, 'provider-a'); assert.equal(a.modified, false);
  assert.equal(b.provider, 'Provider B'); assert.equal(b.providerId, 'provider-b'); assert.equal(b.modified, true);
  const custom = prepareChartData(result, {name: 'Custom', pricing: null});
  assert.equal(custom.provider, ''); assert.equal(custom.providerId, '');
});

test('an unlisted cache-read tariff is marked as unavailable, while retaining the calculation fallback', () => {
  const rates = {hit: .5, miss: .5, output: 2};
  const result = evaluateTokenScenario(parameters({apiCacheHit: 0, fallbackCacheHit: 0, prices: {peak: rates, off: rates}}));
  const data = prepareChartData(result, {pricing: {provider: 'No-cache provider', cacheSupported: false}});
  assert.equal(data.cachePriceListed, false);
  assert.equal(data.prices.find(row => row.key === 'hit').peak, rates.miss);
  assert.equal(prepareChartData(result, {pricing: {cacheSupported: true}}).cachePriceListed, true);
  assert.equal(prepareChartData(result, {pricing: null}).cachePriceListed, true);
});

test('cost composition is month-one running cash only and reconciles both routes', () => {
  const r = evaluateTokenScenario(parameters()), d = prepareChartData(r);
  close(d.cost.api, r.firstMonth.fullApiCost);
  close(d.cost.fixed + d.cost.variable + d.cost.fallback, d.cost.localTotal);
  close(d.cost.localTotal, r.firstMonth.hybridOperatingCost);
  close(d.cost.saving, d.cost.api - d.cost.localTotal);
  const differentCapex = prepareChartData(evaluateTokenScenario(parameters({capex: 1e6, terminalResidual: 5000})));
  assert.deepEqual(d.cost, differentCapex.cost);
});

test('cash curves exclude resale while terminal markers reconcile to net cash and NPV', () => {
  const r = evaluateTokenScenario(parameters({terminalResidual: 1e5})), d = prepareChartData(r);
  assert.equal(d.cash[0].ordinary, -r.parameters.capex); assert.equal(d.cash[0].discounted, -r.parameters.capex);
  const final = d.cash.at(-1);
  close(final.ordinary + d.terminal.amount, d.terminal.ordinary);
  close(final.discounted + d.terminal.presentValue, d.terminal.discounted);
  close(d.terminal.discounted, r.NPV);
  const noTerminal = evaluateTokenScenario(parameters({terminalResidual: 0}));
  assert.equal(r.simplePaybackMonth, noTerminal.simplePaybackMonth);
  assert.deepEqual(d.cash, prepareChartData(noTerminal).cash);
});

test('sensitivity scales input and output together, holds all other inputs, and reproduces baseline exactly', () => {
  const p = parameters({inputAnnualGrowth: .2, outputAnnualGrowth: -.1}), r = evaluateTokenScenario(p), d = prepareChartData(r);
  assert.equal(d.sensitivity[0].multiplier, .25); assert.equal(d.sensitivity.at(-1).multiplier, 3);
  assert.equal(d.sensitivity.find(row => row.multiplier === 1).NPV, r.NPV);
  for (const point of [d.sensitivity[0], d.sensitivity[10], d.sensitivity.at(-1)]) {
    close(point.monthlyInputTokens, p.monthlyInputTokens * point.multiplier);
    close(point.monthlyOutputTokens, p.monthlyOutputTokens * point.multiplier);
    const expected = evaluateTokenScenario({...p, monthlyInputTokens: point.monthlyInputTokens, monthlyOutputTokens: point.monthlyOutputTokens});
    close(point.NPV, expected.NPV);
  }
  assert.deepEqual(d.parameters, r.parameters);
});

test('a new valid input produces new cash and sensitivity rather than stale plot values', () => {
  const a = prepareChartData(evaluateTokenScenario(parameters()));
  const b = prepareChartData(evaluateTokenScenario(parameters({monthlyOutputTokens: 50e6, inputTokensPerSecond: 50, horizon: 24, annualDiscount: .2})));
  assert.equal(b.cash.length, 25); assert.notDeepEqual(a.cost, b.cost);
  assert.notEqual(a.sensitivity.find(row => row.multiplier === 1).NPV, b.sensitivity.find(row => row.multiplier === 1).NPV);
  assert.equal(b.terminal.month, 24);
});
