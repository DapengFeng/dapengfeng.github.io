import test from 'node:test';
import assert from 'node:assert/strict';
import {prepareTokenExamples, teachingExample} from '../content/assets/llm-economics/token-examples.js';
import {customDefaults} from '../content/assets/llm-economics/defaults.js';
import {billTokens} from '../content/assets/llm-economics/token-model.js';

const near = (actual, expected) => assert.ok(Math.abs(actual - expected) <= 1e-10 * Math.max(1, Math.abs(expected)), `${actual} != ${expected}`);
const params = overrides => ({
 ...customDefaults,
 prices: {peak: {hit: .5, miss: 2, output: 8}, off: {hit: .5, miss: 2, output: 8}},
 ...overrides,
});
const metadata = {currency: 'USD', name: 'Example model', pricing: {
 id: 'test-provider:test-model', provider: 'Example provider', providerId: 'test-provider',
 modelId: 'test-model', sourceUrl: 'https://models.dev/api.json',
 fetchedAt: '2026-10-20T00:00:00Z', cacheSupported: true,
 sourceLabel: {en: 'Test catalog', zh: '测试目录'},
}};

test('teaching mixes hold total tokens constant and produce the six stated bills', () => {
 assert.equal(teachingExample.currency, 'CNY');
 assert.equal(teachingExample.source.kind, 'teaching');
 assert.equal(teachingExample.totalTokens, 100e6);
 assert.equal(teachingExample.usesEnteredTotal, false);
 assert.deepEqual(teachingExample.ratios.map(r => [r.uncached.total, r.cached.total]), [[260, 152], [500, 440], [740, 728]]);
 assert.deepEqual(teachingExample.ratios.map(r => r.cacheSaving), [108, 60, 12]);
 for (const ratio of teachingExample.ratios) {
  assert.equal(ratio.inputTokens + ratio.outputTokens, 100e6);
  for (const bill of [ratio.uncached, ratio.cached]) {
   near(bill.hit.tokens + bill.miss.tokens, ratio.inputTokens);
   near(bill.output.tokens, ratio.outputTokens);
   near(bill.hit.amount + bill.miss.amount + bill.output.amount, bill.total);
  }
 }
 assert.ok(Object.isFrozen(teachingExample.ratios[0].cached.hit));
});

test('separate input/output time shares weight each rate and reconcile with billTokens', () => {
 const p = params({
  inputPeakShare: .25, outputPeakShare: .75,
  prices: {peak: {hit: .5, miss: 2, output: 8}, off: {hit: .1, miss: 1, output: 4}},
 });
 const result = prepareTokenExamples(p, metadata);
 near(result.effectivePrices.hit, .2);
 near(result.effectivePrices.miss, 1.25);
 near(result.effectivePrices.output, 7);
 for (const [i, expected] of [[182.5, 106.9], [412.5, 370.5], [642.5, 634.1]].entries()) {
  const r = result.ratios[i];
  near(r.uncached.total, expected[0]); near(r.cached.total, expected[1]);
  for (const bill of [r.uncached, r.cached]) {
   const direct = ['peak', 'off'].map(bucket => {
    const isPeak = bucket === 'peak';
    return billTokens(r.inputTokens * (isPeak ? p.inputPeakShare : 1 - p.inputPeakShare),
     r.outputTokens * (isPeak ? p.outputPeakShare : 1 - p.outputPeakShare),
     {...p, apiCacheHit: bill.cacheHit}, bucket);
   });
   for (const category of ['hit', 'miss', 'output']) {
    near(bill[category].amount, direct.reduce((sum, b) => sum + b[category].amount, 0));
   }
  }
 }
});

test('entered input plus output controls mix total, while volume matrix remains explicit', () => {
 const result = prepareTokenExamples(params({monthlyInputTokens: 240e6, monthlyOutputTokens: 60e6}), metadata);
 assert.equal(result.totalTokens, 300e6);
 assert.equal(result.usesEnteredTotal, true);
 near(result.ratios[0].uncached.total, 780);
 near(result.ratios[0].cached.total, 456);
 assert.deepEqual(result.scaleMatrix.totals, [10e6, 100e6, 1000e6]);
 assert.deepEqual(result.scaleMatrix.rows.map(r => r.cells.map(c => c.total)), [[26, 260, 2600], [50, 500, 5000], [74, 740, 7400]]);
 for (const row of result.scaleMatrix.rows) for (const cell of row.cells) {
  assert.equal(cell.cacheHit, 0);
  near(cell.inputTokens + cell.outputTokens, cell.totalTokens);
  near(cell.hit.amount + cell.miss.amount + cell.output.amount, cell.total);
 }
});

test('a catalog without cache-read prices yields no invented cache bill or saving', () => {
 const result = prepareTokenExamples(params(), {...metadata, pricing: {...metadata.pricing, cacheSupported: false}});
 assert.equal(result.cacheAvailable, false);
 assert.equal(result.effectivePrices.hit, null);
 for (const r of result.ratios) {
  assert.equal(r.cached, null);
  assert.equal(r.cacheSaving, null);
  assert.equal(r.uncached.hit.tokens, 0);
 }
 assert.deepEqual(result.ratios.map(r => r.uncached.total), [260, 500, 740]);
});

test('changing local hardware assumptions cannot change these API-only comparisons', () => {
 const base = prepareTokenExamples(params(), metadata);
 const changed = prepareTokenExamples(params({
  inputTokensPerSecond: 98765, outputTokensPerSecond: 12345, capex: 999999,
  monthlyOpex: 888888, availability: .01, headroom: .2, horizon: 120,
  apiAnnualChange: -.5, inputAnnualGrowth: 2, outputAnnualGrowth: -0.9,
 }), metadata);
 assert.deepEqual(changed, base);
 assert.equal(base.assumptions.apiOnly, true);
 assert.equal(base.assumptions.equalTokensAreNotEquivalentTasks, true);
 assert.equal(base.assumptions.month, 1);
});

test('inputs remain unchanged and provenance stays tied to the selected provider/model', () => {
 const p = params(), m = structuredClone(metadata), originalP = structuredClone(p), originalM = structuredClone(m);
 const result = prepareTokenExamples(p, m);
 assert.deepEqual(p, originalP); assert.deepEqual(m, originalM);
 assert.equal(result.source.providerId, 'test-provider');
 assert.equal(result.source.modelId, 'test-model');
 assert.equal(result.source.url, 'https://models.dev/api.json');
 assert.equal(result.source.fetchedAt, metadata.pricing.fetchedAt);
 assert.equal(result.currency, 'USD');
 result.source.label.en = 'Changed output';
 assert.equal(m.pricing.sourceLabel.en, 'Test catalog');
 const doubled = prepareTokenExamples(params({
  prices: {peak: {hit: 1, miss: 4, output: 16}, off: {hit: 1, miss: 4, output: 16}},
 }), {...m, modified: true});
 assert.equal(doubled.source.modified, true);
 near(doubled.ratios[0].uncached.total, result.ratios[0].uncached.total * 2);
});

test('invalid profiles and overflowing totals fail instead of displaying invalid costs', () => {
 for (const change of [{monthlyInputTokens: -1}, {monthlyOutputTokens: NaN}, {monthlyInputTokens: Number.MAX_VALUE, monthlyOutputTokens: Number.MAX_VALUE}, {inputPeakShare: 1.1}]) {
  assert.throws(() => prepareTokenExamples(params(change), metadata), RangeError);
 }
 assert.throws(() => prepareTokenExamples({prices: {}}, metadata), RangeError);
});
