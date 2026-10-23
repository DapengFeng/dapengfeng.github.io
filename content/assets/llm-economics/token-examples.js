/**
 * API-bill comparisons for controlled token mixes, never hardware benchmarks.
 *
 * prepareTokenExamples(completeParameters, metadata) returns:
 * - totalTokens: the entered input + output total, or a 100M reference total;
 * - ratios: 90/10, 50/50 and 10/90 input/output mixes at that same total;
 *   uncached / cached contain hit, miss and output {tokens, amount} plus total;
 *   cached and cacheSaving are null when the selected catalog has no cache rate;
 * - scaleMatrix: the same ratios at 10M, 100M and 1,000M total tokens, uncached;
 * - effectivePrices: month-one prices weighted by each stream's peak/off shares.
 *
 * Every comparison holds the provider, currency and stream-specific time shares
 * constant. Equal token totals do not represent equivalent tasks. None of these
 * examples uses local throughput, predicts local cost or estimates payback.
 */
import {billTokens} from './token-model.js';
import {customDefaults} from './defaults.js';

const REFERENCE_TOTAL = 100_000_000;
const SCALE_TOTALS = [10_000_000, 100_000_000, 1_000_000_000];
const RATIOS = [
 {id: 'input-heavy', inputShare: .9, label: {en: '90% input · 10% output', zh: '90% 输入 · 10% 输出'}},
 {id: 'balanced', inputShare: .5, label: {en: '50% input · 50% output', zh: '50% 输入 · 50% 输出'}},
 {id: 'output-heavy', inputShare: .1, label: {en: '10% input · 90% output', zh: '10% 输入 · 90% 输出'}},
];

function apiBill(inputTokens, outputTokens, parameters, cacheHit) {
 const p = {...parameters, apiCacheHit: cacheHit};
 const peak = billTokens(inputTokens * p.inputPeakShare, outputTokens * p.outputPeakShare, p, 'peak');
 const off = billTokens(inputTokens * (1 - p.inputPeakShare), outputTokens * (1 - p.outputPeakShare), p, 'off');
 const categories = Object.fromEntries(['hit', 'miss', 'output'].map(key => [key, {
  tokens: peak[key].tokens + off[key].tokens,
  amount: peak[key].amount + off[key].amount,
 }]));
 const total = categories.hit.amount + categories.miss.amount + categories.output.amount;
 if (!Number.isFinite(total)) throw new RangeError('Example API bill is not finite');
 return {cacheHit, ...categories, total};
}

export function prepareTokenExamples(parameters, metadata = {}) {
 if (parameters === null || typeof parameters !== 'object' || metadata === null || typeof metadata !== 'object') {
  throw new RangeError('Token parameters and metadata must be objects');
 }
 const enteredTotal = parameters.monthlyInputTokens + parameters.monthlyOutputTokens;
 if (!Number.isFinite(enteredTotal)) throw new RangeError('Monthly token total must be finite');
 const usesEnteredTotal = enteredTotal > 0;
 const totalTokens = usesEnteredTotal ? enteredTotal : REFERENCE_TOTAL;
 const cacheAvailable = metadata.pricing?.cacheSupported !== false;
 const sourceQuote = metadata.pricing;
 const source = {
  kind: sourceQuote ? 'catalog' : 'manual',
  label: sourceQuote?.sourceLabel ? {...sourceQuote.sourceLabel} : {
   en: sourceQuote ? 'Catalog reference prices' : 'Entered prices',
   zh: sourceQuote ? '目录参考价' : '所填单价',
  },
  provider: String(sourceQuote?.provider || metadata.provider || ''),
  providerId: String(sourceQuote?.providerId || metadata.providerId || ''),
  modelName: String(metadata.name || ''),
  modelId: String(sourceQuote?.modelId || ''),
  quoteId: sourceQuote?.id ?? null,
  url: sourceQuote?.sourceUrl ?? null,
  fetchedAt: sourceQuote?.fetchedAt ?? null,
  modified: Boolean(metadata.modified),
 };
 source.identity = {provider: source.provider, name: source.modelName};
 const ratios = RATIOS.map(ratio => {
  const inputTokens = totalTokens * ratio.inputShare;
  const outputTokens = totalTokens - inputTokens;
  const uncached = apiBill(inputTokens, outputTokens, parameters, 0);
  const cached = cacheAvailable ? apiBill(inputTokens, outputTokens, parameters, .8) : null;
  return {
   ...ratio, label: {...ratio.label}, inputTokens, outputTokens, uncached, cached,
   cacheSaving: cached ? uncached.total - cached.total : null,
  };
 });
 const scaleMatrix = {
  totals: [...SCALE_TOTALS],
  rows: RATIOS.map(ratio => ({
   ...ratio, label: {...ratio.label},
   cells: SCALE_TOTALS.map(total => {
    const inputTokens = total * ratio.inputShare, outputTokens = total - inputTokens;
    return {totalTokens: total, inputTokens, outputTokens, ...apiBill(inputTokens, outputTokens, parameters, 0)};
   }),
  })),
 };
 const weighted = (type, share) => parameters.prices.peak[type] * share + parameters.prices.off[type] * (1 - share);
 return {
  currency: String(metadata.currency || 'CURRENCY'), source,
  totalTokens, usesEnteredTotal, referenceTotalTokens: REFERENCE_TOTAL,
  cacheAvailable, cacheHitComparison: .8,
  effectivePrices: {
   hit: cacheAvailable ? weighted('hit', parameters.inputPeakShare) : null,
   miss: weighted('miss', parameters.inputPeakShare),
   output: weighted('output', parameters.outputPeakShare),
  },
  assumptions: {
   month: 1, priceUnit: 'per million tokens', totalIncludesInputAndOutput: true,
   inputPeakShare: parameters.inputPeakShare, outputPeakShare: parameters.outputPeakShare,
   equalTokensAreNotEquivalentTasks: true, apiOnly: true,
  },
  ratios, scaleMatrix,
 };
}

function deepFreeze(value) {
 if (value && typeof value === 'object') {
  for (const nested of Object.values(value)) deepFreeze(nested);
  Object.freeze(value);
 }
 return value;
}

// Deliberately model-free teaching rates, not a supplier quote or measured load.
const teachingData = prepareTokenExamples({
 ...customDefaults,
 monthlyInputTokens: 0, monthlyOutputTokens: 0,
 prices: {peak: {hit: .5, miss: 2, output: 8}, off: {hit: .5, miss: 2, output: 8}},
}, {currency: 'CNY'});
export const teachingExample = deepFreeze({
 ...teachingData,
 source: {
  ...teachingData.source,
  kind: 'teaching',
  label: {en: 'Illustrative rates, not a provider quote', zh: '教学假设单价，不是供应商报价'},
 },
});
