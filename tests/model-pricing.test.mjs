import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { PRICE_SOURCE, parseModelsDev, updateCatalog, validateCatalog } from '../scripts/lib/model-pricing.mjs';

const stamp = '2026-10-20T08:40:00.000Z';
const later = '2026-10-21T08:40:00.000Z';
const fixture = JSON.parse(await readFile(new URL('fixtures/model-pricing.json', import.meta.url), 'utf8'));
const { catalog } = await updateCatalog(null, async () => fixture, stamp);
const parsed = mutation => {
  const data = structuredClone(fixture);
  mutation?.(data.first.models['shared/model'], data);
  return parseModelsDev(data, stamp);
};

test('CC Switch upstream uses USD per million tokens without multiplying or relabelling currency', () => {
  assert.equal(PRICE_SOURCE.url, 'https://models.dev/api.json');
  assert.match(PRICE_SOURCE.referenceUrl, /farion1231\/cc-switch/);
  assert.equal(validateCatalog(catalog), catalog);
  assert.equal(catalog.schemaVersion, 2);
  assert.equal(catalog.quotes.length, 4);
  const q = catalog.quotes.find(q => q.providerId === 'first' && q.modelId === 'shared/model');
  assert.deepEqual(q.prices.peak, { hit: 0.1, miss: 1, output: 4 });
  assert.deepEqual(q.prices.off, q.prices.peak);
  assert.equal(q.currency, 'USD');
  assert.equal(q.sourceUrl, 'https://models.dev/api.json');
  assert.match(catalog.sources[0].conditionsDefault.en, /not an official quotation/);
  assert.match(catalog.sources[0].conditionsDefault.en, /time-of-day/);
});

test('identical model ids retain each provider price and provenance', () => {
  const quotes = catalog.quotes.filter(q => q.modelId === 'shared/model');
  assert.equal(quotes.length, 2);
  assert.notEqual(quotes[0].id, quotes[1].id);
  assert.deepEqual(quotes.map(q => q.prices.peak.miss).sort(), [1, 2]);
  const q = quotes.find(q => q.providerId === 'first');
  assert.equal(q.canonicalModelId, 'original/shared');
  assert.equal(q.upstreamUpdatedAt, '2026-10-20');
  assert.equal(q.contextLimit, 131072);
  assert.equal(q.providerDocs, 'https://example.com/first/pricing');
});

test('unknown cache price is not free; explicit zero is a valid listed rate', () => {
  const missing = catalog.quotes.find(q => q.modelId === 'no-cache');
  assert.equal(missing.cacheSupported, false);
  assert.equal(missing.prices.peak.hit, 0.2);
  const free = catalog.quotes.find(q => q.modelId === 'free');
  assert.equal(free.cacheSupported, true);
  assert.deepEqual(free.prices.peak, { hit: 0, miss: 0, output: 0 });
});

test('missing input/output, malformed prices and unsupported billing dimensions are excluded', () => {
  const cases = [
    ['missingOrInvalidPrice', q => { delete q.cost.input; }],
    ['missingOrInvalidPrice', q => { q.cost.output = null; }],
    ['missingOrInvalidPrice', q => { q.cost.input = '1'; }],
    ['missingOrInvalidPrice', q => { q.cost.output = -1; }],
    ['invalidCachePrice', q => { q.cost.cache_read = Infinity; }],
    ['cacheWritePricing', q => { q.cost.cache_write = 1.25; }],
    ['separateReasoningPricing', q => { q.cost.reasoning = 8; }],
    ['separateAudioPricing', q => { q.cost.input_audio = 0.1; }],
    ['contextTierPricing', q => { q.cost.tiers = [{ input: 2, output: 8, tier: { type: 'context', size: 200000 } }]; }],
    ['contextTierPricing', q => { q.cost.context_over_200k = { input: 2, output: 8 }; }],
    ['unknownPriceDimension', q => { q.cost.cache_storage = 0.01; }],
  ];
  for (const [reason, mutate] of cases) {
    const { quotes, coverage } = parsed(mutate);
    assert.equal(quotes.some(q => q.providerId === 'first' && q.modelId === 'shared/model'), false, reason);
    assert.equal(coverage.excludedCounts[reason], 1, reason);
  }
});

test('only upstream-tagged open weights and non-deprecated text outputs enter the catalogue', () => {
  assert.equal(catalog.coverage.excludedCounts.notOpenWeights, 1);
  for (const [reason, mutate] of [
    ['notOpenWeights', q => { delete q.open_weights; }],
    ['deprecated', q => { q.status = 'deprecated'; }],
    ['notTextOutput', q => { q.modalities.output = ['text', 'audio']; }],
    ['notTextOutput', q => { q.type = 'embedding'; }],
  ]) assert.ok(parsed(mutate).coverage.excludedCounts[reason] >= 1);
});

test('network or malformed data retains validated prior quotes and successful retrieval timestamp', async () => {
  for (const fetchData of [async () => { throw new Error('offline'); }, async () => ({ wrong: {} })]) {
    const { catalog: retained, failures } = await updateCatalog(catalog, fetchData, later);
    assert.equal(failures.length, 1);
    assert.equal(retained.sources[0].status, 'retained');
    assert.equal(retained.sources[0].fetchedAt, stamp);
    assert.equal(retained.sources[0].lastAttemptAt, later);
    assert.deepEqual(retained.quotes, catalog.quotes);
    assert.deepEqual(retained.coverage, catalog.coverage);
  }
  await assert.rejects(updateCatalog(null, async () => { throw new Error('offline'); }, later), /no validated snapshot/);
});

test('an unexpectedly truncated catalogue cannot replace a broad snapshot', async () => {
  const full = structuredClone(fixture);
  for (let n = 0; n < 120; n++) full.first.models[`extra-${n}`] = structuredClone(full.first.models['no-cache']);
  const { catalog: broad } = await updateCatalog(null, async () => full, stamp);
  const { catalog: retained, failures } = await updateCatalog(broad, async () => fixture, later);
  assert.match(failures[0].error, /half/);
  assert.equal(retained.quotes.length, broad.quotes.length);
});

test('snapshot schema rejects spoofed sources, currency changes, duplicate ids and unsupported costs', async () => {
  for (const mutate of [
    c => { c.schemaVersion = 1; },
    c => { c.sources[0].url = 'https://example.com'; },
    c => { c.quotes[0].currency = 'CNY'; },
    c => { c.quotes[0].prices.peak.miss = NaN; },
    c => { c.quotes[0].prices.off.output += 1; },
    c => { c.quotes[0].fetchedAt = later; },
    c => { c.quotes.push(c.quotes[0]); },
    c => { c.quotes[0].upstreamCost.cache_write = 3; },
    c => { c.coverage.importedQuotes += 1; },
  ]) {
    const broken = structuredClone(catalog); mutate(broken);
    assert.throws(() => validateCatalog(broken));
    await assert.rejects(updateCatalog(broken, async () => fixture, later));
  }
});
