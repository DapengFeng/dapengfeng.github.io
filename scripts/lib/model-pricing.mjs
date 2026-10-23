// CC Switch's modelsDev.ts uses this same public catalogue. The catalogue is
// community-maintained; these rates are not fetched directly from providers.
export const PRICE_SOURCE = Object.freeze({
  id: 'models-dev',
  url: 'https://models.dev/api.json',
  referenceUrl: 'https://github.com/farion1231/cc-switch/blob/main/src/lib/modelsDev.ts',
  unitsUrl: 'https://github.com/anomalyco/models.dev#schema',
  label: { en: 'CC Switch upstream · models.dev', zh: 'CC Switch 同源 · models.dev' },
  tierDefault: { en: 'Catalogue base rate · USD / million tokens', zh: '目录单档价 · 美元／百万 token' },
  conditionsDefault: {
    en: 'Community catalogue, not an official quotation. One listed rate is applied to both periods; time-of-day pricing and contract/access/tool charges are not encoded. Verify the provider tariff and weight licence before a decision.',
    zh: '社区目录，不是官方报价。上游单档价暂用于两个时段，未编码峰谷、合同、接入及工具费用。决策前核对供应商价目与权重许可证。',
  },
  conditionsNoCache: {
    en: 'Community catalogue, not an official quotation. No cache-read rate is listed: use zero cache hits. One rate is used in both periods; time-of-day/contract/access/tool charges are not encoded. Verify tariff and licence.',
    zh: '社区目录，不是官方报价。未列缓存读取价，缓存命中比例设为零。单档价暂用于两个时段，未编码峰谷、合同、接入及工具费用；决策前核对价目与许可证。',
  },
});

const finitePrice = value => typeof value === 'number' && Number.isFinite(value) && value >= 0 && value < 1e12;
const isObject = value => value !== null && typeof value === 'object' && !Array.isArray(value);
const nonempty = value => typeof value === 'string' && value.trim().length > 0 && value.length <= 1024;
const knownCostKeys = new Set(['input', 'output', 'cache_read', 'cache_write', 'reasoning', 'input_audio', 'output_audio', 'tiers', 'context_over_200k']);
const quoteId = (providerId, modelId) => `models-dev:${encodeURIComponent(providerId)}:${encodeURIComponent(modelId)}`;

function requireCondition(condition, message) {
  if (!condition) throw new Error(message);
}
function timestamp(value) {
  return typeof value === 'string' && /^\d{4}-\d{2}-\d{2}T/.test(value) && Number.isFinite(Date.parse(value));
}
function bilingual(value) {
  return value && ['en', 'zh'].every(language => nonempty(value[language]));
}
function safeHttps(value) {
  if (!nonempty(value)) return null;
  try {
    const url = new URL(value);
    return url.protocol === 'https:' && !url.username && !url.password ? url.href : null;
  } catch { return null; }
}
function dateOnly(value) {
  return typeof value === 'string' && /^\d{4}-\d{2}(?:-\d{2})?$/.test(value) ? value : null;
}

export function parseModelsDev(data, fetchedAt) {
  requireCondition(isObject(data) && timestamp(fetchedAt), 'Invalid models.dev response or retrieval time');
  const quotes = [];
  const excludedCounts = {};
  let scannedModels = 0;
  const exclude = reason => { excludedCounts[reason] = (excludedCounts[reason] ?? 0) + 1; };
  for (const [providerId, provider] of Object.entries(data)) {
    requireCondition(nonempty(providerId) && isObject(provider) && isObject(provider.models), 'Invalid provider structure');
    for (const [modelId, model] of Object.entries(provider.models)) {
      scannedModels++;
      requireCondition(scannedModels <= 50_000, 'Catalogue model count exceeds limit');
      if (!nonempty(modelId) || !isObject(model)) { exclude('invalidModel'); continue; }
      // "Open weights" is an upstream tag, not a licence audit by this site.
      if (model.open_weights !== true) { exclude('notOpenWeights'); continue; }
      if (model.status?.toLowerCase?.() === 'deprecated') { exclude('deprecated'); continue; }
      const output = model.modalities?.output;
      if (!Array.isArray(output) || output.length !== 1 || output[0] !== 'text' || model.type === 'embedding') { exclude('notTextOutput'); continue; }
      const cost = model.cost;
      if (!isObject(cost) || !finitePrice(cost.input) || !finitePrice(cost.output)) { exclude('missingOrInvalidPrice'); continue; }
      if (Object.keys(cost).some(key => !knownCostKeys.has(key))) { exclude('unknownPriceDimension'); continue; }
      if (cost.tiers !== undefined || cost.context_over_200k !== undefined) { exclude('contextTierPricing'); continue; }
      if (cost.cache_write !== undefined && (!finitePrice(cost.cache_write) || cost.cache_write > 0)) { exclude('cacheWritePricing'); continue; }
      if (cost.input_audio !== undefined || cost.output_audio !== undefined) { exclude('separateAudioPricing'); continue; }
      if (cost.reasoning !== undefined && (!finitePrice(cost.reasoning) || cost.reasoning !== cost.output)) { exclude('separateReasoningPricing'); continue; }
      if (cost.cache_read !== undefined && !finitePrice(cost.cache_read)) { exclude('invalidCachePrice'); continue; }
      const cacheSupported = cost.cache_read !== undefined;
      const rates = { hit: cacheSupported ? cost.cache_read : cost.input, miss: cost.input, output: cost.output };
      const contextLimit = Number.isSafeInteger(model.limit?.context) && model.limit.context > 0 ? model.limit.context : null;
      quotes.push({
        id: quoteId(providerId, modelId),
        sourceId: PRICE_SOURCE.id,
        sourceUrl: PRICE_SOURCE.url,
        providerId,
        provider: nonempty(provider.name) ? provider.name : providerId,
        providerDocs: safeHttps(provider.doc),
        modelId,
        modelName: nonempty(model.name) ? model.name : modelId,
        revision: modelId,
        canonicalModelId: nonempty(model.canonical_model_id) ? model.canonical_model_id : null,
        openWeights: true,
        contextLimit,
        releaseDate: dateOnly(model.release_date),
        upstreamUpdatedAt: dateOnly(model.last_updated),
        currency: 'USD',
        billing: 'flat',
        prices: { peak: { ...rates }, off: { ...rates } },
        cacheSupported,
        upstreamCost: { ...cost },
        fetchedAt,
      });
    }
  }
  requireCondition(quotes.length > 0, 'No supported open-weight token quotes in models.dev response');
  quotes.sort((a, b) => a.provider.localeCompare(b.provider, 'en') || (b.releaseDate ?? '').localeCompare(a.releaseDate ?? '') || a.modelName.localeCompare(b.modelName, 'en') || a.id.localeCompare(b.id, 'en'));
  return {
    quotes,
    coverage: {
      scannedProviders: Object.keys(data).length,
      scannedModels,
      importedQuotes: quotes.length,
      importedProviders: new Set(quotes.map(q => q.providerId)).size,
      excludedCounts,
    },
  };
}

export function validateCatalog(catalog) {
  requireCondition(catalog?.schemaVersion === 2 && timestamp(catalog.generatedAt), 'Invalid model catalogue metadata');
  requireCondition(Array.isArray(catalog.sources) && catalog.sources.length === 1 && Array.isArray(catalog.quotes) && catalog.quotes.length > 0, 'Invalid model catalogue arrays');
  const source = catalog.sources[0];
  requireCondition(source.id === PRICE_SOURCE.id && source.url === PRICE_SOURCE.url && source.referenceUrl === PRICE_SOURCE.referenceUrl, 'Unexpected price source');
  requireCondition(['label', 'tierDefault', 'conditionsDefault', 'conditionsNoCache'].every(key => bilingual(source[key])) && ['ok', 'retained'].includes(source.status) && timestamp(source.fetchedAt) && timestamp(source.lastAttemptAt), 'Invalid source metadata');
  requireCondition(Date.parse(source.lastAttemptAt) >= Date.parse(source.fetchedAt) && Date.parse(catalog.generatedAt) >= Date.parse(source.lastAttemptAt), 'Invalid source chronology');
  const ids = new Set();
  for (const q of catalog.quotes) {
    requireCondition(q.sourceId === source.id && q.sourceUrl === source.url && q.currency === 'USD', 'Quote source or currency mismatch');
    for (const field of ['providerId', 'provider', 'modelId', 'modelName', 'revision']) requireCondition(nonempty(q[field]), `Invalid quote ${field}`);
    requireCondition(q.id === quoteId(q.providerId, q.modelId) && !ids.has(q.id), 'Invalid or duplicate provider/model quote');
    ids.add(q.id);
    requireCondition(q.fetchedAt === source.fetchedAt && q.openWeights === true && q.billing === 'flat' && typeof q.cacheSupported === 'boolean', 'Invalid quote snapshot/billing');
    requireCondition(q.providerDocs === null || safeHttps(q.providerDocs) === q.providerDocs, 'Invalid provider documentation URL');
    requireCondition(q.contextLimit === null || (Number.isSafeInteger(q.contextLimit) && q.contextLimit > 0), 'Invalid context limit');
    for (const bucket of ['peak', 'off']) for (const key of ['hit', 'miss', 'output']) requireCondition(finitePrice(q.prices?.[bucket]?.[key]), 'Invalid token price');
    for (const key of ['hit', 'miss', 'output']) requireCondition(q.prices.peak[key] === q.prices.off[key], 'Single catalogue rate differs between periods');
    if (!q.cacheSupported) requireCondition(q.prices.peak.hit === q.prices.peak.miss, 'Missing cache price must use ordinary input price');
    requireCondition(isObject(q.upstreamCost) && q.upstreamCost.input === q.prices.peak.miss && q.upstreamCost.output === q.prices.peak.output, 'Upstream price provenance mismatch');
    requireCondition(q.cacheSupported ? q.upstreamCost.cache_read === q.prices.peak.hit : q.upstreamCost.cache_read === undefined, 'Upstream cache provenance mismatch');
    requireCondition(Object.keys(q.upstreamCost).every(key => ['input', 'output', 'cache_read', 'cache_write', 'reasoning'].includes(key)), 'Unsupported price dimension in snapshot');
    requireCondition(q.upstreamCost.cache_write === undefined || q.upstreamCost.cache_write === 0, 'Unsupported cache-write cost in snapshot');
    requireCondition(q.upstreamCost.reasoning === undefined || q.upstreamCost.reasoning === q.upstreamCost.output, 'Unsupported reasoning cost in snapshot');
  }
  const coverage = catalog.coverage;
  requireCondition(isObject(coverage) && coverage.importedQuotes === catalog.quotes.length && coverage.importedProviders === new Set(catalog.quotes.map(q => q.providerId)).size, 'Invalid coverage metadata');
  requireCondition(Number.isSafeInteger(coverage.scannedProviders) && coverage.scannedProviders >= coverage.importedProviders && Number.isSafeInteger(coverage.scannedModels) && coverage.scannedModels >= coverage.importedQuotes && isObject(coverage.excludedCounts), 'Invalid source coverage');
  requireCondition(Object.values(coverage.excludedCounts).every(v => Number.isSafeInteger(v) && v >= 0) && Object.values(coverage.excludedCounts).reduce((a, b) => a + b, 0) + coverage.importedQuotes === coverage.scannedModels, 'Invalid exclusion counts');
  return catalog;
}

export async function updateCatalog(previous, fetchData, now = new Date().toISOString()) {
  if (previous) validateCatalog(previous);
  requireCondition(timestamp(now), 'Invalid update timestamp');
  try {
    const { quotes, coverage } = parseModelsDev(await fetchData(PRICE_SOURCE), now);
    // An unexpectedly truncated catalogue must not replace a broad working snapshot.
    if (previous?.quotes.length >= 100) requireCondition(quotes.length >= previous.quotes.length / 2, 'Catalogue lost more than half its supported quotes');
    const source = { ...PRICE_SOURCE, fetchedAt: now, lastAttemptAt: now, status: 'ok' };
    return { catalog: validateCatalog({ schemaVersion: 2, generatedAt: now, sources: [source], quotes, coverage }), failures: [] };
  } catch (error) {
    if (!previous) throw new Error(`models.dev: no validated snapshot is available (${error.message})`);
    const source = { ...previous.sources[0], lastAttemptAt: now, status: 'retained' };
    return {
      catalog: validateCatalog({ ...previous, generatedAt: now, sources: [source] }),
      failures: [{ sourceId: PRICE_SOURCE.id, error: String(error.message) }],
    };
  }
}
