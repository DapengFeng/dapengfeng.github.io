const buckets = ['peak', 'off'], rateTypes = ['hit', 'miss', 'output'];
const priceNames = buckets.flatMap(b => rateTypes.map(t => `prices.${b}.${t}`));
const moneyNames = ['capex', 'monthlyOpex', 'terminalResidual', 'perMillionLocalInputCost', 'perMillionLocalOutputCost', ...priceNames];
const text = (value, limit = 1000) => typeof value === 'string' && value.length <= limit;
const date = value => typeof value === 'string' && Number.isFinite(Date.parse(value));
export function safeCatalogUrl(value) {
 try {const u = new URL(value); return u.href === 'https://models.dev/api.json' ? u.href : null;} catch {return null;}
}
/** Treat the shipped catalog as data: validate tariff units and provenance before use. */
export function validateCatalog(data) {
 if (!data || data.schemaVersion !== 2 || !date(data.generatedAt) || !Array.isArray(data.sources) || !Array.isArray(data.quotes) || data.quotes.length > 10000) throw new Error('Invalid pricing catalog');
 const sources = new Map();
 for (const s of data.sources) {
  if (!text(s.id, 120) || sources.has(s.id) || !safeCatalogUrl(s.url) || !date(s.fetchedAt) || !date(s.lastAttemptAt) || !['ok', 'retained'].includes(s.status)) throw new Error('Invalid pricing source');
  for (const key of ['label', 'tierDefault', 'conditionsDefault', 'conditionsNoCache']) if (!text(s[key]?.en) || !text(s[key]?.zh)) throw new Error('Invalid source description');
  sources.set(s.id, s);
 }
 const ids = new Set();
 for (const q of data.quotes) {
  const source = sources.get(q.sourceId);
  if (!text(q.id, 2000) || !q.id || ids.has(q.id) || !source || q.sourceUrl !== source.url || q.fetchedAt !== source.fetchedAt || !safeCatalogUrl(q.sourceUrl) || !date(q.fetchedAt) || q.currency !== 'USD' || q.billing !== 'flat' || q.openWeights !== true || typeof q.cacheSupported !== 'boolean') throw new Error('Invalid quote');
  for (const k of ['providerId', 'provider', 'modelId', 'modelName', 'revision']) if (!text(q[k], 1024) || !q[k].trim()) throw new Error('Invalid model identity');
  for (const k of ['tier', 'conditions']) if (q[k] && (!text(q[k].en) || !text(q[k].zh))) throw new Error('Invalid quote conditions');
  for (const b of buckets) for (const type of rateTypes) if (!Number.isFinite(q.prices?.[b]?.[type]) || q.prices[b][type] < 0) throw new Error('Invalid token price');
  if (!q.cacheSupported && buckets.some(b => q.prices[b].hit !== q.prices[b].miss)) throw new Error('Unsupported cache tariff differs from uncached price');
  if (q.billing === 'flat' && rateTypes.some(t => q.prices.peak[t] !== q.prices.off[t])) throw new Error('Flat tariff differs by time');
  ids.add(q.id);
 }
 return {...data, quotes: data.quotes.map(q => {
  const source = sources.get(q.sourceId);
  return {...q, sourceLabel: source.label, tier: q.tier ?? source.tierDefault, conditions: q.conditions ?? (q.cacheSupported ? source.conditionsDefault : source.conditionsNoCache)};
 })};
}

function bilingual(parent, en, zh) {
 for (const [lang, value] of [['en', en], ['zh', zh]]) {
  const span = document.createElement('span');
  span.lang = lang === 'zh' ? 'zh-CN' : 'en'; span.dataset.lang = lang; span.textContent = value; parent.append(span);
 }
}
const identity = q => q ? `${q.provider}|${q.modelId}|${q.revision}` : '';
const dateLabel = s => String(s).slice(0, 10);

export function createPricingUI(form, {onChange, resetParameters}) {
 const select = form.elements.namedItem('presetId'), currency = form.elements.namedItem('currency');
 const provider = form.elements.namedItem('modelProvider'), search = form.elements.namedItem('modelSearch');
 const name = form.elements.namedItem('modelName'), version = form.elements.namedItem('modelVersion');
 const sourceEl = form.querySelector('[data-source-message]'), statusEl = form.querySelector('[data-catalog-status]');
 let catalog = null, quote = null, catalogFailed = false, notice = null;
 const field = key => form.elements.namedItem(key);
 const currentLanguage = () => document.documentElement.dataset.language || 'both';
 const localize = (en, zh) => currentLanguage() === 'both' ? `${en} / ${zh}` : currentLanguage() === 'zh' ? zh : en;
 function matches(q) {
  const query = search?.value.trim().toLocaleLowerCase() || '';
  return q.currency === currency.value && (!provider?.value || q.providerId === provider.value) &&
   (!query || `${q.provider} ${q.providerId} ${q.modelName} ${q.modelId} ${q.canonicalModelId ?? ''}`.toLocaleLowerCase().includes(query));
 }
 function say(el, en, zh) {el.replaceChildren(); bilingual(el, en, zh);}
 function updateStatus() {
  if (!statusEl) return;
  if (catalogFailed) {say(statusEl, 'The price catalog is unavailable. Enter prices manually.', '报价目录暂不可用，可手动填写单价。'); return;}
  if (!catalog) {say(statusEl, 'Loading the price catalog…', '正在读取报价目录…'); return;}
  const count = catalog.quotes.filter(matches).length;
  if (currency.value !== 'USD') {say(statusEl, 'This catalog lists USD prices. Enter CNY prices manually.', '目录提供美元报价，人民币单价请手动填写。'); return;}
  const retainedSelection = quote && !matches(quote);
  say(statusEl, `${count} matching model/provider entries · USD${retainedSelection ? ' · current selection retained' : ''}`, `${count} 项匹配的模型与服务商报价 · USD${retainedSelection ? '；保留当前选择' : ''}`);
 }
 function updateProviders() {
  if (!provider) return;
  const previous = provider.value, providers = new Map();
  for (const q of catalog?.quotes || []) if (q.currency === currency.value) providers.set(q.providerId, q.provider);
  provider.replaceChildren(new Option(localize('All providers', '全部服务商'), ''));
  for (const [id, name] of [...providers].sort((a, b) => a[1].localeCompare(b[1]))) provider.append(new Option(name, id));
  provider.value = providers.has(previous) ? previous : '';
 }
 function updateOptions() {
  const previous = select.value;
  select.replaceChildren();
  const filtered = catalog?.quotes.filter(matches) || [];
  const selected = catalog?.quotes.find(q => q.id === previous && q.currency === currency.value);
  const visible = selected && !filtered.includes(selected) ? [selected, ...filtered] : filtered;
  const options = [
   {value: '', en: 'Choose a model', zh: '选择模型'},
   {value: 'custom', en: 'Enter prices manually', zh: '手动填写价格'},
   ...visible.map(q => ({value: q.id, en: `${q.provider} · ${q.modelId}`, zh: `${q.provider} · ${q.modelId}`})),
  ];
  for (const item of options) {
   const option = document.createElement('option'); option.value = item.value;
   option.dataset.en = item.en; option.dataset.zh = item.zh;
   option.textContent = item.en === item.zh ? item.en : localize(item.en, item.zh);
   if (!item.value) option.disabled = true;
   select.append(option);
  }
  select.value = options.some(o => o.value === previous) ? previous : '';
  updateStatus();
 }
 function overrides() {
  if (!quote) return {prices: false, identity: false};
  return {
   prices: buckets.some(b => rateTypes.some(t => field(`prices.${b}.${t}`).valueAsNumber !== quote.prices[b][t])),
   identity: name.value.trim() !== quote.modelName || version.value.trim() !== quote.revision,
  };
 }
 function provenance() {
  if (!quote) return null;
  const source = catalog.sources.find(s => s.id === quote.sourceId);
  return {
   ...quote, sourceStatus: source.status, lastAttemptAt: source.lastAttemptAt,
   retained: source.status === 'retained', stale: Date.now() - Date.parse(quote.fetchedAt) > 3 * 86400000,
   overrides: overrides(),
  };
 }
 function showSource() {
  const wasOpen = Boolean(sourceEl.querySelector('details')?.open);
  sourceEl.replaceChildren();
  if (quote) {
   const p = provenance(), changed = p.overrides.prices || p.overrides.identity;
   const first = document.createElement('p');
   bilingual(first, `models.dev · snapshot ${dateLabel(quote.fetchedAt)}${changed ? ' · user modified' : ''}`, `models.dev · 快照 ${dateLabel(quote.fetchedAt)}${changed ? ' · 用户已修改' : ''}`);
   sourceEl.append(first);
   const details = document.createElement('details'); details.open = wasOpen;
   const summary = document.createElement('summary'); bilingual(summary, 'Price source and limits', '报价来源与适用范围'); details.append(summary);
   const a = document.createElement('a'); a.href = safeCatalogUrl(quote.sourceUrl); a.target = '_blank'; a.rel = 'noopener noreferrer'; bilingual(a, 'CC Switch upstream catalog', 'CC Switch 同源目录'); details.append(a);
   const conditions = document.createElement('p'); bilingual(conditions, quote.conditions.en, quote.conditions.zh); details.append(conditions);
   const key = document.createElement('p'); key.textContent = `${quote.provider} · ${quote.modelId}`; details.append(key);
   if (Number.isSafeInteger(quote.contextLimit) && quote.contextLimit > 0) {
    const context = document.createElement('p'), count = new Intl.NumberFormat('en-US').format(quote.contextLimit);
    bilingual(context, `Catalog context limit: ${count} tokens.`, `目录标注的上下文上限：${count} token。`); details.append(context);
   }
   sourceEl.append(details);
   if (p.retained || p.stale) {
    const warning = document.createElement('p'); warning.className = 'dr-price-warning';
    bilingual(warning, p.retained ? `The latest update failed; this earlier quote is retained. Last attempt: ${dateLabel(p.lastAttemptAt)}.` : 'This price snapshot is more than three days old.', p.retained ? `本轮更新未成功，沿用此前报价；最近尝试 ${dateLabel(p.lastAttemptAt)}。` : '此价格快照已超过三天。'); sourceEl.append(warning);
   }
   if (!quote.cacheSupported) {const p = document.createElement('p'); bilingual(p, 'No cache-read tariff is listed; both cache-hit rates are fixed at 0%.', '未列缓存读取价，两项命中率固定为 0%。'); details.append(p);}
  } else if (select.value === 'custom') {
   const p = document.createElement('p'); bilingual(p, 'Manually entered prices', '手动填写单价'); sourceEl.append(p);
  }
  if (notice) {const p = document.createElement('p'); p.className = 'dr-price-warning'; bilingual(p, notice.en, notice.zh); sourceEl.append(p);}
 }
 function setCachePolicy() {
  const unsupported = quote && !quote.cacheSupported;
  for (const key of ['apiCacheHit', 'fallbackCacheHit']) {field(key).readOnly = Boolean(unsupported); if (unsupported) field(key).value = '0';}
 }
 function clearCapacity() {field('inputTokensPerSecond').value = '0'; field('outputTokensPerSecond').value = '0';}
 function clearPrices() {for (const key of priceNames) field(key).value = '0';}
 function applyQuote(id) {
  const next = catalog?.quotes.find(q => q.id === id && q.currency === currency.value);
  if (!next) {
   quote = null; clearPrices(); clearCapacity(); name.value = ''; version.value = ''; setCachePolicy();
   const identityDetails = name.closest('details');
   if (identityDetails) identityDetails.open = true;
   notice = {en: 'Enter the model, its prices and measured throughput.', zh: '请填写模型、单价与实测吞吐。'}; updateStatus(); showSource(); onChange(); return;
  }
  const changedModel = !quote || identity(quote) !== identity(next) || overrides().identity;
  quote = next;
  for (const b of buckets) for (const t of rateTypes) field(`prices.${b}.${t}`).value = String(next.prices[b][t]);
  name.value = next.modelName; version.value = next.revision;
  if (changedModel) clearCapacity();
  setCachePolicy();
  notice = changedModel ? {en: 'Model changed: both throughput fields are cleared. Enter measurements for this model.', zh: '模型已改变，两项吞吐已清零，请填入该模型的实测值。'} : null;
  updateStatus(); showSource(); onChange();
 }
 function reset() {
  resetParameters(); select.value = ''; quote = null; name.value = ''; version.value = ''; notice = null;
  if (provider) provider.value = '';
  if (search) search.value = '';
  updateProviders(); updateOptions(); setCachePolicy(); showSource(); onChange();
 }
 select.addEventListener('change', () => applyQuote(select.value));
 provider?.addEventListener('change', updateOptions);
 search?.addEventListener('input', updateOptions);
 currency.addEventListener('change', () => {
  quote = null; for (const key of moneyNames) field(key).value = '0';
  select.value = name.value.trim() ? 'custom' : '';
  if (provider) provider.value = '';
  if (search) search.value = '';
  updateProviders(); updateOptions(); setCachePolicy();
  notice = {en: 'Currency changed: all prices and monetary costs were cleared. Enter amounts in the selected currency.', zh: '币种已改变，单价和金额已清零，请用所选币种重新填写。'};
  showSource(); onChange();
 });
 form.addEventListener('input', event => {
  const key = event.target.name;
  if (!['modelName', 'modelVersion', ...priceNames].includes(key)) return;
  if (!select.value) select.value = 'custom';
  if (key === 'modelName' || key === 'modelVersion') {
   clearCapacity(); notice = {en: 'Model identity changed: enter its measured throughput again.', zh: '模型标识已修改，请重新填写该模型的实测吞吐。'};
  }
  showSource(); onChange();
 });
 window.addEventListener('languagechange', () => {updateProviders(); updateOptions(); showSource();});
 async function load() {
  updateStatus(); const controller = new AbortController(), timer = setTimeout(() => controller.abort(), 8000);
  try {
   const response = await fetch(new URL('./model-prices.json', import.meta.url), {signal: controller.signal, cache: 'no-cache', credentials: 'same-origin'});
   if (!response.ok) throw new Error('Catalog unavailable');
   const raw = await response.text(); if (raw.length > 15000000) throw new Error('Catalog too large');
   catalog = validateCatalog(JSON.parse(raw)); catalogFailed = false; updateProviders(); updateOptions(); showSource();
  } catch {catalogFailed = true; updateOptions(); showSource();} finally {clearTimeout(timer);}
 }
 updateProviders(); updateOptions();
 return {
  reset, load,
  isManagedInput: key => ['modelProvider', 'modelSearch', 'modelName', 'modelVersion', ...priceNames].includes(key),
  metadata: () => ({presetId: quote ? quote.id : 'custom', name: name.value.trim(), version: version.value.trim(), currency: currency.value, pricing: provenance(), modified: quote ? Object.values(overrides()).some(Boolean) : false}),
  hasSelection: () => Boolean(select.value),
  enforceCachePolicy: setCachePolicy,
 };
}
