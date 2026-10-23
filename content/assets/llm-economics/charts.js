/** Live, model-neutral views of the same token scenario used by the calculator.
 * No DOM is bound at import time. Caller supplies a validated result and identity.
 * SVG uses actual container width so labels retain their CSS-pixel size on phones.
 */
import {evaluateTokenScenario} from './token-model.js';

const BLUE = '#94d4ec', MINT = '#bcebc5', AMBER = '#efbd88', SLATE = '#b3becb';
const instances = new WeakMap();
const svgNS = 'http://www.w3.org/2000/svg';
const numeric = n => {
  if (!Number.isFinite(n)) return '—';
  const v = Math.abs(n);
  if ((v > 0 && v < .0001) || v >= 1e12) return n.toExponential(3);
  return new Intl.NumberFormat('en-US', {maximumSignificantDigits: 6}).format(Object.is(n, -0) ? 0 : n);
};
const compact = n => {
  const v = Math.abs(n);
  if (v >= 1e12 || (v > 0 && v < .001)) return n.toExponential(1);
  const [div, unit] = v >= 1e9 ? [1e9, 'B'] : v >= 1e6 ? [1e6, 'M'] : v >= 1e3 ? [1e3, 'k'] : [1, ''];
  return `${new Intl.NumberFormat('en-US', {maximumFractionDigits: 2}).format(Object.is(n, -0) ? 0 : n / div)}${unit}`;
};
const percent = n => `${numeric(100 * n)}%`;
const labels = {
  prices: ['API prices used in this comparison', '本次比较采用的 API 单价'],
  cost: ['Month 1 · how much does local save?', '首月：本地能省多少？'],
  cash: ['When does local deployment pay back?', '何时收回初投？'],
  sensitivity: ['If demand falls, is local still worth it?', '需求下降后，还值得部署吗？'],
  waiting: ['Enter monthly input or output token demand to calculate this view.', '填写月输入或输出 token 用量后，这里显示对应计算。'],
};

/** Pure data preparation is exported to check charts against the financial model. */
export function prepareChartData(result, metadata = {}) {
  const p = result.parameters;
  const prices = ['hit', 'miss', 'output'].map(key => ({key, peak: p.prices.peak[key], off: p.prices.off[key]}));
  const hasDemand = p.monthlyInputTokens > 0 || p.monthlyOutputTokens > 0;
  const demandMultipliers = [...new Set([...Array.from({length: 33}, (_, i) => .25 + 2.75 * i / 32), 1])].sort((a, b) => a - b);
  const sensitivity = hasDemand ? demandMultipliers.map(multiplier => ({
    multiplier,
    monthlyInputTokens: p.monthlyInputTokens * multiplier,
    monthlyOutputTokens: p.monthlyOutputTokens * multiplier,
    NPV: multiplier === 1 ? result.NPV : evaluateTokenScenario({
      ...p, monthlyInputTokens: p.monthlyInputTokens * multiplier, monthlyOutputTokens: p.monthlyOutputTokens * multiplier,
    }).NPV,
  })) : [];
  const m = result.firstMonth;
  const final = result.months.at(-1);
  return {
    currency: String(metadata.currency || 'CURRENCY'),
    model: metadata.modelLabel || metadata.name || metadata.label || {en: 'Current model', zh: '当前模型'},
    provider: String(metadata.pricing?.provider || metadata.provider || ''),
    providerId: String(metadata.pricing?.providerId || metadata.providerId || ''),
    modified: Boolean(metadata.modified),
    cachePriceListed: metadata.pricing?.cacheSupported !== false,
    version: String(metadata.version || ''),
    hasDemand, prices, flatTariff: prices.every(row => row.peak === row.off),
    cost: {api: m.fullApiCost, fixed: m.fixedOpex, variable: m.variableLocalCost, fallback: m.overflowApiCost, localTotal: m.hybridOperatingCost, saving: m.operatingSavings},
    cash: result.months.map(row => ({month: row.month, ordinary: row.cumulativeOperatingNet, discounted: row.cumulativeDiscountedOperatingNet})),
    terminal: {month: p.horizon, amount: p.terminalResidual, presentValue: result.presentValueResidual, ordinary: final.cumulativeNet, discounted: final.cumulativeDiscountedNet},
    sensitivity, parameters: p,
  };
}

function element(doc, tag, cls, text) {
  const el = doc.createElement(tag);
  if (cls) el.className = cls;
  if (text !== undefined) el.textContent = String(text);
  return el;
}
function pair(doc, en, zh, cls = '') {
  const el = element(doc, 'span', cls);
  for (const [lang, text] of [['en', en], ['zh', zh]]) {
    const span = element(doc, 'span', '', text);
    span.dataset.lang = lang; span.lang = lang === 'zh' ? 'zh-CN' : 'en'; el.append(span);
  }
  return el;
}
function pairedValue(value) {
  return typeof value === 'object' && value !== null ? [String(value.en || value.zh || ''), String(value.zh || value.en || '')] : [String(value), String(value)];
}
function svgElement(doc, tag, attrs = {}, text) {
  const el = doc.createElementNS(svgNS, tag);
  for (const [key, value] of Object.entries(attrs)) el.setAttribute(key, String(value));
  if (text !== undefined) el.textContent = String(text);
  return el;
}
function note(doc, en, zh, cls = 'lc-note') {
  const p = element(doc, 'p', cls); p.append(pair(doc, en, zh)); return p;
}
function panel(doc, id, number) {
  const el = element(doc, 'section', 'lc-panel'); el.dataset.chart = id;
  const header = element(doc, 'div', 'lc-panel-heading');
  header.append(element(doc, 'span', 'lc-number', String(number).padStart(2, '0')));
  const title = element(doc, 'h3', 'lc-title'); title.append(pair(doc, ...labels[id])); header.append(title); el.append(header);
  return el;
}
function money(value, data) { return `${data.currency} ${numeric(value)}`; }
function valueLine(doc, en, zh, value, data, extraClass = '') {
  const row = element(doc, 'div', `lc-value-line ${extraClass}`);
  row.append(pair(doc, en, zh), element(doc, 'strong', '', money(value, data)));
  row.dataset.value = String(value); return row;
}
function legendItem(doc, en, zh, color, dashed = false, point = false) {
  const item = element(doc, 'span', 'lc-legend-item');
  const mark = element(doc, 'i', `lc-key${dashed ? ' lc-key-dashed' : ''}${point ? ' lc-key-point' : ''}`);
  mark.style.setProperty('--lc-key-color', color); mark.setAttribute('aria-hidden', 'true');
  item.append(mark, pair(doc, en, zh)); return item;
}

function renderPrices(doc, data) {
  const out = panel(doc, 'prices', 1);
  out.append(note(doc, `${data.currency} per million tokens · input and output are billed separately`, `每百万 token，单位 ${data.currency} · 输入、输出分别计费`));
  const grid = element(doc, 'div', 'lc-price-grid');
  const names = {hit: ['Cached input', '命中缓存的输入'], miss: ['Uncached input', '未命中缓存的输入'], output: ['Output', '输出']};
  const max = Math.max(...data.prices.filter(row => row.key !== 'hit' || data.cachePriceListed).flatMap(row => [row.peak, row.off]));
  for (const row of data.prices) {
    const item = element(doc, 'div', 'lc-price-category'); item.dataset.category = row.key;
    const name = element(doc, 'strong', 'lc-category-title'); name.append(pair(doc, ...names[row.key])); item.append(name);
    if (row.key === 'hit' && !data.cachePriceListed) {
      const unavailable = note(doc, 'No cache-read price listed. Cache hits are 0%; input uses the ordinary input price.', '未列缓存读取价。命中率为 0%，输入按普通输入价计费。');
      unavailable.dataset.price = 'peak.hit'; unavailable.dataset.value = row.peak; unavailable.dataset.unlisted = 'true'; item.append(unavailable);
      if (!data.flatTariff) {
        const fallback = element(doc, 'span'); fallback.hidden = true; fallback.dataset.price = 'off.hit'; fallback.dataset.value = row.off; fallback.dataset.unlisted = 'true'; item.append(fallback);
      }
      grid.append(item); continue;
    }
    for (const bucket of data.flatTariff ? ['peak'] : ['peak', 'off']) {
      const label = data.flatTariff ? ['Flat tariff', '统一单价'] : bucket === 'peak' ? ['Peak', '高峰'] : ['Off-peak', '空闲'];
      const entry = element(doc, 'div', 'lc-price-entry'); entry.dataset.price = `${bucket}.${row.key}`; entry.dataset.value = row[bucket];
      const line = element(doc, 'div', 'lc-price-label'); line.append(pair(doc, ...label), element(doc, 'b', '', numeric(row[bucket])));
      const bar = element(doc, 'div', 'lc-track'); bar.setAttribute('aria-hidden', 'true');
      const fill = element(doc, 'div', 'lc-price-bar'); fill.style.width = `${max > 0 ? row[bucket] / max * 100 : 0}%`; fill.style.background = bucket === 'peak' ? BLUE : MINT; bar.append(fill);
      entry.append(line, bar); item.append(entry);
    }
    grid.append(item);
  }
  out.append(grid);
  out.append(note(doc, max > 0 ? 'All price bars share one scale; 0 is the left edge.' : 'All displayed prices are zero.', max > 0 ? '各价格条使用同一标尺，左端均为 0。' : '展示的各项单价均为零。'));
  return out;
}
function renderCost(doc, data) {
  const out = panel(doc, 'cost', 2);
  if (!data.hasDemand) {out.append(note(doc, ...labels.waiting, 'lc-empty')); return out;}
  out.append(note(doc, `A is the full API bill. L includes local running costs and fallback API. Amounts are in ${data.currency}; the initial investment is separate.`, `A 是全用 API 的账单，L 包含本地运行与回退 API 的全部费用。金额单位为 ${data.currency}，购机初投另计。`));
  const max = Math.max(data.cost.api, data.cost.localTotal);
  for (const definition of [
    {id: 'api', name: ['A · All API', 'A · 全用 API'], total: data.cost.api, parts: [{id: 'api', value: data.cost.api, color: BLUE}]},
    {id: 'local', name: ['L · Local running + API fallback', 'L · 本地运行＋回退 API'], total: data.cost.localTotal, parts: [{id: 'fixed', value: data.cost.fixed, color: AMBER}, {id: 'variable', value: data.cost.variable, color: MINT}, {id: 'fallback', value: data.cost.fallback, color: BLUE}]},
  ]) {
    const row = element(doc, 'div', 'lc-cost-row'); row.dataset.series = definition.id;
    row.append(valueLine(doc, ...definition.name, definition.total, data));
    const bar = element(doc, 'div', 'lc-cost-bar'); bar.setAttribute('aria-hidden', 'true');
    for (const part of definition.parts) {
      const segment = element(doc, 'span', 'lc-cost-segment'); segment.style.width = `${max > 0 ? part.value / max * 100 : 0}%`; segment.style.background = part.color;
      segment.dataset.component = part.id; segment.dataset.value = part.value; bar.append(segment);
    }
    row.append(bar); out.append(row);
  }
  const scale = element(doc, 'div', 'lc-scale'); scale.append(element(doc, 'span', '', '0'), element(doc, 'span', '', money(max, data))); out.append(scale);
  const components = element(doc, 'dl', 'lc-cost-components');
  for (const [key, en, zh, color] of [['fixed', 'Fixed operations', '固定运营', AMBER], ['variable', 'Extra local token costs', '额外本地 token 成本', MINT], ['fallback', 'Fallback API bill', '回退 API 账单', BLUE]]) {
    const item = element(doc, 'div'); const title = element(doc, 'dt'); title.append(legendItem(doc, en, zh, color));
    const amount = element(doc, 'dd', '', money(data.cost[key], data)); amount.dataset.componentTotal = key; amount.dataset.value = data.cost[key]; item.append(title, amount); components.append(item);
  }
  out.append(components);
  const saving = valueLine(doc, 'S = A − L · first-month saving', 'S = A − L · 首月节省', data.cost.saving, data, 'lc-saving');
  saving.dataset.firstMonthSaving = '';
  saving.dataset.negative = String(data.cost.saving < 0); out.append(saving);
  out.append(note(doc, 'A positive S saves money; a negative S means local costs more that month.', 'S 为正代表当月省钱，为负代表当月比全用 API 多花钱。')); return out;
}

function chartGeometry(width, values, xRange) {
  const w = Math.max(160, width), right = 18, top = 22, bottom = 38, height = w < 480 ? 268 : 326;
  let low = Math.min(0, ...values), high = Math.max(0, ...values);
  if (low === high) {low = -1; high = 1;}
  const rough = (high - low) / 4, power = 10 ** Math.floor(Math.log10(rough));
  const normalized = rough / power;
  const step = (normalized <= 1 ? 1 : normalized <= 2 ? 2 : normalized <= 2.5 ? 2.5 : normalized <= 5 ? 5 : 10) * power;
  const tickLow = Math.floor(low / step) * step, tickHigh = Math.ceil(high / step) * step;
  const ticks = Array.from({length: Math.round((tickHigh - tickLow) / step) + 1}, (_, i) => tickLow + step * i);
  const left = Math.min(w * .42, Math.max(w < 480 ? 56 : 72, ...ticks.map(value => compact(value).length * 7.6 + 12)));
  low = tickLow - step * .13; high = tickHigh + step * .13;
  const x = value => left + (value - xRange[0]) / (xRange[1] - xRange[0]) * (w - left - right);
  const y = value => top + (high - value) / (high - low) * (height - top - bottom);
  return {w, height, left, right, top, bottom, low, high, x, y, ticks};
}

function svgBase(doc, geometry, title, xTicks) {
  const {w, height, left, right, bottom, x, y, ticks} = geometry;
  const svg = svgElement(doc, 'svg', {viewBox: `0 0 ${w} ${height}`, width: '100%', height, role: 'img', 'aria-label': title, preserveAspectRatio: 'xMidYMid meet'});
  svg.classList.add('lc-svg');
  for (const value of ticks) {
    svg.append(svgElement(doc, 'line', {x1: left, x2: w - right, y1: y(value), y2: y(value), class: 'lc-grid-line'}));
    svg.append(svgElement(doc, 'text', {x: left - 9, y: y(value), 'text-anchor': 'end', 'dominant-baseline': 'middle', class: 'lc-tick'}, compact(value)));
  }
  svg.append(svgElement(doc, 'line', {x1: left, x2: w - right, y1: y(0), y2: y(0), class: 'lc-zero-line', 'data-zero': ''}));
  for (const value of xTicks) svg.append(svgElement(doc, 'text', {x: x(value), y: height - bottom + 25, 'text-anchor': 'middle', class: 'lc-tick'}, numeric(value)));
  return svg;
}
function seriesPath(doc, points, geometry, {id, color, dashed = false}) {
  const path = svgElement(doc, 'path', {
    d: points.map(([x, y], i) => `${i === 0 ? 'M' : 'L'}${geometry.x(x).toFixed(3)},${geometry.y(y).toFixed(3)}`).join(' '),
    fill: 'none', stroke: color, 'stroke-width': 2.5, 'stroke-linejoin': 'round', 'stroke-linecap': 'round', 'data-series': id,
  });
  if (dashed) path.setAttribute('stroke-dasharray', '7 5');
  path.dataset.values = JSON.stringify(points); return path;
}
function marker(doc, geometry, month, value, color, id) {
  const circle = svgElement(doc, 'circle', {cx: geometry.x(month), cy: geometry.y(value), r: 4.5, fill: '#111e2b', stroke: color, 'stroke-width': 2, 'data-series': id, 'data-value': value, 'data-month': month});
  circle.append(svgElement(doc, 'title', {}, `${month}: ${numeric(value)}`)); return circle;
}
function renderCash(doc, data, drawingJobs) {
  const out = panel(doc, 'cash', 3);
  if (!data.hasDemand) {out.append(note(doc, ...labels.waiting, 'lc-empty')); return out;}
  out.append(note(doc, 'Start below zero by the initial investment, then add each month’s saving S. The first zero crossing recovers the investment; a later drop below zero reverses that recovery.', '起点是负的初始投资，之后逐月加上节省 S。曲线首次达到零线时收回初投；若后来又落回零线下，则已收回的资金会再次被消耗。'));
  const legend = element(doc, 'div', 'lc-legend');
  legend.append(legendItem(doc, 'Ordinary cash', '未折现现金', BLUE), legendItem(doc, 'Discounted cash', '折现现金', MINT, true), legendItem(doc, 'With terminal value', '含期末出售所得', SLATE, false, true));
  out.append(legend, note(doc, `Cumulative balance · ${data.currency}`, `累计余额 · ${data.currency}`, 'lc-axis-unit'));
  const plot = element(doc, 'div', 'lc-plot'); out.append(plot);
  const xTitle = element(doc, 'div', 'lc-x-title'); xTitle.append(pair(doc, 'Month', '月份')); out.append(xTitle);
  drawingJobs.push(() => {
    const final = data.cash.at(-1), t = data.terminal;
    const g = chartGeometry(plot.clientWidth, [...data.cash.flatMap(r => [r.ordinary, r.discounted]), t.ordinary, t.discounted], [0, t.month]);
    const ticks = [...new Set([0, Math.round(t.month / 3), Math.round(2 * t.month / 3), t.month])];
    const svg = svgBase(doc, g, `Cumulative cash / 累计现金 · ${data.currency}`, ticks);
    svg.append(seriesPath(doc, data.cash.map(r => [r.month, r.ordinary]), g, {id: 'ordinary', color: BLUE}));
    svg.append(seriesPath(doc, data.cash.map(r => [r.month, r.discounted]), g, {id: 'discounted', color: MINT, dashed: true}));
    for (const [key, color] of [['ordinary', BLUE], ['discounted', MINT]]) {
      if (t.amount !== 0) svg.append(svgElement(doc, 'line', {x1: g.x(t.month), x2: g.x(t.month), y1: g.y(final[key]), y2: g.y(t[key]), stroke: color, 'stroke-width': 1.5, 'stroke-dasharray': '2 3', 'data-terminal-connector': key}));
      svg.append(marker(doc, g, t.month, t[key], color, `terminal-${key}`));
    }
    plot.replaceChildren(svg);
  });
  const terminal = element(doc, 'div', 'lc-terminal-note');
  terminal.append(valueLine(doc, `Terminal cash at month ${data.terminal.month}`, `第 ${data.terminal.month} 月期末现金`, data.terminal.amount, data));
  terminal.append(valueLine(doc, 'Present value of terminal cash', '期末现金的现值', data.terminal.presentValue, data)); out.append(terminal);
  out.append(note(doc, 'Lines exclude terminal value; open dots add it at the final month. Operating payback uses the lines, not the dots.', '曲线不计期末出售所得，空心点在最后一个月加上该笔现金；运营回本按曲线判定，不使用空心点。'));
  return out;
}
function fixedParameters(doc, data) {
  const p = data.parameters, details = element(doc, 'details', 'lc-assumptions');
  const summary = element(doc, 'summary'); summary.append(pair(doc, 'Inputs held fixed in this sweep', '本轮敏感度计算中保持不变的参数')); details.append(summary);
  const list = element(doc, 'dl', 'lc-parameter-list');
  const fields = [
    ['inputTokensPerSecond', 'Paired input rate', '成对输入速率', `${numeric(p.inputTokensPerSecond * p.capacityMultiplier)} token/s`],
    ['outputTokensPerSecond', 'Paired output rate', '成对输出速率', `${numeric(p.outputTokensPerSecond * p.capacityMultiplier)} token/s`],
    ['availability', 'Availability', '可用率', percent(p.availability)], ['headroom', 'Scheduling fraction', '允许调度比例', percent(p.headroom)],
    ['apiCacheHit', 'Baseline cache hits', '基线缓存命中率', percent(p.apiCacheHit)], ['fallbackCacheHit', 'Fallback cache hits', '回退缓存命中率', percent(p.fallbackCacheHit ?? p.apiCacheHit)],
    ['inputPeakShare', 'Peak input share', '高峰输入占比', percent(p.inputPeakShare)], ['outputPeakShare', 'Peak output share', '高峰输出占比', percent(p.outputPeakShare)],
    ['peakHours', 'Peak hours/month', '月高峰小时', numeric(p.peakHours)], ['offHours', 'Off-peak hours/month', '月空闲小时', numeric(p.offHours)],
    ['capex', 'Initial investment', '初始投资', money(p.capex, data)], ['monthlyOpex', 'Fixed operations/month', '月固定运营', money(p.monthlyOpex, data)],
    ['perMillionLocalInputCost', 'Extra cost/million local input', '每百万本地输入额外成本', money(p.perMillionLocalInputCost, data)],
    ['perMillionLocalOutputCost', 'Extra cost/million local output', '每百万本地输出额外成本', money(p.perMillionLocalOutputCost, data)],
    ['horizon', 'Horizon (months)', '评估期限（月）', numeric(p.horizon)], ['annualDiscount', 'Annual discount', '年折现率', percent(p.annualDiscount)],
    ['terminalResidual', 'Net terminal cash', '期末净现金', money(p.terminalResidual, data)],
    ['inputAnnualGrowth', 'Input growth/year', '输入年增长', percent(p.inputAnnualGrowth)], ['outputAnnualGrowth', 'Output growth/year', '输出年增长', percent(p.outputAnnualGrowth)],
    ['apiAnnualChange', 'API price change/year', 'API 年变价', percent(p.apiAnnualChange)], ['opexAnnualGrowth', 'Operating-cost growth/year', '运营成本年增长', percent(p.opexAnnualGrowth)],
  ];
  for (const [key, en, zh, value] of fields) {
    const item = element(doc, 'div'); const name = element(doc, 'dt'); name.append(pair(doc, en, zh)); const amount = element(doc, 'dd', '', value);
    amount.dataset.parameter = key; amount.dataset.value = p[key] ?? p.apiCacheHit; item.append(name, amount); list.append(item);
  }
  details.append(list, note(doc, 'The current token prices above also stay fixed. This varies demand only; it does not assume faster hardware or a different task mix.', '上方当前 token 单价也保持不变。这里只改变需求，不假定硬件变快或任务组合变化。'));
  return details;
}
function renderSensitivity(doc, data, drawingJobs) {
  const out = panel(doc, 'sensitivity', 4);
  if (!data.hasDemand) {out.append(note(doc, ...labels.waiting, 'lc-empty')); return out;}
  out.append(note(doc, 'Left of 1× means lower demand. Below the zero line, local deployment has negative NPV at the entered discount rate. Input and output scale together; all other inputs stay fixed.', '1 倍左侧代表需求下降；落在零线下，表示按所填折现率计算，本地部署的净现值为负。输入与输出同比缩放，其他参数不变。'));
  const base = element(doc, 'div', 'lc-demand-baseline');
  base.append(pair(doc, `1× = ${compact(data.parameters.monthlyInputTokens)} input + ${compact(data.parameters.monthlyOutputTokens)} output tokens/month`, `1×＝每月 ${compact(data.parameters.monthlyInputTokens)} 输入＋${compact(data.parameters.monthlyOutputTokens)} 输出 token`));
  base.dataset.inputTokens = data.parameters.monthlyInputTokens; base.dataset.outputTokens = data.parameters.monthlyOutputTokens; out.append(base);
  out.append(note(doc, `NPV · ${data.currency}`, `净现值 · ${data.currency}`, 'lc-axis-unit'));
  const plot = element(doc, 'div', 'lc-plot'); out.append(plot);
  const xTitle = element(doc, 'div', 'lc-x-title'); xTitle.append(pair(doc, 'Token-demand multiplier', 'Token 需求倍率')); out.append(xTitle);
  drawingJobs.push(() => {
    const g = chartGeometry(plot.clientWidth, data.sensitivity.map(row => row.NPV), [.25, 3]);
    const svg = svgBase(doc, g, `NPV versus token demand / NPV 随 token 需求变化 · ${data.currency}`, [.25, 1, 2, 3]);
    svg.append(svgElement(doc, 'line', {x1: g.x(1), x2: g.x(1), y1: g.top, y2: g.height - g.bottom, class: 'lc-reference-line'}));
    svg.append(seriesPath(doc, data.sensitivity.map(r => [r.multiplier, r.NPV]), g, {id: 'demand-npv', color: MINT}));
    const baseline = data.sensitivity.find(row => row.multiplier === 1);
    const dot = marker(doc, g, 1, baseline.NPV, MINT, 'current-demand'); dot.removeAttribute('data-month'); dot.dataset.multiplier = '1'; svg.append(dot);
    plot.replaceChildren(svg);
  });
  const baseline = data.sensitivity.find(row => row.multiplier === 1);
  const current = valueLine(doc, 'Current inputs · 1× NPV', '当前输入 · 1 倍需求 NPV', baseline.NPV, data, 'lc-current-npv'); current.dataset.currentNpv = ''; out.append(current);
  out.append(fixedParameters(doc, data)); return out;
}

function dispose(container) {
  const state = instances.get(container);
  if (state) {state.observer?.disconnect(); if (state.frame !== null) state.win.cancelAnimationFrame(state.frame); instances.delete(container);}
}
export function clearCharts(container, message = {en: 'Complete valid inputs to update the charts.', zh: '填写有效参数后更新图表。'}) {
  if (!container) return;
  dispose(container); container.classList.add('dr-live-charts'); container.removeAttribute('data-model'); container.removeAttribute('data-provider'); container.removeAttribute('data-currency');
  container.replaceChildren(note(container.ownerDocument, ...pairedValue(message), 'lc-empty'));
}
export function renderCharts(container, result, metadata = {}) {
  if (!container) return null;
  dispose(container);
  const data = prepareChartData(result, metadata), doc = container.ownerDocument, win = doc.defaultView;
  container.classList.add('dr-live-charts'); container.dataset.model = pairedValue(data.model)[0]; container.dataset.provider = data.providerId || data.provider; container.dataset.currency = data.currency;
  const identity = element(doc, 'div', 'lc-identity');
  const identityLabel = element(doc, 'span', 'lc-identity-label'); identityLabel.append(pair(doc, 'API option compared', '本次比较的 API 方案')); identity.append(identityLabel);
  if (data.provider) identity.append(element(doc, 'strong', 'lc-provider-name', data.provider), element(doc, 'span', 'lc-identity-separator', '·'));
  const model = element(doc, 'strong', 'lc-model-name');
  if (typeof data.model === 'object') model.append(pair(doc, ...pairedValue(data.model))); else model.textContent = String(data.model);
  identity.append(model); if (data.version) identity.append(element(doc, 'span', 'lc-model-version', data.version));
  if (data.modified) identity.append(pair(doc, 'Edited inputs', '已手动修改', 'lc-identity-label'));
  const drawingJobs = [];
  container.replaceChildren(identity, renderPrices(doc, data), renderCost(doc, data), renderCash(doc, data, drawingJobs), renderSensitivity(doc, data, drawingJobs));
  const draw = () => {for (const job of drawingJobs) job();};
  draw();
  const state = {observer: null, frame: null, win, width: container.clientWidth}; instances.set(container, state);
  if (win?.ResizeObserver) {
    state.observer = new win.ResizeObserver(() => {
      const width = container.clientWidth;
      if (Math.abs(width - state.width) < 1) return;
      state.width = width;
      if (state.frame !== null) win.cancelAnimationFrame(state.frame);
      state.frame = win.requestAnimationFrame(() => {state.frame = null; draw();});
    });
    state.observer.observe(container);
  }
  return data;
}
