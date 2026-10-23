import {evaluateTokenScenario} from './token-model.js';
import {customDefaults} from './defaults.js';
import {createPricingUI} from './pricing-ui.js';
import {renderCharts, clearCharts} from './charts.js';
import {prepareTokenExamples} from './token-examples.js';
import {renderTokenExamples, clearTokenExamples} from './token-examples-view.js';

const form = document.getElementById('dr-calculator');
if (form) {
 const output = form.querySelector('#dr-calc-output');
 const download = form.querySelector('[data-action="download"]');
 const original = output.innerHTML;
 const currencySelect = form.elements.namedItem('currency');
 const modelName = form.elements.namedItem('modelName');
 const modelVersion = form.elements.namedItem('modelVersion');
 const chartHost = document.getElementById('dr-live-charts');
 const tokenExamplesHost = document.getElementById('dr-token-examples-live');
 let pricing;
 const number = (n, digits = 0) => new Intl.NumberFormat('en-US', {
  maximumFractionDigits: digits, minimumFractionDigits: digits,
 }).format(n);
 const tokens = n => new Intl.NumberFormat('en-US', {maximumFractionDigits: 6}).format(n / 1e6);
 const amount = n => new Intl.NumberFormat('en-US', {maximumFractionDigits: 2}).format(n);
 const pair = (en, zh) => `<span lang="en" data-lang="en">${en}</span><span lang="zh-CN" data-lang="zh">${zh}</span>`;
 // Explicit paths and scales keep the editable schema bounded; no dynamic object traversal.
 const scales = new Map([
  ['monthlyInputTokens', 1e6], ['monthlyOutputTokens', 1e6],
  ['capex', 1], ['monthlyOpex', 1], ['horizon', 1], ['annualDiscount', .01],
  ['inputTokensPerSecond', 1], ['outputTokensPerSecond', 1],
  ['peakHours', 1], ['offHours', 1], ['inputPeakShare', .01], ['outputPeakShare', .01],
  ['availability', .01], ['headroom', .01], ['apiCacheHit', .01], ['fallbackCacheHit', .01],
  ['apiAnnualChange', .01], ['terminalResidual', 1], ['inputAnnualGrowth', .01],
  ['outputAnnualGrowth', .01], ['opexAnnualGrowth', .01],
  ['perMillionLocalInputCost', 1], ['perMillionLocalOutputCost', 1],
  ['prices.peak.hit', 1], ['prices.peak.miss', 1], ['prices.peak.output', 1],
  ['prices.off.hit', 1], ['prices.off.miss', 1], ['prices.off.output', 1],
 ]);
 const pricePaths = new Map([
  ['prices.peak.hit', ['peak', 'hit']], ['prices.peak.miss', ['peak', 'miss']],
  ['prices.peak.output', ['peak', 'output']], ['prices.off.hit', ['off', 'hit']],
  ['prices.off.miss', ['off', 'miss']], ['prices.off.output', ['off', 'output']],
 ]);
 let snapshot = null;
 const account = (caption, headings, rows) => `<p class="dr-account-title">${caption}</p><table class="dr-account"><thead><tr>${headings.map(h => `<th scope="col">${h}</th>`).join('')}</tr></thead><tbody>${rows.map(row => `<tr><th scope="row">${row[0]}</th>${row.slice(1).map((cell, i) => `<td><span class="dr-account-label">${headings[i + 1]}</span>${cell}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
 function fail(en, zh) {
  snapshot = null;
  download.disabled = true;
  output.innerHTML = `<p class="dr-error">${pair(en, zh)}</p>`;
  if (chartHost) clearCharts(chartHost, {en, zh});
  clearTokenExamples(tokenExamplesHost, {en, zh});
 }
 function render() {
  const currency = currencySelect.value;
  if (!['CNY', 'USD'].includes(currency)) {
   fail('Choose CNY or USD and enter all amounts in that currency.', '请选择 CNY 或 USD，并以同一币种填写全部金额。');
   return;
  }
  for (const el of form.querySelectorAll('[data-currency]')) el.textContent = currency;
  if (!pricing.hasSelection()) {
   snapshot = null; download.disabled = true;
   const en = 'Choose a catalog model, then enter your usage, budget and measured throughput.';
   const zh = '选择目录中的模型，再填写用量、预算与实测吞吐。';
   output.innerHTML = `<p class="dr-result-detail">${pair(en, zh)}</p>`;
   if (chartHost) clearCharts(chartHost, {en, zh});
   clearTokenExamples(tokenExamplesHost);
   return;
  }
  pricing.enforceCachePolicy();
  const validName = modelName.validity.valid && Boolean(modelName.value.trim());
  modelName.setAttribute('aria-invalid', String(!validName));
  if (!validName || !modelVersion.validity.valid) {
   fail('Enter a model name for this scenario; version is optional.', '请为该场景填写模型名称；版本可选填。');
   return;
  }
  const metadata = pricing.metadata();
  const input = {capacityMultiplier: 1, prices: {peak: {}, off: {}}};
  let invalid = false;
  for (const [name, scale] of scales) {
   const field = form.elements.namedItem(name);
   const valid = field && field.validity.valid && Number.isFinite(field.valueAsNumber * scale);
   field?.setAttribute('aria-invalid', String(!valid));
   if (!valid) {invalid = true; continue;}
   const value = field.valueAsNumber * scale;
   const path = pricePaths.get(name);
   if (path) input.prices[path[0]][path[1]] = value;
   else input[name] = value;
  }
  if (invalid) {
   fail('Complete every field within its stated range; the horizon must be a whole number of months.', '请在允许范围内填全各项；评估期限须为整数月。');
   return;
  }
  if (input.peakHours + input.offHours <= 0 || input.peakHours + input.offHours > 744) {
   fail('Peak and off-peak hours must total more than 0 and no more than 744 hours per month.', '每月峰谷时长之和须大于 0，且不超过 744 小时。');
   return;
  }
  if ((input.peakHours === 0 && (input.inputPeakShare > 0 || input.outputPeakShare > 0)) ||
      (input.offHours === 0 && (input.inputPeakShare < 1 || input.outputPeakShare < 1))) {
   fail('A time bucket with zero hours cannot carry token demand. Adjust its hours or token shares.', '时长为零的时段不能承担 token 需求，请调整时长或输入、输出占比。');
   return;
  }
  try {
   const r = evaluateTokenScenario(input), m = r.firstMonth, horizon = r.parameters.horizon;
   const tokenExamples = prepareTokenExamples(r.parameters, metadata);
   output.innerHTML = original;
   const set = (key, value, negative = false) => {
    const e = output.querySelector(`[data-result="${key}"]`);
    if (!e) return;
    e.innerHTML = value;
    e.dataset.negative = String(negative);
   };
   set('api-label', pair(`All API · ${currency}/month`, `全用 API · ${currency}/月`));
   set('hybrid-label', pair(`Local + fallback API · ${currency}/month`, `本地运营 + 回退 API · ${currency}/月`));
   set('saving-label', pair(`Month 1 saving · ${currency}`, `首月节省 · ${currency}`));
   set('api', amount(m.fullApiCost));
   set('hybrid', amount(m.hybridOperatingCost));
   set('saving', amount(m.operatingSavings), m.operatingSavings < 0);
   set('monthly-story', pair(
    `Local running costs ${currency} ${amount(m.opex)}, plus ${amount(m.overflowApiCost)} for tokens still sent to the API. The difference from all-API is ${m.operatingSavings >= 0 ? 'a saving of' : 'an extra cost of'} ${currency} ${amount(Math.abs(m.operatingSavings))} in month 1. Initial investment is separate.`,
    `本地运行花费 ${amount(m.opex)} ${currency}，仍需发给 API 的 token 花费 ${amount(m.overflowApiCost)} ${currency}。相比全用 API，首月${m.operatingSavings >= 0 ? '省下' : '多花'} ${amount(Math.abs(m.operatingSavings))} ${currency}；购机初投另计。`
   ));
   set('capex-label', pair(`Initial investment · ${currency}`, `一次性投入 · ${currency}`));
   set('capex', amount(input.capex));
   set('npv-label', pair(`${horizon}-month NPV · ${currency}`, `${horizon} 月净现值 · ${currency}`));
   for (const el of output.querySelectorAll('[data-currency]')) el.textContent = currency;
   set('roi-label', pair(`${horizon}-month ROI`, `${horizon} 月投资回报率`));
   set('npv', number(r.NPV), r.NPV < 0);
   set('ceiling', number(r.maxInitialInvestment), r.maxInitialInvestment < 0);
   set('ceiling-note', pair(`Net resale is held at the entered amount. With zero resale, the ceiling is ${currency} ${number(r.maxInitialInvestmentWithoutResale)}.`, `净出售所得固定为所填金额；若出售所得为零，上限为 ${number(r.maxInitialInvestmentWithoutResale)} ${currency}。`));
   set('roi', r.capitalROI === null ? pair('Not applicable', '不适用') : number(r.capitalROI * 100, 1) + '%', r.capitalROI < 0);
   set('payback', input.capex === 0 ? pair('No initial capital', '无初投需回收') : r.simplePaybackMonth === null ? pair(`Not within ${horizon} months`, `${horizon} 月内未回本`) : pair(`Month ${r.simplePaybackMonth}`, `第 ${r.simplePaybackMonth} 月`), input.capex > 0 && r.simplePaybackMonth === null);
   const discounted = input.capex === 0 ? pair('no initial capital', '无初投需回收') : r.discountedPaybackMonth === null ? pair(`not within ${horizon} months`, `${horizon} 月内未回本`) : pair(`month ${r.discountedPaybackMonth}`, `第 ${r.discountedPaybackMonth} 月`);
   set('investment-story', pair(
    `Payback uses the cumulative monthly savings to recover ${currency} ${amount(input.capex)}; it excludes resale. NPV discounts ${horizon} months of savings plus the entered resale proceeds, then subtracts this investment.`,
    `回本看每月省下的钱何时累计抵回 ${amount(input.capex)} ${currency} 初投，不靠卖机器凑回本。净现值则把 ${horizon} 个月的节省和所填出售所得折算到今天，再减去初投。`
   ));
   set('detail', pair(`Month 1 local running cost = fixed costs ${currency} ${amount(m.fixedOpex)} + extra token costs ${amount(m.variableLocalCost)}.`, `首月本地运行费用 = 固定费用 ${amount(m.fixedOpex)} ${currency} + 额外 token 成本 ${amount(m.variableLocalCost)} ${currency}。`) + '<br>' + pair('Payback after discounting future savings: ', '考虑资金时间价值后的回本：') + discounted + '<br>' + pair(`${horizon}-month total spending · all API ${currency} ${number(r.fullApiTco)}; local deployment ${currency} ${number(r.hybridTco)} including purchase, running costs and fallback API, less resale.`, `${horizon} 个月总支出：全用 API 为 ${number(r.fullApiTco)} ${currency}；本地部署为 ${number(r.hybridTco)} ${currency}，包含购机、运行和回退 API，已扣除出售所得。`));
   const tokenHeadings = [pair('Token stream', 'Token 流'), pair('Local', '本地承接'), pair('Fallback API', '回退 API')];
   set('tokens', account(pair('Month 1 token routing · millions', '首月 token 分流 · 百万'), tokenHeadings, [
    [pair(`Input · total ${tokens(m.demandInputTokens)}`, `输入 · 共 ${tokens(m.demandInputTokens)}`), tokens(m.localInputTokens), tokens(m.overflowInputTokens)],
    [pair(`Output · total ${tokens(m.demandOutputTokens)}`, `输出 · 共 ${tokens(m.demandOutputTokens)}`), tokens(m.localOutputTokens), tokens(m.overflowOutputTokens)],
   ]));
   const billed = b => `${number(b.amount, 2)}<span class="dr-account-small">${pair(`${tokens(b.tokens)} million tokens`, `${tokens(b.tokens)} 百万 token`)}</span>`;
   set('billing', account(pair(`Month 1 API bill · ${currency}`, `首月 API 账单 · ${currency}`), [pair('Billed category', '计费项目'), pair('All-API baseline', '全 API 基线'), pair('Fallback only', '仅回退 API')], [
    [pair('Cached input', '缓存命中输入'), billed(m.baselineBillBreakdown.hit), billed(m.fallbackBillBreakdown.hit)],
    [pair('Uncached input', '缓存未命中输入'), billed(m.baselineBillBreakdown.miss), billed(m.fallbackBillBreakdown.miss)],
    [pair('Output', '输出'), billed(m.baselineBillBreakdown.output), billed(m.fallbackBillBreakdown.output)],
    [pair('Total', '合计'), number(m.fullApiCost, 2), number(m.overflowApiCost, 2)],
   ]));
   let en = r.NPV > 0 ? 'Under these assumptions, local deployment costs less over the chosen period after discounting. Validate throughput and lower-demand scenarios before purchasing.' : 'Under these assumptions, buying hardware does not beat all-API over the chosen period after discounting.';
   let zh = r.NPV > 0 ? '按当前假设，考虑时间价值后，本地部署在评估期内更省钱。采购前还要验证实测吞吐与用量下降的情况。' : '按当前假设，考虑时间价值后，购买硬件在评估期内没有比全用 API 更省钱。';
   if (input.monthlyInputTokens === 0 && input.monthlyOutputTokens === 0) {en = 'Enter monthly token usage to calculate the two bills.'; zh = '先填写每月 token 用量，才能比较两笔账。';}
   else if ((input.monthlyInputTokens > 0 && input.inputTokensPerSecond === 0) || (input.monthlyOutputTokens > 0 && input.outputTokensPerSecond === 0)) {en = 'Enter measured throughput for the token streams you serve locally. The current calculation sends any unserved workload to the API.'; zh = '请填写本地承接所需的实测吞吐；当前无法承接的负载全部按回退 API 计费。';}
   set('verdict', pair(en, zh));
   const caveats = [];
   if (r.simplePaybackReversalMonths.length) caveats.push(pair(`Cumulative operating savings fall below the investment again in month ${r.simplePaybackReversalMonths.join(', ')}.`, `累计运营节省在第 ${r.simplePaybackReversalMonths.join('、')} 月再次低于初投。`));
   if (input.terminalResidual > input.capex) caveats.push(pair('Entered resale proceeds exceed the initial investment; verify this assumption.', '所填出售所得超过初投，请核实这一假设。'));
   if (m.cachePenalty !== 0) caveats.push(pair(`Different cache-hit rates change the month 1 fallback bill by ${currency} ${amount(m.cachePenalty)} versus using the all-API cache rate for those same tokens.`, `相对回退 token 沿用全 API 命中率，两项命中率的差别让首月回退账单改变 ${amount(m.cachePenalty)} ${currency}。`));
   set('caveats', caveats.join('<br>'));
   snapshot = {schemaVersion: 2, model: 'open-model-token-investment', metadata, currency, priceUnit: `${currency} per million tokens`, parameters: r.parameters, results: r, tokenExamples};
   download.disabled = false;
   if (chartHost) renderCharts(chartHost, r, metadata);
   renderTokenExamples(tokenExamplesHost, tokenExamples);
  } catch {
   fail('These inputs cannot be calculated. Check the values and reset if needed.', '这组输入无法计算，请核对数值，或重置表单。');
  }
 }
 function resetParameters() {
  for (const [name, scale] of scales) {
   const path = pricePaths.get(name);
   const value = path ? customDefaults.prices[path[0]][path[1]] : name === 'fallbackCacheHit' ? (customDefaults.fallbackCacheHit ?? customDefaults.apiCacheHit) : customDefaults[name];
   form.elements.namedItem(name).value = String(value / scale);
   form.elements.namedItem(name).removeAttribute('aria-invalid');
  }
 }
 pricing = createPricingUI(form, {onChange: render, resetParameters});
 form.hidden = false;
 form.addEventListener('input', event => {
  if (event.target.tagName === 'SELECT' || pricing.isManagedInput(event.target.name)) return;
  render();
 });
 form.addEventListener('reset', event => {event.preventDefault(); pricing.reset();});
 form.addEventListener('submit', e => e.preventDefault());
 download.addEventListener('click', () => {
  if (!snapshot) return;
  const timestamp = new Date().toISOString();
  const blob = new Blob([JSON.stringify({...snapshot, exportedAt: timestamp}, null, 2) + '\n'], {type: 'application/json'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `token-investment-${timestamp.slice(0, 10)}.json`;
  a.hidden = true;
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
 });
 pricing.reset();
 pricing.load();
}

const toc = document.querySelector('.article-page:has(.dr-article) .article-toc details');
if (toc) {
 const narrow = matchMedia('(max-width:1023px)');
 const adapt = () => {toc.open = !narrow.matches;};
 adapt();
 narrow.addEventListener('change', adapt);
}
