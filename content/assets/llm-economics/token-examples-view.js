/**
 * Accessible, script-free markup for controlled API token-cost comparisons.
 * The static article and live calculator use exactly the same view.
 */
const escape = value => String(value ?? '').replace(/[&<>"']/g, char => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[char]));
const number = value => new Intl.NumberFormat('en-US', {maximumSignificantDigits: 6}).format(value);
const tokens = value => `${number(value / 1e6)}M`;
const pair = (en, zh, cls = '') => `<span${cls ? ` class="${escape(cls)}"` : ''}><span lang="en" data-lang="en">${escape(en)}</span><span lang="zh-CN" data-lang="zh">${escape(zh)}</span></span>`;
const note = (en, zh, cls = '') => `<p class="te-note${cls ? ` ${escape(cls)}` : ''}">${pair(en, zh)}</p>`;
const money = (value, currency) => `${currency} ${number(value)}`;
const parts = [
  {key: 'miss', en: 'Ordinary input', zh: '普通输入', cls: 'te-miss'},
  {key: 'hit', en: 'Cached input', zh: '缓存输入', cls: 'te-hit'},
  {key: 'output', en: 'Output', zh: '输出', cls: 'te-output'},
];

function scaleMaximum(values) {
  const max = Math.max(0, ...values);
  if (!max) return 0;
  const power = 10 ** Math.floor(Math.log10(max));
  return [1, 2, 2.5, 5, 8, 10].find(step => step * power >= max) * power;
}

function caseHTML(bill, currency, maximum, cached) {
  if (!bill) return `<div class="te-case te-unavailable" data-cache-hit="0.8" data-cache-unavailable="true">${pair('80% input cached: no cache-read price listed', '输入缓存 80%：未提供缓存报价')}</div>`;
  const percent = number(bill.cacheHit * 100);
  const label = cached
    ? [`${percent}% input cached · ${tokens(bill.hit.tokens)}`, `输入缓存 ${percent}% · ${tokens(bill.hit.tokens)}`]
    : ['0% input cached', '输入缓存 0%'];
  const description = parts.map(part => `${part.en} / ${part.zh}: ${tokens(bill[part.key].tokens)} tokens, ${money(bill[part.key].amount, currency)}`).join('; ');
  const segments = parts.map(part => {
    const item = bill[part.key], width = maximum > 0 ? item.amount / maximum * 100 : 0;
    return `<span class="te-segment ${part.cls}" style="width:${Number.isFinite(width) ? width : 0}%" data-component="${part.key}" data-tokens="${escape(item.tokens)}" data-amount="${escape(item.amount)}" title="${escape(`${part.en} / ${part.zh}: ${tokens(item.tokens)} tokens · ${money(item.amount, currency)}`)}"></span>`;
  }).join('');
  return `<div class="te-case" data-cache-hit="${escape(bill.cacheHit)}" data-bill="${escape(bill.total)}">
    <div class="te-case-heading">${pair(...label, 'te-case-label')}<strong class="te-bill">${escape(money(bill.total, currency))}</strong></div>
    <div class="te-bar" role="img" aria-label="${escape(description)}"><span class="te-segments" aria-hidden="true">${segments}</span></div>
  </div>`;
}

function matrixHTML(data) {
  const headers = data.scaleMatrix.totals.map(total => `<th scope="col">${escape(tokens(total))}</th>`).join('');
  const rows = data.scaleMatrix.rows.map(row => {
    const ratio = `${number(row.inputShare * 100)} / ${number((1 - row.inputShare) * 100)}`;
    const cells = row.cells.map(cell => {
      const exact = number(cell.total), label = `${money(cell.total, data.currency)}; ${tokens(cell.inputTokens)} input / 输入 + ${tokens(cell.outputTokens)} output / 输出`;
      return `<td data-ratio="${escape(row.id)}" data-total-tokens="${escape(cell.totalTokens)}" data-bill="${escape(cell.total)}" title="${escape(label)}"><span class="te-cell-total" aria-hidden="true">${escape(tokens(cell.totalTokens))}</span><span class="te-matrix-value">${escape(exact)}</span></td>`;
    }).join('');
    return `<tr><th scope="row">${escape(ratio)}<span class="te-mobile-ratio">${pair('input / output %', '输入 / 输出 %')}</span></th>${cells}</tr>`;
  }).join('');
  return `<div class="te-matrix">
    <table>
      <caption>${pair('Change the total: API bill with no cached input', '改变总量：无输入缓存时的 API 账单')}<span class="te-matrix-unit">${escape(data.currency)} · ${pair('columns = total tokens; prices stay fixed', '各列为总 token 数，单价固定')}</span></caption>
      <thead><tr><th scope="col">${pair('In / out %', '输入 / 输出 %')}</th>${headers}</tr></thead>
      <tbody>${rows}</tbody>
    </table>
  </div>`;
}

/** Returns complete safe HTML; no browser globals, scripts or event handlers. */
export function tokenExamplesHTML(data, options = {}) {
  if (!data) return `<div class="token-examples te-empty">${note('Choose a model to compare token mixes at its current prices.', '选择模型后，比较不同 token 组合的 API 账单。')}</div>`;
  const mode = options.mode === 'teaching' || (options.mode !== 'live' && data.source?.kind === 'teaching') ? 'teaching' : 'live';
  const showMatrix = options.includeMatrix ?? mode === 'live';
  const source = data.source || {}, identity = [source.provider, source.modelName].filter(Boolean).join(' · ');
  const maximum = scaleMaximum(data.ratios.flatMap(row => [row.uncached.total, ...(row.cached ? [row.cached.total] : [])]));
  const title = mode === 'teaching'
    ? ['The same token total can produce different bills', '同样的 token 总量，账单也会不同']
    : ['Compare token mixes at the selected prices', '用当前单价比较 token 组合'];
  const totalLabel = mode === 'teaching'
    ? ['Teaching total', '教学总量']
    : data.usesEnteredTotal ? ['Entered total', '当前填写总量'] : ['Reference total', '参考总量'];
  const rates = parts.map(part => {
    const value = data.effectivePrices[part.key];
    return `<div><dt>${pair(part.en, part.zh)}</dt><dd>${value === null ? pair('Not listed', '未提供') : escape(number(value))}</dd></div>`;
  }).join('');
  const legend = parts.map(part => `<span class="te-legend-item"><i class="te-key ${part.cls}" aria-hidden="true"></i>${pair(`${part.en} cost`, `${part.zh}费用`)}</span>`).join('');
  const groups = data.ratios.map(row => `<div class="te-group" data-ratio="${escape(row.id)}" data-input-tokens="${escape(row.inputTokens)}" data-output-tokens="${escape(row.outputTokens)}">
    <div class="te-group-heading"><strong>${pair(row.label.en, row.label.zh)}</strong><span class="te-token-count">${pair(`${tokens(row.inputTokens)} input + ${tokens(row.outputTokens)} output`, `${tokens(row.inputTokens)} 输入＋${tokens(row.outputTokens)} 输出`)}</span></div>
    ${caseHTML(row.uncached, data.currency, maximum, false)}
    ${caseHTML(row.cached, data.currency, maximum, true)}
  </div>`).join('');
  const context = mode === 'teaching'
    ? note('Illustrative prices, not a provider quote. M = one million tokens.', '以下是教学假设单价，不是供应商报价。M 表示百万 token。')
    : note(`Month 1 prices weighted by the current input and output time shares. M = one million tokens.${data.usesEnteredTotal ? '' : ' The 100M reference is not a workload forecast.'}`, `按首月单价及当前输入、输出各自的峰谷占比加权。M 表示百万 token。${data.usesEnteredTotal ? '' : '100M 参考量不代表实际用量预测。'}`);
  return `<figure class="token-examples" data-token-examples="${mode}" data-total-tokens="${escape(data.totalTokens)}" data-currency="${escape(data.currency)}" data-provider="${escape(source.providerId || source.provider || '')}" data-model="${escape(source.modelName || '')}">
    <div class="te-heading"><h3>${pair(...title)}</h3><div class="te-total"><span class="te-total-label">${pair(...totalLabel)}</span><strong>${escape(tokens(data.totalTokens))}</strong><span>tokens</span></div></div>
    ${identity && mode === 'live' ? `<p class="te-identity">${escape(identity)}${source.modified ? ` · ${pair('edited inputs', '已手动修改')}` : ''}</p>` : ''}
    ${context}
    <div class="te-rates-heading">${pair('Prices per million tokens', '每百万 token 单价')}<span>${escape(data.currency)}</span></div>
    <dl class="te-rates">${rates}</dl>
    <div class="te-legend">${legend}</div>
    <div class="te-scale"><span>0</span><span>${escape(money(maximum, data.currency))}</span></div>
    <div class="te-groups">${groups}</div>
    ${note('Every bar starts at zero and uses the same money scale. Cache percentages refer only to input tokens.', '所有横条从零起算，共用金额标尺；缓存比例只针对输入 token。', 'te-scale-note')}
    ${showMatrix ? matrixHTML(data) : ''}
    <figcaption>${pair('API bills only. Equal token totals do not mean equal tasks or quality; these are not local-deployment costs.', '这里只比较 API 账单。相同 token 总量不代表相同任务或质量，也不代表本地部署成本。')}</figcaption>
  </figure>`.replace(/[ \t]+$/gm, '');
}

export function clearTokenExamples(host, message = {en: 'Choose a model to compare token mixes at its current prices.', zh: '选择模型后，比较不同 token 组合的 API 账单。'}) {
  if (!host) return;
  const text = typeof message === 'object' && message !== null ? message : {en: String(message), zh: String(message)};
  host.innerHTML = `<div class="token-examples te-empty">${note(text.en || text.zh, text.zh || text.en)}</div>`;
}

export function renderTokenExamples(host, data) {
  if (!host) return;
  if (!data) {clearTokenExamples(host); return;}
  host.innerHTML = tokenExamplesHTML(data, {mode: 'live'});
}
