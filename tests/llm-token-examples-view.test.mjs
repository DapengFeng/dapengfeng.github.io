import test from 'node:test';
import assert from 'node:assert/strict';
import {load} from 'cheerio';
import {teachingExample, prepareTokenExamples} from '../content/assets/llm-economics/token-examples.js';
import {tokenExamplesHTML, renderTokenExamples, clearTokenExamples} from '../content/assets/llm-economics/token-examples-view.js';
import {customDefaults} from '../content/assets/llm-economics/defaults.js';

const close = (a, b) => assert.ok(Math.abs(a - b) <= 1e-10 * Math.max(1, Math.abs(b)), `${a} != ${b}`);
const parameters = overrides => ({
 ...customDefaults,
 prices: {peak: {hit: .5, miss: 2, output: 8}, off: {hit: .5, miss: 2, output: 8}},
 ...overrides,
});
const metadata = {
 name: 'Example model', currency: 'USD', modified: true,
 pricing: {
  provider: 'Example provider', providerId: 'provider-one', modelId: 'model-one',
  id: 'provider-one:model-one', sourceUrl: 'https://models.dev/api.json',
  fetchedAt: '2026-10-20T00:00:00Z', cacheSupported: true,
  sourceLabel: {en: 'Catalog reference prices', zh: '目录参考价'},
 },
};
const parse = (data, options) => load(tokenExamplesHTML(data, options));

test('the teaching figure shows all six bills and their exact token/amount decomposition', () => {
 const $ = parse(teachingExample);
 assert.equal($('figure[data-token-examples=teaching]').length, 1);
 assert.equal($('.te-group').length, 3);
 assert.equal($('.te-case[data-bill]').length, 6);
 assert.deepEqual($('.te-case[data-bill]').map((_, node) => Number($(node).attr('data-bill'))).get(), [260, 152, 500, 440, 740, 728]);
 for (const ratio of teachingExample.ratios) {
  const group = $(`.te-group[data-ratio="${ratio.id}"]`);
  close(Number(group.attr('data-input-tokens')), ratio.inputTokens);
  close(Number(group.attr('data-output-tokens')), ratio.outputTokens);
  for (const bill of [ratio.uncached, ratio.cached]) {
   const row = group.find(`.te-case[data-cache-hit="${bill.cacheHit}"]`);
   assert.equal(row.find('.te-bill').text(), `CNY ${bill.total}`);
   assert.equal(row.find('.te-segment').length, 3);
   let total = 0;
   for (const part of ['hit', 'miss', 'output']) {
    const segment = row.find(`[data-component="${part}"]`);
    close(Number(segment.attr('data-tokens')), bill[part].tokens);
    close(Number(segment.attr('data-amount')), bill[part].amount);
    total += Number(segment.attr('data-amount'));
   }
   close(total, Number(row.attr('data-bill')));
   assert.ok(row.find('.te-bar[role=img]').attr('aria-label').includes('Ordinary input'));
  }
 }
 assert.equal($('.te-matrix').length, 0, 'the compact static figure omits the optional volume matrix');
 assert.match($('figcaption').text(), /相同 token 总量不代表相同任务或质量/);
});

test('every component uses the same money scale rather than normalizing each bar', () => {
 const $ = parse(teachingExample);
 assert.equal($('.te-scale').children().last().text(), 'CNY 800');
 for (const row of $('.te-case[data-bill]').toArray()) {
  let sumWidth = 0;
  for (const part of $(row).find('.te-segment').toArray()) {
   const match = /^width:([\d.eE+-]+)%$/.exec($(part).attr('style'));
   assert.ok(match, 'segment width is a finite numeric percentage');
   const width = Number(match[1]), amount = Number($(part).attr('data-amount'));
   close(width, amount / 800 * 100);
   sumWidth += width;
  }
  close(sumWidth, Number($(row).attr('data-bill')) / 800 * 100);
 }
 const outputHeavyWidth = $('.te-group[data-ratio=output-heavy] .te-case[data-cache-hit="0"] .te-segment')
  .toArray().reduce((sum, part) => sum + Number(/^width:([\d.eE+-]+)%$/.exec($(part).attr('style'))[1]), 0);
 close(outputHeavyWidth, 92.5);
});

test('a missing cache-read tariff has no fabricated cached bill, segment or zero-price label', () => {
 const data = prepareTokenExamples(parameters(), {...metadata, pricing: {...metadata.pricing, cacheSupported: false}});
 const $ = parse(data, {mode: 'live'});
 assert.equal($('.te-case[data-bill]').length, 3);
 assert.equal($('.te-case[data-cache-unavailable=true]').length, 3);
 for (const unavailable of $('.te-case[data-cache-unavailable=true]').toArray()) {
  assert.equal($(unavailable).attr('data-bill'), undefined);
  assert.equal($(unavailable).find('.te-bar, .te-segment, .te-bill').length, 0);
  assert.match($(unavailable).text(), /未提供缓存报价/);
 }
 assert.match($('.te-rates > div').eq(1).text(), /Not listed/);
 assert.doesNotMatch($('.te-rates > div').eq(1).find('dd').text(), /^0(?:\.0+)?$/);
 assert.equal($('.te-matrix td[data-bill]').length, 9, 'ordinary-input volume comparisons remain available');
});

test('the live matrix has nine correctly mapped bills and preserves the selected identity', () => {
 const data = prepareTokenExamples(parameters({monthlyInputTokens: 240e6, monthlyOutputTokens: 60e6}), metadata);
 const $ = parse(data, {mode: 'live'}), figure = $('figure');
 assert.equal(figure.attr('data-token-examples'), 'live');
 assert.equal(figure.attr('data-total-tokens'), '300000000');
 assert.equal(figure.attr('data-currency'), 'USD');
 assert.equal(figure.attr('data-provider'), 'provider-one');
 assert.equal(figure.attr('data-model'), 'Example model');
 assert.match($('.te-identity').text(), /Example provider · Example model/);
 assert.match($('.te-identity').text(), /已手动修改/);
 assert.match($('.te-total').text(), /Entered total/);
 assert.match($('.te-total').text(), /300M/);
 assert.equal($('.te-matrix th[scope=col]').length, 4);
 assert.equal($('.te-matrix tbody th[scope=row]').length, 3);
 assert.equal($('.te-matrix td').length, 9);
 for (const ratio of data.scaleMatrix.rows) for (const bill of ratio.cells) {
  const cell = $(`.te-matrix td[data-ratio="${ratio.id}"][data-total-tokens="${bill.totalTokens}"]`);
  assert.equal(cell.length, 1);
  close(Number(cell.attr('data-bill')), bill.total);
  assert.ok(cell.attr('title').includes('USD'));
  assert.ok(cell.attr('title').includes('input / 输入'));
 }
 assert.match($('.te-matrix caption').text(), /单价固定/);
 assert.match($('.te-note').text(), /输入、输出各自的峰谷占比加权/);
});

test('model, provider, currency, labels and attributes are escaped rather than executed', () => {
 const payload = `"><img src=x onerror="alert('test')"><script>evil()</script>&`;
 const data = prepareTokenExamples(parameters(), {
  ...metadata, name: payload, currency: payload,
  pricing: {...metadata.pricing, provider: payload, providerId: payload},
 });
 data.ratios[0].label = {en: payload, zh: payload};
 data.ratios[0].id = payload;
 data.scaleMatrix.rows[0].id = payload;
 const html = tokenExamplesHTML(data, {mode: 'live'}), $ = load(html);
 assert.equal($('script, img, [onerror], [onload]').length, 0);
 assert.equal($('figure').attr('data-model'), payload);
 assert.equal($('figure').attr('data-provider'), payload);
 assert.equal($('figure').attr('data-currency'), payload);
 assert.equal($('.te-group').first().attr('data-ratio'), payload);
 assert.equal($('.te-matrix td').first().attr('data-ratio'), payload);
 assert.ok($('.te-identity').text().includes(payload));
 assert.match(html, /&lt;script&gt;evil\(\)&lt;\/script&gt;/);
});

test('zero prices and reference totals remain honest, finite and accessible', () => {
 const zero = {hit: 0, miss: 0, output: 0};
 const $ = parse(prepareTokenExamples(parameters({prices: {peak: zero, off: zero}}), metadata));
 assert.match($('.te-total').text(), /Reference total/);
 assert.match($('.te-note').text(), /100M 参考量不代表实际用量预测/);
 assert.equal($('.te-scale').children().last().text(), 'USD 0');
 for (const part of $('.te-segment').toArray()) assert.equal($(part).attr('style'), 'width:0%');
 assert.doesNotMatch($.html(), /NaN|Infinity/);
 assert.equal($('.te-bar[role=img][aria-label]').length, 6);
});

test('render and clear replace stale examples using the same escaped HTML function', () => {
 const host = {innerHTML: '<p>previous result</p>'};
 renderTokenExamples(host, teachingExample);
 assert.equal(host.innerHTML, tokenExamplesHTML(teachingExample, {mode: 'live'}));
 clearTokenExamples(host, {en: '<img onerror=oops()>', zh: '填写 <模型>'});
 let $ = load(host.innerHTML);
 assert.equal($('figure, .te-case, img').length, 0);
 assert.ok($('.te-empty').text().includes('<img onerror=oops()>'));
 renderTokenExamples(host, null);
 $ = load(host.innerHTML);
 assert.equal($('.te-empty').length, 1);
 assert.match($('.te-empty').text(), /选择模型/);
 assert.doesNotThrow(() => renderTokenExamples(null, teachingExample));
 assert.doesNotThrow(() => clearTokenExamples(null));
});
