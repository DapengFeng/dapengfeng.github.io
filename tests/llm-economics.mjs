import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {chromium} from '@playwright/test';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
import {evaluateTokenScenario} from '../content/assets/llm-economics/token-model.js';
import {prepareTokenExamples, teachingExample} from '../content/assets/llm-economics/token-examples.js';
import {tokenExamplesHTML} from '../content/assets/llm-economics/token-examples-view.js';

const articleSource = await fs.readFile('content/posts/llm-deployment-roi.html', 'utf8');
assert.equal(articleSource.split('<!-- token-examples:start -->')[1].split('<!-- token-examples:end -->')[0].trim(),
  tokenExamplesHTML(teachingExample, {mode: 'teaching'}), 'static figure matches its calculation');

const server = serve(4247);
const browser = await chromium.launch({executablePath: chromiumExecutable(), args: ['--no-sandbox']});
const url = 'http://localhost:4247/blog/llm-deployment-roi.html';
const catalog = JSON.parse(await fs.readFile('dist/assets/content/llm-economics/model-prices.json', 'utf8'));
const first = catalog.quotes.find(q => q.cacheSupported && q.prices.peak.miss > 0 && q.prices.peak.output > 0);
const second = catalog.quotes.find(q => q.providerId !== first?.providerId && q.modelId !== first?.modelId && q.cacheSupported && q.prices.peak.miss > 0 && q.prices.peak.output > 0);
const noCache = catalog.quotes.find(q => q.currency === 'USD' && !q.cacheSupported && q.prices.peak.output > 0);
assert.ok(first && second && noCache);
const errors = [];
const screenshotStyle = '.lab-header, .article-action-dock, .site-skip, .article-progress { visibility: hidden !important; }';
const near = (a, b) => assert.ok(Math.abs(a - b) <= 1e-9 * Math.max(1, Math.abs(b)), `${a} != ${b}`);
await fs.mkdir('test-results/economics', {recursive: true});
try {
  const context = await browser.newContext({viewport: {width: 1440, height: 1000}});
  await context.addInitScript(() => localStorage.setItem('feng-language', 'zh'));
  await context.route('**/*', route => new URL(route.request().url()).hostname === 'localhost'
    ? route.continue() : route.fulfill({status: 200, body: ''}));
  const page = await context.newPage();
  page.on('pageerror', e => errors.push(e.message));
  await page.goto('http://localhost:4247/blog/deepseek-v4-1-flash-roi.html?via=legacy#token-cost-examples');
  await page.waitForURL(`${url}?via=legacy#token-cost-examples`);
  await page.goto(url);
  await page.waitForFunction(() => document.querySelector('#dr-presetId').options.length > 2);
  const field = name => page.locator(`[name="${name}"]`);
  const charts = page.locator('#dr-live-charts');
  const examples = page.locator('#dr-token-examples-live');
  const teaching = page.locator('#dr-token-examples-teaching');
  const select = field('presetId');
  const fill = async (name, value) => field(name).fill(String(value));
  const allDetails = () => page.locator('#dr-calculator details').evaluateAll(nodes => nodes.forEach(n => n.open = true));
  const exportScenario = async () => {
    const pending = page.waitForEvent('download');
    await page.locator('[data-action=download]').click();
    const download = await pending;
    return JSON.parse(await fs.readFile(await download.path(), 'utf8'));
  };
  const checkChartValues = async scenario => {
    assert.equal(await charts.getAttribute('data-model'), scenario.metadata.name);
    assert.equal(await charts.getAttribute('data-currency'), scenario.currency);
    for (const rate of ['hit', 'miss', 'output']) {
      near(Number(await charts.locator(`[data-price="peak.${rate}"]`).getAttribute('data-value')), scenario.parameters.prices.peak[rate]);
    }
    near(Number(await charts.locator('[data-component-total=fixed]').getAttribute('data-value')), scenario.results.firstMonth.fixedOpex);
    near(Number(await charts.locator('[data-series=current-demand]').getAttribute('data-value')), scenario.results.NPV);
    const ordinary = JSON.parse(await charts.locator('[data-series=ordinary]').getAttribute('data-values'));
    near(ordinary.at(-1)[1], scenario.results.months.at(-1).cumulativeOperatingNet);
    assert.equal(ordinary.at(-1)[0], scenario.parameters.horizon);
  };
  const checkTokenExamples = async scenario => {
    const expected = prepareTokenExamples(scenario.parameters, scenario.metadata);
    assert.deepEqual(scenario.tokenExamples, expected, 'exported comparisons use the same parameters');
    const figure = examples.locator('[data-token-examples=live]');
    near(Number(await figure.getAttribute('data-total-tokens')), expected.totalTokens);
    assert.equal(await figure.getAttribute('data-model'), scenario.metadata.name);
    assert.equal(await figure.getAttribute('data-currency'), scenario.currency);
    for (const row of expected.ratios) {
      const group = examples.locator(`.te-group[data-ratio="${row.id}"]`);
      near(Number(await group.getAttribute('data-input-tokens')), row.inputTokens);
      near(Number(await group.getAttribute('data-output-tokens')), row.outputTokens);
      for (const [hit, bill] of [[0, row.uncached], [.8, row.cached]]) {
        const rendered = group.locator(`[data-cache-hit="${hit}"]`);
        if (!bill) {
          assert.equal(await rendered.getAttribute('data-cache-unavailable'), 'true');
          assert.equal(await rendered.locator('.te-bar').count(), 0);
          continue;
        }
        near(Number(await rendered.getAttribute('data-bill')), bill.total);
        for (const key of ['hit', 'miss', 'output']) {
          near(Number(await rendered.locator(`[data-component="${key}"]`).getAttribute('data-amount')), bill[key].amount);
        }
      }
    }
    assert.equal(await examples.locator('td[data-bill]').count(), 9);
    for (const row of expected.scaleMatrix.rows) for (const cell of row.cells) {
      near(Number(await examples.locator(`td[data-ratio="${row.id}"][data-total-tokens="${cell.totalTokens}"]`).getAttribute('data-bill')), cell.total);
    }
  };
  assert.equal(await select.inputValue(), '', 'no fixed provider preset');
  assert.equal(await field('currency').inputValue(), 'USD', 'models.dev prices are USD; no implicit conversion');
  assert.equal(await field('monthlyInputTokens').inputValue(), '0');
  assert.equal(await field('capex').inputValue(), '0');
  assert.equal(await teaching.locator('[data-bill]').count(), 6);
  assert.equal(await examples.locator('[data-token-examples]').count(), 0);
  await field('modelProvider').selectOption(first.providerId);
  await field('modelSearch').fill(first.modelId);
  assert.ok(await select.locator('option').count() <= catalog.quotes.filter(q => q.providerId === first.providerId).length + 2);
  await select.selectOption(first.id);
  await field('modelSearch').fill('no-such-model-qa-check');
  assert.equal(await select.inputValue(), first.id, 'searching does not alter the current calculation');
  await field('modelSearch').fill('');
  await field('modelProvider').selectOption('');
  await allDetails();
  for (const bucket of ['peak', 'off']) for (const type of ['hit', 'miss', 'output']) {
    near(Number(await field(`prices.${bucket}.${type}`).inputValue()), first.prices[bucket][type]);
  }
  assert.ok(await charts.locator('[data-chart=prices] [data-price]').count() > 0, 'rate plot works without invented workload');
  await checkTokenExamples(await exportScenario());
  assert.match(await examples.textContent(), /100M reference/);
  assert.ok(await charts.locator('.lc-track').first().evaluate(el => el.getBoundingClientRect().height >= 8), 'published chart stylesheet is applied');
  assert.match(await charts.textContent(), new RegExp(first.modelName.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')));
  // Synthetic test inputs are supplied only by the test, never by the article.
  for (const [name, value] of Object.entries({
    monthlyInputTokens: 200, monthlyOutputTokens: 50, capex: 10000, monthlyOpex: 100,
    inputTokensPerSecond: 1000, outputTokensPerSecond: 250, availability: 98,
    headroom: 80, horizon: 24, annualDiscount: 10, terminalResidual: 1000,
    apiCacheHit: 40, fallbackCacheHit: 0,
  })) await fill(name, value);
  const oldPlot = await charts.innerHTML();
  const exported = await exportScenario();
  near(exported.results.NPV, evaluateTokenScenario(exported.parameters).NPV);
  near(exported.parameters.monthlyInputTokens, 200e6);
  near(exported.parameters.monthlyOutputTokens, 50e6);
  near(exported.parameters.capex, 10000);
  near(exported.parameters.monthlyOpex, 100);
  const resultOrder = await page.locator('#dr-calc-output [data-result]').evaluateAll(nodes => nodes.map(n => n.dataset.result));
  assert.ok(resultOrder.indexOf('api') < resultOrder.indexOf('payback'), 'monthly bill precedes investment metrics');
  assert.ok(resultOrder.indexOf('hybrid') < resultOrder.indexOf('payback'));
  assert.ok(resultOrder.indexOf('saving') < resultOrder.indexOf('payback'));
  assert.equal(await page.locator('.dr-bill-details').getAttribute('open'), null, 'detailed billing is initially collapsed');
  const displayed = async key => Number((await page.locator(`[data-result="${key}"]`).textContent()).replaceAll(',', ''));
  assert.ok(Math.abs(await displayed('api') - exported.results.firstMonth.fullApiCost) <= .0051);
  assert.ok(Math.abs(await displayed('hybrid') - exported.results.firstMonth.hybridOperatingCost) <= .0051);
  assert.ok(Math.abs(await displayed('saving') - exported.results.firstMonth.operatingSavings) <= .0051);
  assert.ok(Math.abs(await displayed('api') - await displayed('hybrid') - await displayed('saving')) <= .011, 'displayed monthly bills reconcile within cent rounding');
  assert.equal(exported.metadata.pricing.id, first.id);
  assert.equal(exported.metadata.pricing.overrides.prices, false);
  assert.equal(exported.metadata.pricing.fetchedAt, first.fetchedAt);
  assert.equal(await charts.locator('[data-chart]').count(), 4, 'price, cost, cumulative and sensitivity plots render');
  await checkChartValues(exported);
  await checkTokenExamples(exported);
  const oldExamples = await examples.innerHTML();
  await select.selectOption(second.id);
  assert.equal(await field('inputTokensPerSecond').inputValue(), '0', 'new model cannot inherit throughput');
  assert.equal(await field('outputTokensPerSecond').inputValue(), '0');
  assert.equal(await field('capex').inputValue(), '10000', 'budget remains a user input');
  assert.notEqual(await charts.innerHTML(), oldPlot, 'model selection updates plots');
  assert.notEqual(await examples.innerHTML(), oldExamples, 'model selection updates token comparisons');
  await checkTokenExamples(await exportScenario());
  await fill('inputTokensPerSecond', 1200);
  await fill('outputTokensPerSecond', 300);
  const secondPlot = await charts.innerHTML();
  await fill('prices.peak.output', second.prices.peak.output + 5);
  assert.notEqual(await charts.innerHTML(), secondPlot, 'manual price updates plots');
  const edited = await exportScenario();
  assert.equal(edited.metadata.pricing.overrides.prices, true);
  assert.equal(edited.metadata.pricing.prices.peak.output, second.prices.peak.output, 'source quote preserved');
  near(edited.results.NPV, evaluateTokenScenario(edited.parameters).NPV);
  await checkChartValues(edited);
  await checkTokenExamples(edited);
  await examples.screenshot({path: 'test-results/economics/token-mixes-cache-desktop-zh.png', style: screenshotStyle});
  await field('modelName').evaluate(el => {const details = el.closest('details'); if (details) details.open = true;});
  await fill('modelName', 'Manually renamed test model');
  assert.equal(await field('inputTokensPerSecond').inputValue(), '0', 'manual identity changes invalidate the throughput measurement');
  assert.equal(await field('outputTokensPerSecond').inputValue(), '0');
  const renamed = await exportScenario();
  assert.equal(renamed.metadata.pricing.overrides.identity, true);
  assert.equal(renamed.metadata.pricing.id, second.id, 'source identity remains auditable after manual changes');
  await fill('monthlyInputTokens', '');
  assert.equal(await page.locator('[data-action=download]').isDisabled(), true);
  assert.equal(await charts.locator('svg').count(), 0, 'invalid inputs clear previous plots');
  assert.equal(await examples.locator('[data-bill]').count(), 0, 'invalid inputs clear token comparisons');
  await fill('monthlyInputTokens', 200);

  await field('currency').selectOption('CNY');
  for (const name of ['capex', 'monthlyOpex', 'terminalResidual', 'perMillionLocalInputCost', 'prices.peak.output']) {
    assert.equal(await field(name).inputValue(), '0', 'currency change clears money');
  }
  assert.equal(await select.inputValue(), 'custom');
  assert.equal(await select.locator('option').count(), 2, 'no invented CNY conversion of USD quotes');
  await field('currency').selectOption('USD');
  await select.selectOption(noCache.id);
  assert.equal(await field('apiCacheHit').inputValue(), '0');
  assert.ok(await field('apiCacheHit').getAttribute('readonly') !== null);
  await fill('capex', 10000);
  await fill('monthlyOpex', 100);
  await fill('inputTokensPerSecond', 1200);
  await fill('outputTokensPerSecond', 300);
  const usd = await exportScenario();
  assert.equal(usd.currency, 'USD');
  assert.equal(usd.metadata.pricing.id, noCache.id);
  await checkChartValues(usd);
  await checkTokenExamples(usd);
  for (const width of [320, 390, 768, 1440]) {
    await page.setViewportSize({width, height: 1000});
    for (const lang of ['en', 'zh', 'both']) {
      await page.locator(`[data-language-choice=${lang}]`).click();
      await page.waitForTimeout(60);
      assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), `${width}/${lang}: page overflow`);
      assert.ok(await charts.evaluate(el => el.scrollWidth <= el.clientWidth + 1), `${width}/${lang}: chart overflow`);
      for (const area of [examples, teaching]) {
        assert.ok(await area.evaluate(el => el.scrollWidth <= el.clientWidth + 1), `${width}/${lang}: token comparison overflow`);
      }
      assert.ok(await examples.locator('.te-matrix caption').evaluate(el =>
        el.getBoundingClientRect().width >= el.closest('table').getBoundingClientRect().width * .95
      ), `${width}/${lang}: matrix caption must retain the full table width`);
      assert.doesNotMatch(await charts.textContent(), /NaN|Infinity/);
      if (lang !== 'both' && [390, 1440].includes(width)) {
        await teaching.screenshot({path: `test-results/economics/token-teaching-${width}-${lang}.png`, style: screenshotStyle});
        await examples.screenshot({path: `test-results/economics/token-mixes-${width}-${lang}.png`, style: screenshotStyle});
        await charts.screenshot({path: `test-results/economics/charts-${width}-${lang}.png`, style: screenshotStyle});
        await charts.locator('[data-chart=cash]').screenshot({path: `test-results/economics/cash-${width}-${lang}.png`, style: screenshotStyle});
        await page.locator('#dr-calculator details').evaluateAll(nodes => nodes.forEach(n => n.open = false));
        await page.locator('#framework').evaluate(el => el.scrollIntoView({block: 'start'}));
        await page.screenshot({path: `test-results/economics/opening-${width}-${lang}.png`, style: screenshotStyle});
        await page.locator('.dr-worked').screenshot({path: `test-results/economics/worked-ledger-${width}-${lang}.png`, style: screenshotStyle});
        await page.locator('.dr-metadata').screenshot({path: `test-results/economics/model-filter-${width}-${lang}.png`, style: screenshotStyle});
        await page.locator('#dr-calc-output').screenshot({path: `test-results/economics/results-${width}-${lang}.png`, style: screenshotStyle});
      }
    }
  }
  await page.setViewportSize({width: 1440, height: 1000});
  await page.locator('[data-language-choice=zh]').click();
  await page.locator('#dr-calculator details').evaluateAll(nodes => nodes.forEach(n => n.open = false));
  await page.locator('#dr-calculator').screenshot({path: 'test-results/economics/calculator-desktop.png', style: screenshotStyle});
  await page.locator('#dr-calculator button[type=reset]').click();
  assert.equal(await select.inputValue(), '');
  assert.equal(await charts.locator('svg').count(), 0);
  assert.equal(await examples.locator('[data-token-examples]').count(), 0);
  assert.equal(await field('capex').inputValue(), '0');

  // Failure uses a local mock; no provider API, publication or paid service is called.
  const offline = await context.newPage();
  await offline.route('**/model-prices.json', route => route.fulfill({status: 503, body: ''}));
  await offline.goto(url);
  await offline.waitForFunction(() => document.querySelector('[data-catalog-status]').textContent.includes('unavailable'));
  await offline.locator('#dr-presetId').selectOption('custom');
  await offline.locator('#dr-modelName').fill('Offline manual model');
  assert.equal(await offline.locator('[data-action=download]').isDisabled(), false);
  await offline.close();

  const nojs = await browser.newPage({javaScriptEnabled: false});
  await nojs.goto(url);
  assert.ok(await nojs.locator('#framework').isVisible());
  assert.ok(await nojs.locator('#dr-generic-example').isVisible());
  assert.equal(await nojs.locator('#dr-token-examples-teaching [data-bill]').count(), 6);
  assert.deepEqual(await nojs.locator('#dr-token-examples-teaching [data-bill]').evaluateAll(nodes => nodes.map(n => Number(n.dataset.bill))), [260, 152, 500, 440, 740, 728]);
  assert.ok(await nojs.locator('noscript').first().isVisible());
  await nojs.close();
  assert.deepEqual(errors, []);
  console.log('Economics: CC Switch upstream imports, provider/search filters, monthly-bill ordering, million-token/currency units, provenance, reset/invalid states, model/currency changes, dynamic plots, token-mix/cache/volume comparisons, 12 responsive language layouts, export and offline/no-JS reading passed.');
  await context.close();
} finally {
  await browser.close();
  server.close();
}
