import {chromium, webkit} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';

// Small shared-UI sweep. Larger scientific demonstrations have dedicated suites.
const server = serve(4225);
const origin = 'http://localhost:4225';
const directory = 'test-results/visual';
await fs.mkdir(directory, {recursive: true});
try {
 for (const [engine, browserType] of [['chromium', chromium], ['webkit', webkit]]) {
  const browser = await browserType.launch({headless: true, ...(engine === 'chromium' ? {executablePath: chromiumExecutable(), args: ['--no-sandbox']} : {})});
  try {
   const context = await browser.newContext({reducedMotion: 'reduce', locale: 'zh-CN'});
   await context.addInitScript(() => {
    localStorage.setItem('feng-language', 'zh');
   });
   await context.route('https://giscus.app/**', route => route.fulfill({contentType: 'text/html', body: '<!doctype html><html lang="zh-CN"><label>评论<textarea></textarea></label></html>'}));
   // Deterministic fallback exercises a real renderer import failure without GPU dependencies.
   await context.route('**/assets/eye-renderer.js', route => route.abort());
   const page = await context.newPage(), errors = [], requests = [];
   const snapshotTime = Date.parse('2026-10-06T08:00:00Z');
   await page.clock.install({time: new Date(snapshotTime)});
   await page.clock.pauseAt(new Date(snapshotTime + 3600000));
   page.on('pageerror', error => errors.push(error.message));
   page.on('request', request => requests.push(request.url()));
   for (const width of [1440, 390]) {
    await page.setViewportSize({width, height: 900});
    for (const [name, path, body] of [
     ['home', '/', '#featured'],
     ['technical', '/blog/human-visual-system.html', '#vision-eye-3d'],
     ['journal', '/blog/chaoshan-streets-and-sea.html', '#shantou-arcades']
    ]) {
     await page.goto(origin + path);
     await page.evaluate(() => document.fonts.ready);
     await page.clock.runFor(32);
     if (name === 'home') {
      await page.waitForSelector('[data-rendered="true"]');
      assert.equal(await page.locator('.featured-card').count(), 3);
     }
     assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1), `${engine}: ${path} fits ${width}`);
     await page.screenshot({path: `${directory}/${engine}-${name}-${width}-top.png`});
     await page.locator(body).scrollIntoViewIfNeeded();
     await page.clock.runFor(32);
     if (name === 'technical') {
      await page.waitForSelector('[data-eye-state="fallback"]');
      assert.ok(await page.locator('.eye-poster').isVisible(), 'anatomy remains available without WebGL');
     }
     if (name === 'journal') assert.equal(await page.locator('script[src$="reading-demos.js"],script[src$="syntax.js"]').count(), 0);
     await page.screenshot({path: `${directory}/${engine}-${name}-${width}-body.png`});
    }
   }
   await page.clock.resume();
   assert.ok(!requests.some(url => url.includes('api.country.is')), 'no IP lookup');
   await page.goto(origin + '/blog/chaoshan-streets-and-sea.html');
   await page.locator('.mobile-menu').click();
   assert.equal(await page.locator('.mobile-menu').getAttribute('aria-expanded'), 'true');
   await page.keyboard.press('Escape');
   assert.equal(await page.locator('.mobile-menu').getAttribute('aria-expanded'), 'false');
   await page.locator('[data-share-open]').click();
   assert.ok(await page.locator('#article-share-dialog').isVisible());
   await page.locator('[data-share-platform=wechat]').click();
   await page.locator('.share-qr img').waitFor({state: 'visible'});
   await page.keyboard.press('Escape');
   assert.ok(await page.locator('[data-share-open]').evaluate(node => node === document.activeElement));
   if (await page.locator('[data-support-open]').count()) {
    await page.locator('[data-support-open]').click();
    await page.locator('[data-support-custom]').fill('12.50');
    assert.match(await page.locator('.support-paypal-link').getAttribute('href'), /\/12\.5USD$/);
    await page.keyboard.press('Escape');
   }
   await page.locator('.discussion-launcher').click();
   await page.frameLocator('.discussion-frame').locator('textarea').fill('Local test draft');
   await page.locator('.discussion-close').click();
   assert.deepEqual(errors, [], `${engine} shared UI errors`);
   await context.close();
   console.log(`${engine}: layout, fallback, navigation, native dialog, sharing and discussion passed; fixed screenshots saved.`);
  } finally { await browser.close(); }
 }
} finally { server.close(); }
