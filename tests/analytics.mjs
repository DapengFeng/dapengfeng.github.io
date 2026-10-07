import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {chromium} from '@playwright/test';
import {chromiumExecutable} from './browser-options.mjs';
import {shell} from '../scripts/templates.mjs';
import {site} from '../scripts/config.mjs';
import {localizePage} from '../scripts/i18n.mjs';
import {optimizePage} from '../scripts/seo.mjs';

const token = '0123456789abcdef0123456789abcdef';
const variable = 'CLOUDFLARE_WEB_ANALYTICS_TOKEN';
const previous = process.env[variable];
let enabled, disabled;
try {
  delete process.env[variable];
  disabled = shell({title:'Analytics', body:'<main id="main">Article</main>'});
  assert.ok(!disabled.includes('/assets/analytics.js'));
  process.env[variable] = '  ';
  assert.ok(!shell({title:'Empty', body:''}).includes('/assets/analytics.js'));
  for (const invalid of [variable, '`token`', '<script>alert(1)</script>', '" onload="alert(1)']) {
    process.env[variable] = invalid;
    assert.throws(() => shell({title:'Invalid', body:''}), /32-character token/);
  }
  process.env[variable] = ` ${token}\n`;
  enabled = optimizePage(localizePage(shell({title:'Analytics', body:'<main id="main">Article</main>'}), []), '/');
} finally {
  if (previous === undefined) delete process.env[variable];
  else process.env[variable] = previous;
}

const loader = await fs.readFile('dist/assets/analytics.js', 'utf8');
const browser = await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try {
  for (const {name, url, html=enabled, privacy={}, expected=0} of [
    {name:'production', url:site.url+'/', expected:1},
    {name:'missing token', url:site.url+'/', html:disabled},
    {name:'localhost preview', url:'http://localhost:4173/'},
    {name:'HTTPS preview', url:'https://preview.example/'},
    {name:'HTTP', url:site.url.replace('https:', 'http:')+'/'},
    {name:'Do Not Track', url:site.url+'/', privacy:{doNotTrack:'1'}},
    {name:'Global Privacy Control', url:site.url+'/', privacy:{globalPrivacyControl:true}},
  ]) {
    const context = await browser.newContext();
    await context.addInitScript(values => {
      for (const [key, value] of Object.entries(values)) Object.defineProperty(navigator, key, {get:()=>value});
    }, privacy);
    const page = await context.newPage(), errors=[];
    let requests=0;
    page.on('pageerror', error=>errors.push(error.message));
    // No traffic leaves the browser: exercise the real module with a mock beacon.
    await context.route('**/*', async route => {
      const request = route.request(), address = new URL(request.url());
      if (request.isNavigationRequest()) return route.fulfill({contentType:'text/html',body:html});
      if (address.pathname === '/assets/analytics.js') return route.fulfill({contentType:'text/javascript',body:loader});
      if (address.hostname === 'static.cloudflareinsights.com') {
        requests++;
        return route.fulfill({contentType:'text/javascript',headers:{'Access-Control-Allow-Origin':'*'},body:'window.beaconLoaded = true;'});
      }
      return route.fulfill({contentType:request.resourceType()==='stylesheet'?'text/css':'text/javascript',body:''});
    });
    await page.goto(url);
    await page.waitForLoadState('networkidle');
    assert.equal(requests, expected, name);
    assert.deepEqual(errors, [], name);
    const beacon = page.locator('script[data-cf-beacon]');
    assert.equal(await beacon.count(), expected, name);
    if (expected) {
      assert.deepEqual(JSON.parse(await beacon.getAttribute('data-cf-beacon')), {token});
      assert.equal(await beacon.getAttribute('type'), 'module');
      assert.equal(await page.evaluate(()=>window.beaconLoaded), true);
    }
    await context.close();
  }
  console.log('Analytics checks passed: optional build configuration, token validation, production-only loading, and privacy opt-outs. No live analytics requests sent.');
} finally {
  await browser.close();
}
