import {chromiumExecutable} from './browser-options.mjs';
import { chromium } from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {serve} from '../scripts/serve.mjs';
const server=serve(4175);
const executablePath=chromiumExecutable();
const browser=await chromium.launch({...(executablePath?{executablePath}:{}),headless:true,args:['--no-sandbox']});
const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[];page.on('pageerror',e=>errors.push(e.message));
try{
 await page.goto('http://localhost:4175/');assert.equal(await page.locator('.featured-card').count(),3);await page.locator('[data-language-choice="en"]').click();assert.equal(await page.locator('.hero-copy h1 [data-lang="zh"]').isVisible(),false);await page.reload();assert.equal(await page.locator('[data-language-choice="en"]').getAttribute('aria-pressed'),'true');assert.equal(await page.locator('select#site-language').count(),0);
 await page.locator('#surface-frequency').fill('2.4');await page.locator('#surface-frequency').dispatchEvent('input');assert.equal(await page.locator('#frequency-value').textContent(),'2.4');await page.locator('#surface-toggle').click();assert.equal(await page.locator('#surface-toggle').getAttribute('aria-pressed'),'true');
 await page.locator('.search-trigger').click();await page.locator('#global-search').fill('EventProp');await page.waitForSelector('.search-result');assert.ok((await page.locator('.search-results').innerText()).includes('spikes'));await page.keyboard.press('Escape');assert.equal(await page.locator('dialog').isVisible(),false);
 await page.goto('http://localhost:4175/blog/?category=physics');assert.equal(await page.locator('[data-library] .knowledge-card:visible').count(),1);await page.locator('[data-category-filter="all"]').click();await page.locator('#library-search').fill('CUDA');assert.equal(await page.locator('[data-library] .knowledge-card:visible').count(),1);await page.locator('#library-search').fill('no-such-result-1234');assert.equal(await page.locator('.empty-state').isVisible(),true);await page.locator('[data-reset-filters]').click();assert.equal(await page.locator('[data-library] .knowledge-card:visible').count(),10);
 await page.goto('http://localhost:4175/blog/spike_notes.html');assert.equal(await page.locator('.paired-article').count(),1);assert.equal(await page.locator('.article-toc nav a').count(),8);assert.equal(await page.locator('.parallel-text').count(),173);
 const first=page.locator('.chapter-lead.parallel-text').first();assert.equal(await first.locator('[data-lang="en"]').isVisible(),true);assert.equal(await first.locator('[data-lang="zh"]').isVisible(),false);
 await page.locator('[data-language-choice="both"]').click();assert.equal(await first.locator('[data-lang="zh"]').isVisible(),true);const englishBox=await first.locator('.parallel-en').boundingBox(),chineseBox=await first.locator('.parallel-zh').boundingBox();assert.ok(chineseBox.y>=englishBox.y+englishBox.height);
 const before=await page.locator('#lif-detail').innerText();await page.locator('#lif-input').fill('0.5');await page.locator('#lif-input').dispatchEvent('input');assert.notEqual(await page.locator('#lif-detail').innerText(),before);assert.ok((await page.locator('#lif-detail').innerText()).includes('no threshold crossing'));
 await page.locator('[data-language-choice="zh"]').click();assert.equal(await first.locator('[data-lang="en"]').isVisible(),false);assert.equal(await first.locator('[data-lang="zh"]').isVisible(),true);assert.ok(await page.evaluate(()=>Boolean(window.SpikeNotes)));
 await page.goto('http://localhost:4175/blog/rust-vs-cpp-blog.html');await page.locator('#step-next').click();await page.locator('[data-code-mode="fixed"]').first().click();assert.equal(errors.length,0,errors.join('\n'));
 await page.goto('http://localhost:4175/blog/benchmark-with-evidence.html');await page.locator('[data-language-choice="en"]').click();const lab=page.locator('[data-benchmark-lab]');await lab.locator('[data-benchmark-run]').click();await page.waitForFunction(()=>document.querySelector('[data-benchmark-lab] [data-benchmark-result]').textContent.includes('Verified'));assert.equal(await lab.locator('[data-benchmark-export]').isEnabled(),true);
 const downloadPromise=page.waitForEvent('download');await lab.locator('[data-benchmark-export]').click();const download=await downloadPromise;const data=JSON.parse(await fs.readFile(await download.path(),'utf8'));assert.equal(data.samples.plain.length,15);assert.equal(data.samples.unrolled.length,15);assert.ok(data.expected>0);
 // A single experiment survives all language changes without duplicated controls.
 const numbers=text=>text.replace('×4','').match(/\d+(?:\.\d+)?/g);
 const englishResult=await lab.locator('[data-benchmark-result]').innerText();
 await page.locator('[data-language-choice="zh"]').click();
 const chineseLab=page.locator('[data-benchmark-lab]');
 assert.equal(await chineseLab.count(),1);
 assert.equal(await chineseLab.locator('[data-benchmark-run]').count(),1);
 assert.deepEqual(numbers(await chineseLab.locator('[data-benchmark-result]').innerText()),numbers(englishResult));
 assert.equal(await chineseLab.locator('[data-benchmark-export]').isEnabled(),true);
 await page.locator('[data-language-choice="both"]').click();
 assert.equal(await page.locator('[data-benchmark-lab]:visible').count(),1);
 assert.equal(await chineseLab.locator('.benchmark-row').count(),6);
 assert.match(await chineseLab.locator('[data-benchmark-result]').innerText(),/Verified/);
 assert.match(await chineseLab.locator('[data-benchmark-result]').innerText(),/校验通过/);
 await chineseLab.locator('[data-benchmark-size]').selectOption('65536');
 assert.equal(await chineseLab.locator('[data-benchmark-export]').isEnabled(),false);
 assert.match(await chineseLab.locator('[data-benchmark-result]').innerText(),/Not run yet/);
 await page.goto('http://localhost:4175/blog/waves-and-phase.html');
 await page.locator('[data-language-choice="en"]').click();
 await page.locator('[data-wave-phase]').fill('2.5');
 await page.locator('[data-wave-k]').fill('2.1');
 await page.locator('[data-wave-pause]').click();
 const paused=await page.locator('[data-wave-pause]').getAttribute('aria-pressed');
 await page.locator('[data-language-choice="zh"]').click();
 assert.equal(await page.locator('[data-wave-phase]').inputValue(),'2.5');
 assert.equal(await page.locator('[data-wave-pause]').getAttribute('aria-pressed'),paused);
 assert.equal(await page.locator('[data-wave-k]').inputValue(),'2.1');
 await page.locator('[data-language-choice="both"]').click();
 assert.equal(await page.locator('[data-wave-lab]:visible').count(),1);
 assert.equal(await page.locator('[data-wave-lab] canvas').count(),1);
 assert.equal(await page.locator('[data-wave-phase]').inputValue(),'2.5');
 assert.equal(await page.locator('[data-wave-pause]').getAttribute('aria-pressed'),paused);
 assert.equal(await page.locator('[data-wave-pause] [data-lang="en"]').isVisible(),true);
 assert.equal(await page.locator('[data-wave-pause] [data-lang="zh"]').isVisible(),true);
 await page.goto('http://localhost:4175/blog/rust-vs-cpp-blog.html');
 await page.locator('[data-language-choice="en"]').click();
 for(const key of ['unreal','qt','cuda','embedded','legacy']){
  await page.locator(`[data-scenario="${key}"]`).click();
  await page.waitForFunction(()=>document.querySelector('#scenario-result [data-lang="en"]'));
  assert.doesNotMatch(await page.locator('#scenario-result').innerText(),/[\u3400-\u9fff]/);
 }
 await page.locator('[data-code-mode="fixed"]').click();
 await page.waitForFunction(()=>document.querySelector('#rust-result [data-lang="en"]'));
 assert.match(await page.locator('#rust-result').innerText(),/independent integer/);
 await page.setViewportSize({width:390,height:844});
 for(const url of ['/','/blog/','/categories/','/archive/','/about/','/blog/spike_notes.html','/blog/rust-vs-cpp-blog.html','/blog/cuda-rust-two-tracks-blog.html','/blog/waves-and-phase.html','/blog/benchmark-with-evidence.html']){
  await page.goto('http://localhost:4175'+url);for(const mode of ['both','en','zh']){await page.locator(`[data-language-choice="${mode}"]`).click();await page.evaluate(()=>new Promise(requestAnimationFrame));assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${url} overflow in ${mode}`);}
 }
 assert.equal(errors.length,0,errors.join('\n'));console.log('Browser checks passed: bilingual persistence, search, filters, experiments, export, mobile layout, no script errors.');
}finally{await browser.close();server.close();}
