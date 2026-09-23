import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const server=serve(4195),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage();await page.addInitScript(()=>localStorage.setItem('feng-language','en'));
 await page.goto('http://localhost:4195/blog/');await page.locator('#article-sort').selectOption('title');
 for(const mode of ['en','zh','both','en']){
  await page.locator(`[data-language-choice=${mode}]`).click();
  const lang=mode==='zh'?'zh':'en';
  const titles=await page.locator(`[data-library] .knowledge-card h3 [data-lang=${lang}]`).allTextContents();
  const collator=new Intl.Collator(lang==='zh'?'zh-CN':'en',{numeric:true,sensitivity:'base'});
  assert.deepEqual(titles,[...titles].sort(collator.compare),`${mode}: alphabetical order must match the visible language`);
 }
 await page.locator('[data-category-filter=systems]').click();
 await page.locator('[data-language-choice=zh]').click();
 assert.equal(await page.locator('[data-category-filter=systems]').getAttribute('aria-pressed'),'true');
 await page.locator('#article-sort').selectOption('oldest');await page.locator('[data-language-choice=en]').click();
 const dates=await page.locator('[data-library] .knowledge-card:visible').evaluateAll(nodes=>nodes.map(n=>n.dataset.date));
 assert.deepEqual(dates,[...dates].sort(),'language changes preserve date sorting');
 await page.goto('http://localhost:4195/categories/');
 for(const mode of ['en','zh','both']){
  await page.locator(`[data-language-choice=${mode}]`).click();
  const headings=await page.locator('.category-group h2').allInnerTexts();
  if(mode==='en')assert.ok(headings.every(t=>!/[\u4e00-\u9fff]/.test(t)));
  if(mode==='zh')assert.ok(headings.every(t=>!/[A-Za-z]/.test(t)));
  if(mode==='both')assert.ok(headings.every(t=>/[A-Za-z]/.test(t)&&/[\u4e00-\u9fff]/.test(t)));
 }
 await page.locator('.search-trigger').click();await page.locator('#global-search').fill('cuTile');
 await page.waitForFunction(()=>document.querySelectorAll('.search-result').length===1);
 for(const mode of ['en','zh','both']){
  await page.evaluate(mode=>document.querySelector(`[data-language-choice=${mode}]`).click(),mode);
  const text=await page.locator('.search-result small').innerText();
  assert.equal(text.includes('SYSTEMS'),mode!=='zh');assert.equal(text.includes('系统与性能'),mode!=='en');
 }
 console.log('Reading regressions passed: language-aware sorting, filter persistence, category headings, and search metadata.');
}finally{await browser.close();server.close();}
