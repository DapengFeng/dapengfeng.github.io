import {chromiumExecutable} from './browser-options.mjs';
import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import {serve} from '../scripts/serve.mjs';
const server=serve(4194);
const browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));
 await page.addInitScript(()=>{
  localStorage.setItem('feng-language','both');
  Object.defineProperty(navigator,'clipboard',{value:{writeText:async text=>{if(window.failCopy)throw Error('Denied');window.copiedMath=text;}}});
 });
 for(const slug of ['band-storage-gaxpy','spike_notes']){
  await page.goto(`http://localhost:4194/blog/${slug}.html`);
  await page.evaluate(()=>document.querySelectorAll('details').forEach(n=>n.open=true));
  const formula=page.locator('.formula-block').first(), button=formula.locator('.formula-copy');
  const latex=await formula.getAttribute('data-latex');
  assert.equal(await page.locator('.formula-inline .formula-copy').count(),0);
  for(const mode of ['en','zh','both']){
   await page.locator(`[data-language-choice=${mode}]`).click();
   await button.focus();await page.keyboard.press('Enter');
   await page.waitForFunction(()=>window.copiedMath!==undefined);
   assert.equal(await page.evaluate(()=>window.copiedMath),latex);
   assert.equal(await button.getAttribute('title'),mode==='en'?'Copied':mode==='zh'?'已复制':'Copied / 已复制');
   assert.equal(await button.locator('svg').count(),1);
   assert.ok(await button.locator('.copied-symbol').isVisible());
   assert.equal(await button.getAttribute('data-copied'),'true');
  }
  await page.evaluate(()=>window.failCopy=true);await button.click();
  assert.equal(await formula.locator('textarea').inputValue(),latex);
  assert.equal(await formula.locator('textarea').evaluate(n=>n.selectionEnd-n.selectionStart),latex.length);
  await page.evaluate(()=>window.failCopy=false);await button.click();
  assert.equal(await formula.locator('textarea').count(),0);
  for(const width of [1440,768,390,320]){
   await page.setViewportSize({width,height:1000});
   await page.evaluate(()=>new Promise(requestAnimationFrame));
   assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${slug}/${width}: page overflow`);
   assert.ok(await page.locator('.formula-scroll').evaluateAll(nodes=>nodes.every(n=>getComputedStyle(n).overflowX==='auto')));
   if([1440,390].includes(width)){
    await formula.scrollIntoViewIfNeeded();
    await page.screenshot({path:`/tmp/feng-math-${slug}-${width}.png`});
   }
  }
 }
 const offline=await browser.newPage({javaScriptEnabled:false});
 await offline.goto('http://localhost:4194/blog/spike_notes.html');
 assert.ok(await offline.locator('.formula-block svg').count()>60);
 assert.ok(await offline.locator('.formula-inline svg').count()>200);
 assert.deepEqual(errors,[]);
 console.log('AMS SVG, bilingual copy feedback, keyboard copying, clipboard fallback, mobile overflow, and no-JavaScript rendering passed.');
}finally{await browser.close();server.close();}
