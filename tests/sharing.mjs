import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {createRequire} from 'node:module';
const require=createRequire(import.meta.url);
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const server=serve(4219),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
const canonical='https://dapengfeng.github.io/blog/jiuzhaigou-water-and-mountains.html';
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[],assets=[];
 page.on('pageerror',e=>errors.push(e.message));
 page.on('request',r=>{if(r.url().includes('/assets/share/'))assets.push(r.url());});
 await page.route('https://giscus.app/**',r=>r.fulfill({contentType:'text/html',body:'<!doctype html><title>Mock discussion</title>'}));
 await page.addInitScript(()=>{
  if(window!==window.top)return;
  localStorage.setItem('feng-language','en');window.copies=[];window.shares=[];
  Object.defineProperty(navigator,'clipboard',{configurable:true,value:{writeText:async text=>{if(window.denyCopy)throw Error('denied');window.copies.push(text);}}});
  Object.defineProperty(navigator,'share',{configurable:true,value:data=>{window.shares.push({...data,files:data.files?.map(f=>({name:f.name,type:f.type,size:f.size}))});return window.shareError?Promise.reject(new DOMException('mock',window.shareError)):Promise.resolve();}});
  Object.defineProperty(navigator,'canShare',{configurable:true,value:data=>!!data.files?.length});
 });
 await page.goto('http://localhost:4219/blog/jiuzhaigou-water-and-mountains.html?tracking=ignored#reflections');
 assert.equal(assets.length,0,'sharing assets must not delay initial reading');
 await page.addScriptTag({path:require.resolve('axe-core/axe.min.js')});
 const trigger=page.locator('[data-share-open]').first(),dialog=page.locator('#article-share-dialog');
 await trigger.click();assert.ok(await dialog.isVisible());
 await page.locator('[data-share-action=copy-link]').click();
 assert.equal(await page.evaluate(()=>copies.at(-1)),canonical);
 assert.equal(new URL(await page.locator('[data-share-platform=x]').getAttribute('href')).searchParams.get('url'),canonical);
 await page.locator('[data-share-native]').click();
 assert.equal(await page.evaluate(()=>shares.at(-1).url),canonical);
 await page.evaluate(()=>window.shareError='AbortError');await page.locator('[data-share-native]').click();
 await page.waitForFunction(()=>!document.querySelector('[data-share-native]').disabled);
 assert.equal(await page.locator('[data-share-status]').innerText(),'');
 await page.evaluate(()=>window.shareError='NotAllowedError');await page.locator('[data-share-native]').click();
 await page.waitForFunction(()=>document.querySelector('[data-share-status]').textContent.includes('unavailable'));
 await page.evaluate(()=>window.shareError='');
 await page.locator('#article-share-text').fill('My edited introduction');
 await page.locator('[data-share-action=copy-text]').click();assert.equal(await page.evaluate(()=>copies.at(-1)),'My edited introduction');
 await page.evaluate(()=>window.denyCopy=true);await page.locator('[data-share-action=copy-link]').click();
 assert.equal(await page.locator('#share-manual-text').inputValue(),canonical);
 assert.ok(await page.locator('.share-manual').isVisible());
 await page.evaluate(()=>window.denyCopy=false);
 await page.keyboard.press('Escape');assert.ok(!(await dialog.isVisible()));assert.ok(await trigger.evaluate(n=>n===document.activeElement));
 for(const mode of ['en','zh','both']){
  await page.locator(`[data-language-choice=${mode}]`).click();await trigger.click();
  const content=await page.locator('#article-share-text').inputValue();
  if(mode==='en')assert.equal(content,'My edited introduction');
  if(mode==='zh'){assert.ok(content.includes('九寨沟'));assert.ok(!content.includes('Jiuzhaigou'));}
  if(mode==='both')assert.ok(content.indexOf('Jiuzhaigou')<content.indexOf('九寨沟'));
  await page.locator('[data-share-platform=wechat]').click();await page.locator('[data-share-qr] img').evaluate(n=>n.decode());
  assert.ok(await page.locator('[data-share-qr]').isVisible());
  await page.locator('[data-share-platform=rednote]').click();await page.locator('[data-share-poster] img').evaluate(n=>n.decode());
  assert.ok(!(await page.locator('[data-share-qr]').isVisible()));
  assert.ok((await page.locator('[data-share-download]').getAttribute('href')).endsWith(`-${mode}.png`));
  await page.locator('[data-share-file]').waitFor({state:'visible'});await page.locator('[data-share-file]').click();
  assert.ok((await page.evaluate(()=>shares.at(-1).files[0].name)).endsWith(`-${mode}.png`));
  assert.equal(await page.evaluate(()=>shares.at(-1).files[0].type),'image/png');
  const violations=await page.evaluate(async()=>(await axe.run(document.getElementById('article-share-dialog'),{runOnly:{type:'tag',values:['wcag2a','wcag2aa','wcag21aa']}})).violations.map(v=>({id:v.id,nodes:v.nodes.map(n=>n.target)})));
  assert.deepEqual(violations,[],`sharing accessibility: ${mode}`);
  await page.keyboard.press('Shift+Tab');assert.ok(await dialog.evaluate(n=>n.contains(document.activeElement)));
  for(const width of [1440,768,390,320]){
   await page.setViewportSize({width,height:800});
   const rect=await dialog.boundingBox();assert.ok(rect.x>=0&&rect.x+rect.width<=width+1&&rect.y>=0&&rect.y+rect.height<=800);
   assert.ok(await dialog.evaluate(n=>n.scrollWidth<=n.clientWidth+1),'dialog has no horizontal overflow');
   assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
  }
  await page.keyboard.press('Escape');
 }
 await page.setViewportSize({width:1440,height:1000});await trigger.click();await page.locator('[data-share-platform=rednote]').click();
 await page.screenshot({path:'/tmp/feng-sharing-desktop.png'});
 await page.setViewportSize({width:390,height:844});await page.locator('[data-share-platform=wechat]').click();
 await page.screenshot({path:'/tmp/feng-sharing-mobile.png'});
 const download=await Promise.all([page.waitForEvent('download'),page.locator('[data-share-download]').click()]).then(r=>r[0]);
 const file=await download.path();assert.ok((await fs.stat(file)).size>1000);
 await page.keyboard.press('Escape');
 await page.evaluate(()=>{Object.defineProperty(navigator,'share',{configurable:true,value:undefined});Object.defineProperty(navigator,'clipboard',{configurable:true,value:undefined});});
 // A fresh document without native sharing retains the other sharing choices.
 await page.addInitScript(()=>{Object.defineProperty(navigator,'share',{configurable:true,value:undefined});});
 await page.reload();await trigger.click();assert.ok(!(await page.locator('[data-share-native]').isVisible()));
 await page.locator('[data-share-action=copy-link]').click();assert.equal(await page.evaluate(()=>copies.at(-1)),canonical);
 await page.keyboard.press('Escape');
 const nojs=await browser.newPage({javaScriptEnabled:false});await nojs.goto('http://localhost:4219/blog/jiuzhaigou-water-and-mountains.html');
 assert.equal(await nojs.locator('[data-share-open]:visible').count(),0);assert.ok(await nojs.locator('#article-content').isVisible());await nojs.close();
 assert.deepEqual(errors,[]);console.log('Sharing: canonical links, clipboard fallback, editable drafts, native/file sharing, cancellation, downloads, locales, keyboard and responsive checks passed (no public posts).');
}finally{await browser.close();server.close();}
