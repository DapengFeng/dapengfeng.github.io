import {chromiumExecutable} from './browser-options.mjs';
import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {serve} from '../scripts/serve.mjs';
const server=serve(4197);
const browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1280,height:1000}}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));
 await page.addInitScript(()=>{
  localStorage.setItem('feng-language','both');
  Object.defineProperty(navigator,'clipboard',{value:{writeText:async text=>{if(window.failCopy)throw Error('Denied');window.copiedCode=text;}}});
 });
 for(const file of (await fs.readdir('content/posts')).filter(f=>f.endsWith('.html'))){
  await page.goto(`http://localhost:4197/blog/${file}`);
  assert.equal(await page.locator('#article-content pre:not(.compiler-highlight)').evaluateAll(nodes=>nodes.filter(n=>!n.hidden&&!n.closest('[aria-hidden="true"]')&&!n.parentElement.matches('.code-copy-block')).length),0,file);
  assert.equal(await page.locator('.compiler-editor-surface .code-copy').count(),0,'no duplicate overlay controls');
  assert.equal(await page.locator('#article-content button[data-copy]:visible').count(),0);
  assert.equal(await page.locator('.code-copy-block').count(),await page.locator('.code-language').count(),file+' language labels');
 }
 await page.goto('http://localhost:4197/blog/pytorch-01-what-is-pytorch.html');
 const source=page.locator('#pt-python-main');
 const block=source.locator('xpath=ancestor::div[contains(@class,"code-copy-block")][1]'),button=block.locator('.code-copy');
 const expected=await source.textContent();
 assert.equal(await block.locator('.code-language').textContent(),'Python');
 assert.equal(await page.locator('.compiler-editor-surface').locator('..').locator('.code-language').textContent(),'C++');
 assert.ok(await page.locator('.code-language[data-code-language=yaml]').count()>0);
 for(const mode of ['en','zh','both']){
  await page.locator(`[data-language-choice=${mode}]`).click();
  await button.focus();await page.keyboard.press('Enter');
  await page.waitForFunction(()=>window.copiedCode!==undefined);
  assert.equal(await page.evaluate(()=>window.copiedCode),expected);
  assert.equal(await button.getAttribute('title'),mode==='en'?'Copied':mode==='zh'?'已复制':'Copied / 已复制');
  assert.equal(await button.locator('svg').count(),1);
 }
 await page.evaluate(()=>window.failCopy=true);await button.click();
 assert.equal(await block.locator('textarea').inputValue(),expected);
 assert.equal(await block.locator('textarea').evaluate(n=>n.selectionEnd-n.selectionStart),expected.length);
 await page.evaluate(()=>window.failCopy=false);await button.click();assert.equal(await block.locator('textarea').count(),0);
 const editor=page.locator('.compiler-editor').first(),editable=editor.locator('xpath=ancestor::div[contains(@class,"code-copy-block")][1]');
 const edited='#include <iostream>\nint main() {\n  std::cout << "<edited>";\n}\n';
 await editor.fill(edited);await editable.locator('.code-copy').click();assert.equal(await page.evaluate(()=>window.copiedCode),edited);
 await page.evaluate(()=>{
  const pre=document.createElement('pre');pre.id='dynamic-copy-test';pre.innerHTML='<code><span class="code-line"><span class="line-number">1</span>int x = 1;</span><span class="code-line"><span class="line-number">2</span>  x += 2;</span></code>';
  document.getElementById('article-content').append(pre);
 });
 const dynamic=page.locator('#dynamic-copy-test').locator('..');await dynamic.locator('.code-copy').click();assert.equal(await page.evaluate(()=>window.copiedCode),'int x = 1;\n  x += 2;');
 assert.equal(await dynamic.locator('.code-language').getAttribute('data-code-language'),'text');
 await page.locator('#dynamic-copy-test').evaluate(n=>n.textContent='new output\n  <raw>');
 await dynamic.locator('.code-copy').click();assert.equal(await page.evaluate(()=>window.copiedCode),'new output\n  <raw>');
 for(const width of [1440,768,390,320]){
  await page.setViewportSize({width,height:1000});
  assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
 }
 await page.setViewportSize({width:390,height:1000});await block.scrollIntoViewIfNeeded();await page.screenshot({path:'/tmp/feng-code-copy-mobile.png'});
 assert.deepEqual(errors,[]);
 console.log('All article code blocks covered; clipboard, edited source, dynamic output, line numbers, languages, keyboard, fallback and mobile checks passed.');
}finally{await browser.close();server.close();}
