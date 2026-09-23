import {chromiumExecutable} from './browser-options.mjs';
import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {serve} from '../scripts/serve.mjs';
const live=process.argv.includes('--live'),server=serve(4181);
const executablePath=chromiumExecutable();
const browser=await chromium.launch({...(executablePath?{executablePath}:{}),headless:true,args:['--no-sandbox','--no-proxy-server']});
try {
 const page=await browser.newPage({viewport:{width:1280,height:960},reducedMotion:'reduce'}),requests=[],errors=[];
 page.on('pageerror',e=>errors.push(e.message));
 let offline=false,held=null,release=null,hold=false,customResult=null;
 if(!live)await page.route('https://godbolt.org/api/compiler/*/compile',async route=>{
  const body=route.request().postDataJSON();requests.push(body);
  if(offline){await route.abort('failed');return;}
  if(hold){held=true;await new Promise(resolve=>release=resolve);}
  if(customResult){await route.fulfill({contentType:'application/json',body:JSON.stringify(customResult)});return;}
  const rustFailure=body.lang==='rust'&&body.source.includes('let first = &items[0]');
  const result=body.lang==='python'?{code:0,buildResult:{code:0},didExecute:true,stdout:[{text:'3.782015427334855 3.782015427067353'}],stderr:[]}:{code:rustFailure?1:0,stderr:rustFailure?[{text:'error[E0502]: cannot borrow `items` as mutable'}]:[],stdout:[],asm:rustFailure?[]:[{text:'main:'},{text:'  ret'}],execResult:rustFailure?undefined:{code:body.source.includes('return 3')?3:0,didExecute:true,buildResult:{code:0},stdout:[{text:body.source.includes('edited output')?'edited output':'10'}],stderr:body.source.includes('return 3')?[{text:'runtime note'}]:[]}};
  await route.fulfill({contentType:'application/json',body:JSON.stringify(result)});
 });
 await page.goto('http://localhost:4181/blog/rust-vs-cpp-blog.html');
 await page.locator('[data-language-choice="en"]').click();
 const rust=page.locator('.rust-pane .compiler-check'),cpp=page.locator('.cpp-pane .compiler-check');
 assert.equal(await page.locator('.compiler-check').count(),2);
 assert.equal(requests.length,0,'do not send examples before a click');
 const verifyHighlight=async panel=>{
  assert.equal(await panel.locator('.compiler-highlight code').textContent(),await panel.locator('.compiler-editor').inputValue()+'\n');
  assert.ok(await panel.locator('.syntax-keyword').count()>0);
  assert.ok(await panel.locator('.syntax-string').count()>0);
 };
 await verifyHighlight(rust);await verifyHighlight(cpp);
 const check=async(panel,state)=>{await panel.locator('.compiler-run').click();await panel.locator('.compiler-status').waitFor();await page.waitForFunction(()=>![...document.querySelectorAll('.compiler-check')].some(n=>n.dataset.state==='pending'),{},{timeout:35000});assert.equal(await panel.getAttribute('data-state'),state,await panel.innerText());};
 await check(rust,'rejected');assert.match(await rust.locator('.compiler-output').innerText(),/E0502/);
 const originalDiagnostics=await rust.locator('.compiler-output').innerText();
 for(const mode of ['zh','both','en']){
  await page.locator(`[data-language-choice="${mode}"]`).click();
  assert.equal(await rust.locator('.compiler-output').innerText(),originalDiagnostics,'Godbolt diagnostics remain verbatim in every language');
  assert.equal(await rust.locator('.compiler-output [data-lang]').count(),0);
 }
 await check(cpp,'passed');
 await page.locator('[data-code-mode="fixed"]').click();assert.equal(await rust.locator('.compiler-results').isVisible(),false);
 await check(rust,'passed');await check(cpp,'passed');
 assert.equal((await rust.locator('.compiler-stdout').innerText()).trim(),'10');
 assert.equal((await cpp.locator('.compiler-stdout').innerText()).trim(),'10');
 assert.equal(await cpp.locator('.compiler-assembly').isVisible(),true);
 if(!live){
  assert.ok(requests.at(-2).source.includes('let first = items[0]'));
  assert.ok(!requests.at(-2).source.includes('1fn main'));
  assert.equal(requests.at(-2).options.filters.execute,true);
 }
 const verifyFallback=async(selector,en,zh)=>{
  const output=cpp.locator(selector);
  for(const mode of ['en','zh','both']){
   await page.locator(`[data-language-choice="${mode}"]`).click();
   assert.equal(await output.locator('[data-lang=en]').first().isVisible(),mode!=='zh');
   assert.equal(await output.locator('[data-lang=zh]').first().isVisible(),mode!=='en');
   const visible=await output.innerText();
   assert.equal(visible.includes(en),mode!=='zh');assert.equal(visible.includes(zh),mode!=='en');
  }
 };
 await verifyFallback('.compiler-output','No compiler diagnostics','无编译诊断');
 await verifyFallback('.compiler-stderr','No standard error','无标准错误输出');
 await page.locator('[data-language-choice=en]').click();
  const edited='#include <iostream>\nint main(){std::cout << "edited output";std::cerr << "runtime note";return 3;}';
  await cpp.locator('.compiler-editor').fill(edited);
  assert.equal(await cpp.locator('.compiler-results').isVisible(),false);
  await verifyHighlight(cpp);
  await check(cpp,'rejected');if(!live)assert.equal(requests.at(-1).source,edited);
  assert.match(await cpp.locator('.compiler-status').innerText(),/exit code 3/);
  assert.equal(await cpp.locator('.compiler-stdout').innerText(),'edited output');
  assert.equal(await cpp.locator('.compiler-stderr').innerText(),'runtime note');
  for(const mode of ['zh','both','en']){
   await page.locator(`[data-language-choice="${mode}"]`).click();
   assert.equal(await cpp.locator('.compiler-stdout').innerText(),'edited output');
   assert.equal(await cpp.locator('.compiler-stderr').innerText(),'runtime note');
  }
  await page.locator('[data-language-choice="zh"]').click();assert.equal(await cpp.locator('.compiler-editor').inputValue(),edited);
  await page.evaluate(()=>Object.defineProperty(navigator,'clipboard',{configurable:true,value:{writeText:async text=>window.copiedCode=text}}));
  await cpp.locator('.code-copy-block').filter({has:page.locator('.compiler-editor')}).locator('.code-copy').click();assert.equal(await page.evaluate(()=>window.copiedCode),edited);
  await cpp.locator('.compiler-reset').click();assert.equal(await cpp.locator('.compiler-results').isVisible(),false);
  assert.ok((await cpp.locator('.compiler-editor').inputValue()).includes('const int first'));
  await page.locator('[data-language-choice="en"]').click();
 if(!live){
  offline=true;await check(rust,'error');assert.equal(await rust.locator('.compiler-run').isEnabled(),true);offline=false;
  hold=true;await rust.locator('.compiler-run').click();await page.waitForTimeout(30);assert.ok(held);
  await rust.locator('.compiler-editor').fill('fn main() { println!("new draft"); }');release();await page.waitForTimeout(30);hold=false;
  assert.equal(await rust.locator('.compiler-results').isVisible(),false,'discard old in-flight result on code changes');
  assert.equal(await rust.locator('.compiler-run').isEnabled(),true);
  await page.locator('[data-code-mode="bad"]').click();assert.ok((await rust.locator('.compiler-editor').inputValue()).includes('let first = &items[0]'));
 }
 // A resized textarea and its painted layer must reveal the same source area.
 const resizeSource=Array.from({length:32},(_,i)=>`int visible_line_${i} = ${i};`).join('\n');
 await cpp.locator('.compiler-editor').fill(resizeSource);
 for(const height of [260,700,1100,300]){
  await cpp.locator('.compiler-editor').evaluate((n,height)=>{n.style.height=`${height}px`;n.scrollTop=0;},height);
  await page.waitForFunction(()=>{
   const e=document.querySelector('.cpp-pane .compiler-editor'),p=e.parentElement.querySelector('.compiler-highlight');
   return p.clientHeight===e.clientHeight;
  });
  const dimensions=await cpp.locator('.compiler-editor').evaluate(e=>{
   const p=e.parentElement.querySelector('.compiler-highlight');
   e.scrollTop=e.scrollHeight;e.dispatchEvent(new Event('scroll'));
   const token=p.querySelector('code').lastElementChild,r=token.getBoundingClientRect(),box=p.getBoundingClientRect();
   return {height:e.getBoundingClientRect().height,scroll:e.scrollTop,paintScroll:p.scrollTop,lastVisible:r.top>=box.top&&r.bottom<=box.bottom,overflow:getComputedStyle(p).overflow};
  });
  assert.equal(dimensions.height,height,'manual expansion must not stop at 720px');
  assert.equal(dimensions.scroll,dimensions.paintScroll,'scroll positions must stay aligned');
  assert.equal(dimensions.lastVisible,true,'bottom source line must be visible after resizing and scrolling');
  assert.equal(dimensions.overflow,'hidden','only the textarea owns scrollbars');
 }
 await cpp.locator('.compiler-reset').click();
 await cpp.locator('.compiler-editor').evaluate(n=>n.style.removeProperty('height'));
 await page.setViewportSize({width:390,height:844});
 for(const mode of ['en','zh','both']){
  await page.locator(`[data-language-choice="${mode}"]`).click();
  assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
 }
 if(!live){
  customResult={code:0,stdout:[],stderr:[],execResult:{code:0,didExecute:true,stdout:[],stderr:[],truncated:true}};
  await check(cpp,'passed');
  await verifyFallback('.compiler-stdout','No standard output','无标准输出');
  await verifyFallback('.compiler-stderr .compiler-output-notice','Output truncated by Godbolt','Godbolt 已截断输出');
  customResult={code:0,stdout:[],stderr:[],execResult:{code:0,didExecute:true,stdout:[{text:'<b>English output only</b>'}],stderr:[{text:'English error only'}],truncated:true}};
  await check(cpp,'passed');
  for(const mode of ['en','zh','both']){
   await page.locator(`[data-language-choice="${mode}"]`).click();
   assert.equal(await cpp.locator('.compiler-stdout').innerText(),'<b>English output only</b>');
   assert.equal(await cpp.locator('.compiler-stdout b').count(),0,'output is text, never HTML');
   assert.match(await cpp.locator('.compiler-stderr').innerText(),/^English error only/);
  }
  customResult={};await check(cpp,'error');
  await verifyFallback('.compiler-output','Unexpected API response','接口返回格式异常');
  customResult=null;
 }
 await page.goto('http://localhost:4181/blog/spike_notes.html');
 const python=page.locator('.compiler-check');await page.locator('#event-code').evaluate(n=>{let p=n.parentElement;while(p){if(p.tagName==='DETAILS')p.open=true;p=p.parentElement;}});
 await check(python,'passed');
 await verifyHighlight(python);
 const numbers=(await python.locator('.compiler-stdout').innerText()).trim().split(/\s+/).map(Number);
 assert.equal(numbers.length,2);assert.ok(Math.abs(numbers[0]-numbers[1])<1e-6);assert.ok(Math.abs(numbers[0]-3.782015)<1e-5);
 assert.match(page.url(),/localhost:4181\/blog\/spike_notes.html$/,'stay in the article');
 assert.equal(errors.length,0,errors.join('\n'));
 console.log(`${live?'Live Godbolt':'Isolated API'} checks passed: Rust rejection/fix, C++ execution, editable source, stdout/stderr, Python output, source switching, mobile layout, no navigation.`);
}finally{await browser.close();server.close();}
