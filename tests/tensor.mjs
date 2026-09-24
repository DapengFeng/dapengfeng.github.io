import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {execFileSync} from 'node:child_process';
import {load} from 'cheerio';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const article='pytorch-02-tensor-strides-storage.html';
const source=load(await fs.readFile('content/posts/'+article,'utf8'));
const temp=await fs.mkdtemp(path.join(os.tmpdir(),'feng-tensor-'));
try{
 await fs.writeFile(path.join(temp,'model.cpp'),source('#tensor-cpp-model').text());
 execFileSync(process.env.CXX||'g++',['-std=c++17','-O2','-Wall','-Wextra',path.join(temp,'model.cpp'),'-o',path.join(temp,'model')],{timeout:30000});
 assert.equal(execFileSync(path.join(temp,'model'),{encoding:'utf8'}).trim(),'6 7\n60 60');
}finally{await fs.rm(temp,{recursive:true,force:true});}
assert.equal(source('[data-lang] .tl-demo,[data-lang] .shared-equation').length,0);
// Expected address traces independently verified with PyTorch 2.10.0+cpu.
const cases={base:[0,1,2,3,4,5,6,7,8,9,10,11],transpose:[0,4,8,1,5,9,2,6,10,3,7,11],slice:[1,3,5,7,9,11],copy:[0,1,2,3,4,5,6,7,8,9,10,11],expand:[0,1,2,3,0,1,2,3,0,1,2,3]};
const server=serve(4207),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));await page.addInitScript(()=>localStorage.setItem('feng-language','both'));
 await page.goto('http://localhost:4207/blog/'+article);
 assert.equal(await page.locator('#tensor-layout').getAttribute('data-layout'),'transpose');
 await page.locator('#tensor-layout > .reading-explore > summary').click();
 for(const [mode,expected] of Object.entries(cases)){
  await page.locator(`button[data-layout=${mode}]`).click();
  assert.deepEqual((await page.locator('#tl-trace').textContent()).split(' → ').map(Number),expected);
  for(let i=0;i<expected.length;i++){
   await page.locator(`#tl-matrix [data-index="${i}"]`).click();
   assert.equal(Number(await page.locator('#tensor-layout').getAttribute('data-offset')),expected[i]);
   assert.equal(await page.locator('#tl-storage [data-active=true]').getAttribute('data-slot'),String(expected[i]));
   const value=mode==='copy'?cases.transpose[i]:expected[i];
   assert.equal(await page.locator('#tl-matrix [aria-pressed=true]').textContent(),String(value));
  }
 }
 assert.equal(await page.locator('#tl-matrix [data-alias=true]').count(),2);
 await page.locator('button[data-layout=transpose]').click();
 await page.locator('#tl-matrix [data-index="7"]').click();
 assert.match(await page.locator('#tl-equation').textContent(),/\[2, 1\].*= 6.*24 B/);
 await page.locator('#tl-next').click();assert.equal(await page.locator('#tensor-layout').getAttribute('data-offset'),'10');
 await page.locator('#tl-reset').click();await page.locator('#tl-scan').click();
 await page.waitForFunction(()=>Number(document.querySelector('#tensor-layout').dataset.selected)>=2);
 await page.locator('#tl-scan').click();assert.equal(await page.locator('#tl-scan').getAttribute('aria-pressed'),'false');
 await page.locator('button[data-layout=copy]').click();await page.locator('#tl-matrix [data-index="1"]').click();
 assert.equal(await page.locator('#tl-original [data-active=true]').getAttribute('data-slot'),'4');
 for(const width of [1440,768,390,320]){
  await page.setViewportSize({width,height:1000});
  for(const lang of ['en','zh','both']){
   await page.locator(`[data-language-choice=${lang}]`).click();
   for(const mode of Object.keys(cases)){
    await page.locator(`button[data-layout=${mode}]`).click();
    assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${mode}/${width}/${lang}`);
    if(lang!=='both')assert.equal(await page.locator(`#tl-readout [data-lang=${lang==='en'?'zh':'en'}]`).isVisible(),false);
   }
  }
 }
 for(const width of [1440,390]){
  await page.setViewportSize({width,height:1000});await page.locator('[data-language-choice=both]').click();
  await page.locator('button[data-layout=transpose]').click();await page.locator('#tl-matrix [data-index="7"]').click();
  await page.locator('#tensor-layout').screenshot({path:`/tmp/feng-tensor-${width}.png`});
 }
 const nojs=await browser.newPage({javaScriptEnabled:false});await nojs.goto('http://localhost:4207/blog/'+article);
 assert.equal(await nojs.locator('.tl-fallback').isVisible(),true);assert.equal(await nojs.locator('.tl-interactive').isVisible(),false);await nojs.close();
 assert.deepEqual(errors,[]);console.log('Tensor: C++ output, five address maps, aliasing, scan controls, three languages, four widths and no-JS fallback passed.');
}finally{await browser.close();server.close();}
