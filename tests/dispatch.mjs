import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {execFileSync} from 'node:child_process';
import {load} from 'cheerio';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const article='pytorch-03-operator-dispatch.html';
const source=load(await fs.readFile('content/posts/'+article,'utf8'));
const temp=await fs.mkdtemp(path.join(os.tmpdir(),'feng-dispatch-'));
try{
 const cpp=path.join(temp,'model.cpp'),bin=path.join(temp,'model');
 await fs.writeFile(cpp,source('#dispatch-cpp-model').text());
 execFileSync(process.env.CXX||'g++',['-std=c++17','-O2','-Wall','-Wextra',cpp,'-o',bin],{timeout:30000});
 assert.equal(execFileSync(bin,{encoding:'utf8',timeout:10000}).trim(),'Autograd: record=1\nCPU: twice\nAttach backward rule\nvalue=6\nAutograd: record=0\nCPU: twice\nvalue=6\nMissing CUDA kernel');
}finally{await fs.rm(temp,{recursive:true,force:true});}
const server=serve(4209),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));
 await page.addInitScript(()=>localStorage.setItem('feng-language','both'));
 await page.goto('http://localhost:4209/blog/'+article);
 const lab=page.locator('#dispatch-lab');
 assert.equal(await lab.getAttribute('data-demo-frame'),'complete','reduced motion shows a complete diagram');
 assert.equal(await page.locator('.compiler-check').count(),1,'only standalone C++ uses Godbolt');
 await page.locator('.dk-compare summary').click();
 // Outcomes independently reproduced by the article's PyTorch 2.10 CPU experiment.
 const cases=[['grad','CPU',true,4],['plain','CPU',false,4],['no-grad','CPU',false,4],['inference','CPU',false,2],['cuda','CUDA',true,4]];
 for(const [mode,backend,record,count]of cases){
  await page.locator(`[data-dk-mode="${mode}"]`).click();
  assert.equal(await lab.getAttribute('data-records-grad'),String(record));
  assert.equal(await page.locator('.dk-card:visible').count(),count);
  assert.match(await page.locator('#dk-backend').textContent(),new RegExp('mm_out_'+backend.toLowerCase()));
  assert.match(await page.locator('#dk-result').innerText(),record?/MmBackward0/:/None/);
  const keys=await page.locator('#dk-keys code').allTextContents();
  assert.deepEqual(keys,mode==='inference'?[backend]:['Autograd'+backend,backend]);
  for(const width of [1440,768,390,320]){
   await page.setViewportSize({width,height:1000});
   for(const lang of ['en','zh','both']){
    await page.locator(`[data-language-choice=${lang}]`).click();
    assert.equal(await lab.getAttribute('data-mode'),mode,'language does not reset comparison');
    assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${mode}/${lang}/${width}`);
    if(lang!=='both')assert.equal(await page.locator(`#dk-result [data-lang=${lang==='en'?'zh':'en'}]`).isVisible(),false);
   }
  }
 }
 await page.setViewportSize({width:1440,height:1000});await page.locator('[data-language-choice=both]').click();
 await page.locator('[data-dk-mode=grad]').click();await page.locator('.dk-compare summary').click();
 for(const width of [1440,390]){await page.setViewportSize({width,height:1000});await lab.screenshot({path:`/tmp/feng-dispatch-${width}.png`});}
 // Explicit replay is allowed under reduced motion; it terminates after one pass.
 await page.setViewportSize({width:1440,height:1000});
 await lab.locator('[data-icon=replay]').click();
 await page.waitForFunction(()=>document.querySelector('#dispatch-lab').dataset.demoFrame==='1');
 assert.equal(await lab.locator('[data-active=true]').getAttribute('data-dk-stage'),'1');
 await page.waitForFunction(()=>document.querySelector('#dispatch-lab').dataset.demoFrame==='complete');
 assert.equal(await lab.getAttribute('data-demo-playing'),'false');
 await page.locator('.dk-compare summary').click();await page.locator('[data-dk-mode=inference]').click();
 await lab.locator('[data-icon=replay]').click();
 await page.waitForFunction(()=>document.querySelector('#dispatch-lab').dataset.demoFrame==='1');
 assert.equal(await lab.locator('[data-active=true]').getAttribute('data-dk-stage'),'3','inference bypasses the wrapper');
 await lab.locator('[data-icon=complete]').click();
 const nojs=await browser.newPage({javaScriptEnabled:false});
 await nojs.goto('http://localhost:4209/blog/'+article);
 assert.equal(await nojs.locator('.dk-card:visible').count(),4);
 assert.match(await nojs.locator('#dk-result').innerText(),/11/);
 assert.equal(await nojs.locator('.dk-compare').isVisible(),false);
 await nojs.close();
 assert.deepEqual(errors,[]);
 console.log('Dispatch: C++ wrapper/redispatch/error output, five routes, inference bypass, finite replay, languages, responsive layout and no-JS diagram passed.');
}finally{await browser.close();server.close();}
