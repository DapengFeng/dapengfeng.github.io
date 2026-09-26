import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {execFileSync} from 'node:child_process';
import {load} from 'cheerio';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const article='pytorch-04-autograd-engine.html';
const source=load(await fs.readFile('content/posts/'+article,'utf8'));
const temp=await fs.mkdtemp(path.join(os.tmpdir(),'feng-autograd-'));
try{
 const file=path.join(temp,'model.cpp'),bin=path.join(temp,'model');
 await fs.writeFile(file,source('#autograd-cpp-model').text());
 execFileSync(process.env.CXX||'g++',['-std=c++17','-O2','-Wall','-Wextra',file,'-o',bin],{timeout:30000});
 for(const x of [2,3,-1,0,2.5]){
  const output=execFileSync(bin,[String(x)],{encoding:'utf8',timeout:10000});
  assert.equal((output.match(/^run /gm)||[]).length,5,'each node runs once, including the shared node');
  assert.ok(output.indexOf('u buffer=2, pending=1')<output.indexOf('u buffer=5, pending=0'));
  assert.ok(output.indexOf('u buffer=5, pending=0')<output.indexOf('run u=x*x: 5'));
  const actual=Number(output.match(/x\.grad=([^\n]+)/)[1]);
  const f=z=>2*(z*z)+3*(z*z),epsilon=1e-5;
  const numerical=(f(x+epsilon)-f(x-epsilon))/(2*epsilon);
  assert.ok(Math.abs(actual-numerical)<1e-7,`finite difference at ${x}`);
 }
}finally{await fs.rm(temp,{recursive:true,force:true});}
const server=serve(4221),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));
 await page.addInitScript(()=>localStorage.setItem('feng-language','both'));
 await page.goto('http://localhost:4221/blog/'+article);
 const lab=page.locator('#autograd-lab');
 assert.equal(await lab.getAttribute('data-demo-frame'),'complete');
 assert.equal(await page.locator('.compiler-check').count(),1);
 assert.equal(await page.locator('.ag-edge').count(),6,'keep the two separate edges to the same leaf');
 assert.equal(await page.locator('#ag-grad').textContent(),'20');
 await page.locator('.ag-controls summary').click();
 for(const order of ['left','right']){
  await page.locator(`[data-ag-order=${order}]`).click();
  const states=[[2,0,'None'],[2,0,'None'],[1,order==='left'?2:3,'None'],[0,5,'None'],[0,5,'None'],[0,5,'20']];
  for(let stage=0;stage<states.length;stage++){
   await lab.locator('[data-icon=step]').click();
   assert.equal(await lab.getAttribute('data-stage'),String(stage));
   const [pending,received,grad]=states[stage];
   assert.equal(await lab.getAttribute('data-pending'),String(pending));
   assert.equal(await lab.getAttribute('data-received'),String(received));
   assert.equal(await page.locator('#ag-grad').textContent(),grad);
   if(stage===3)assert.equal(await page.locator('#ag-square').getAttribute('data-active'),'false','ready does not mean already executed');
   if(stage===4)assert.equal(await page.locator('#ag-square').getAttribute('data-active'),'true');
  }
  await lab.locator('[data-icon=complete]').click();
 }
 for(const width of [1440,768,390,320]){
  await page.setViewportSize({width,height:1000});
  for(const lang of ['en','zh','both']){
   await page.locator(`[data-language-choice=${lang}]`).click();
   assert.equal(await lab.getAttribute('data-order'),'right','language preserves scenario');
   assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${width}/${lang}`);
   if(lang!=='both')assert.equal(await page.locator(`#ag-caption [data-lang=${lang==='en'?'zh':'en'}]`).isVisible(),false);
   assert.ok(await page.locator('.ag-node').evaluateAll(nodes=>nodes.every(n=>n.scrollHeight<=n.clientHeight+1)),'node text fits');
  }
 }
 await page.setViewportSize({width:1440,height:1000});await page.locator('[data-language-choice=both]').click();
 await page.locator('.ag-controls summary').click();
 await page.evaluate(()=>document.activeElement?.blur());await page.mouse.move(0,0);
 await lab.screenshot({path:'/tmp/feng-autograd-1440.png'});
 await page.setViewportSize({width:390,height:1000});await page.locator('[data-language-choice=zh]').click();
 await page.evaluate(()=>document.activeElement?.blur());
 await lab.screenshot({path:'/tmp/feng-autograd-390.png'});
 await page.setViewportSize({width:1440,height:1000});
 await lab.locator('[data-icon=replay]').click();
 await page.waitForFunction(()=>document.querySelector('#autograd-lab').dataset.stage==='2');
 assert.equal(await lab.getAttribute('data-received'),'3');
 await page.waitForFunction(()=>document.querySelector('#autograd-lab').dataset.demoFrame==='complete');
 assert.equal(await lab.getAttribute('data-demo-playing'),'false');
 const nojs=await browser.newPage({javaScriptEnabled:false,viewport:{width:390,height:1000}});
 await nojs.goto('http://localhost:4221/blog/'+article);
 assert.equal(await nojs.locator('.ag-node:visible').count(),5);
 assert.equal(await nojs.locator('.ag-edge').count(),6);
 assert.equal(await nojs.locator('#ag-grad').textContent(),'20');
 assert.equal(await nojs.locator('.ag-controls').isVisible(),false);
 assert.ok(await nojs.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
 await nojs.close();assert.deepEqual(errors,[]);
 console.log('Autograd: finite-difference C++ validation, shared-node dependencies, duplicate leaf edges, both orders, playback, language, layout and no-JS diagram passed.');
}finally{await browser.close();server.close();}
