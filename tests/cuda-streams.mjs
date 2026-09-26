import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {execFileSync} from 'node:child_process';
import {load} from 'cheerio';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const article='pytorch-05-cuda-streams-timing.html';
const source=load(await fs.readFile('content/posts/'+article,'utf8'));
const temp=await fs.mkdtemp(path.join(os.tmpdir(),'feng-cuda-streams-'));
try{
 const file=path.join(temp,'model.cpp'),bin=path.join(temp,'model');
 await fs.writeFile(file,source('#cuda-stream-cpp').text());
 execFileSync(process.env.CXX||'g++',['-std=c++17','-O2','-Wall','-Wextra',file,'-o',bin],{timeout:30000});
 // Independent constraints: ordered consumers cannot start before either
 // submission or producer completion; unordered consumers have no such edge.
 for(const [a,b] of [[6,3],[1,3],[10,1],[3,9]]){
  const output=execFileSync(bin,[String(a),String(b)],{encoding:'utf8',timeout:10000});
  const lines=output.trim().split('\n');assert.equal(lines.length,3);
  for(const [i,line] of lines.entries()){
   const [,done,start,end,ordered]=line.match(/A_done=(\d+), B=\[(\d+),(\d+)\], ordered=(true|false)/);
   assert.equal(Number(done),1+a);assert.equal(Number(end)-Number(start),b);
   if(i<2){assert.ok(Number(start)>=Number(done));assert.ok(Number(start)>=2);assert.equal(Number(start),Math.max(2,1+a));}
   else assert.equal(Number(start),2);
   assert.equal(ordered,String(Number(start)>=Number(done)));
  }
 }
 for(const node of source('code.language-python').toArray()){
  const text=source(node).text();
  execFileSync('python3',['-c','import ast,sys; ast.parse(sys.stdin.read())'],{input:text,timeout:10000});
  assert.equal(source(node).attr('data-godbolt'),undefined,'GPU programs are not submitted to a CPU-only compiler');
 }
}finally{await fs.rm(temp,{recursive:true,force:true});}
const server=serve(4223),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));
 await page.addInitScript(()=>localStorage.setItem('feng-language','both'));
 await page.goto('http://localhost:4223/blog/'+article);
 const lab=page.locator('#cuda-streams-lab');
 await page.waitForFunction(()=>document.querySelector('#cuda-streams-lab')?.dataset.demoFrame==='complete');
 assert.equal(await page.locator('.compiler-check').count(),1);
 await page.locator('.cs-controls summary').click();
 for(const mode of ['same','event','unsafe']){
  await page.locator(`[data-cs-mode=${mode}]`).click();
  const expectedStart=mode==='unsafe'?2:7;
  assert.equal(await lab.getAttribute('data-b-start'),String(expectedStart));
  assert.equal(await lab.getAttribute('data-b-end'),String(expectedStart+3));
  assert.equal(await lab.getAttribute('data-ordered'),String(mode!=='unsafe'));
  assert.equal(await lab.locator('.cs-wait').isVisible(),mode==='event');
  assert.equal(await page.locator(`${mode==='same'?'#cs-stream0':'#cs-stream1'} .cs-b`).count(),1);
  if(mode==='unsafe')assert.match(await page.locator('#cs-caption').textContent(),/No dependency/);
  // Keyboard stepping must expose host return while device work remains pending.
  await lab.locator('[data-icon=step]').focus();await page.keyboard.press('Enter');
  assert.equal(await lab.getAttribute('data-time'),'1');
  await page.keyboard.press('Enter');
  assert.equal(await lab.getAttribute('data-time'),'2');
  assert.equal(await page.locator('#cs-host-status').textContent(),'✓ · 2');
  assert.equal(await page.locator('#cs-a-status').textContent(),'… · 7');
  assert.equal(await page.locator('.cs-b').getAttribute('data-state'),mode==='unsafe'?'running':'pending');
  await page.keyboard.press('Enter');
  assert.equal(await lab.getAttribute('data-time'),'7');
  assert.equal(await page.locator('.cs-b').getAttribute('data-state'),mode==='unsafe'?'done':'running');
  await lab.locator('[data-icon=complete]').click();
  for(const width of [1440,768,390,320]){
   await page.setViewportSize({width,height:1000});
   for(const lang of ['en','zh','both']){
    await page.locator(`[data-language-choice=${lang}]`).click();
    assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${mode}/${width}/${lang}`);
    assert.equal(await lab.getAttribute('data-mode'),mode);
    if(lang!=='both')assert.equal(await page.locator(`#cs-caption [data-lang=${lang==='en'?'zh':'en'}]`).isVisible(),false);
    assert.ok(await lab.locator('.cs-status>div').evaluateAll(nodes=>nodes.every(n=>n.scrollWidth<=n.clientWidth+1)));
   }
  }
 }
 await page.setViewportSize({width:1440,height:1000});await page.locator('[data-language-choice=both]').click();
 await page.locator('[data-cs-mode=event]').click();await page.locator('.cs-controls summary').click();
 await page.evaluate(()=>document.activeElement?.blur());await page.mouse.move(0,0);
 await lab.screenshot({path:'/tmp/feng-cuda-streams-1440.png'});
 await page.setViewportSize({width:390,height:1000});await page.locator('[data-language-choice=zh]').click();
 await page.evaluate(()=>document.activeElement?.blur());
 await lab.screenshot({path:'/tmp/feng-cuda-streams-390.png'});
 // Actual motion, finite completion and visibility pause, independent of reduced motion.
 await page.emulateMedia({reducedMotion:'no-preference'});
 await lab.locator('[data-icon=replay]').click();
 await page.waitForFunction(()=>+document.querySelector('#cuda-streams-lab').dataset.time>0.1);
 await page.evaluate(()=>scrollTo(0,0));await page.waitForTimeout(250);
 const pausedTime=Number(await lab.getAttribute('data-time'));
 await page.waitForTimeout(250);assert.equal(Number(await lab.getAttribute('data-time')),pausedTime);
 await lab.scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelector('#cuda-streams-lab').dataset.demoFrame==='complete');
 assert.equal(await lab.getAttribute('data-demo-playing'),'false');
 assert.equal(await lab.getAttribute('data-time'),'10');
 const nojs=await browser.newPage({javaScriptEnabled:false,viewport:{width:390,height:1000}});
 await nojs.goto('http://localhost:4223/blog/'+article);
 assert.equal(await nojs.locator('#cs-stream1 .cs-b').count(),1);
 assert.equal(await nojs.locator('.cs-wait').isVisible(),true);
 assert.equal(await nojs.locator('#cs-b-status').textContent(),'✓ · 10');
 assert.equal(await nojs.locator('.cs-controls').isVisible(),false);
 assert.ok(await nojs.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
 await nojs.close();assert.deepEqual(errors,[]);
 console.log('CUDA streams: C++ dependency bounds, Python syntax, three timeline scenarios, host/device boundaries, playback, visibility pause, languages, responsive layout and no-JS diagram passed.');
}finally{await browser.close();server.close();}
