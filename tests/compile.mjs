import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {execFileSync} from 'node:child_process';
import {load} from 'cheerio';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const article='pytorch-06-compile-and-codegen.html';
const source=load(await fs.readFile('content/posts/'+article,'utf8'));
const temp=await fs.mkdtemp(path.join(os.tmpdir(),'feng-compile-'));
try{
 const file=path.join(temp,'fusion.cpp'),bin=path.join(temp,'fusion');
 await fs.writeFile(file,source('#compile-fusion-cpp').text());
 execFileSync(process.env.CXX||'g++',['-std=c++17','-O2','-Wall','-Wextra',file,'-o',bin],{timeout:30000});
 for(const size of [1,4,7,32]){
  const result=execFileSync(bin,[String(size)],{encoding:'utf8',timeout:10000});
  const values=result.split('\n')[0].slice('output:'.length).trim().split(' ').map(Number);
  const inputs=[-2,-1,1,2];
  assert.deepEqual(values,Array.from({length:Math.min(size,8)},(_,i)=>2*Math.max(inputs[i%4]+1,0)));
  const [,separate,fused]=result.match(/logical accesses: (\d+) -> (\d+)/);
  assert.equal(Number(separate),7*size);assert.equal(Number(fused),3*size);
 }
 for(const node of source('code.language-python[id]').toArray()){
  execFileSync('python3',['-c','import ast,sys;ast.parse(sys.stdin.read())'],{input:source(node).text(),timeout:10000});
  assert.equal(source(node).attr('data-godbolt'),undefined);
 }
}finally{await fs.rm(temp,{recursive:true,force:true});}
const server=serve(4225),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));
 await page.addInitScript(()=>localStorage.setItem('feng-language','both'));
 await page.goto('http://localhost:4225/blog/'+article);
 const lab=page.locator('#compile-route-lab');
 await page.waitForFunction(()=>document.querySelector('#compile-route-lab')?.dataset.demoFrame==='complete');
 assert.equal(await page.locator('.compiler-check').count(),1);
 assert.equal(await lab.locator('.cp-edge').count(),7);
 await page.locator('.cp-controls summary').click();
 for(const [name,hit,variants]of [['first',false,[4]],['reuse',true,[4]],['resize',false,[4,8]],['return',true,[4,8]]]){
  await page.locator(`[data-cp-case=${name}]`).click();
  assert.deepEqual(JSON.parse(await lab.getAttribute('data-variants')),variants);
  const steps=hit?['guards','cache','run','output']:['guards','dynamo','aot','inductor','cache','run','output'];
  for(const [i,node]of steps.entries()){
   await lab.locator('[data-icon=step]').focus();await page.keyboard.press('Enter');
   assert.equal(await lab.getAttribute('data-stage'),String(i));
   assert.equal(await lab.locator('.cp-node[data-active=true]').getAttribute('id'),'cp-'+node);
   assert.equal(await lab.getAttribute('data-compile-visited'),String(!hit&&i>=1));
   if(i===0&&!hit)assert.deepEqual(JSON.parse(await lab.getAttribute('data-variants')),name==='first'?[]:[4]);
  }
  await lab.locator('[data-icon=complete]').click();
  assert.equal(await page.locator('#cp-result').textContent(),name==='resize'?'[8]':'[4]');
  for(const width of [1440,768,390,320]){
   await page.setViewportSize({width,height:1000});
   for(const language of ['en','zh','both']){
    await page.locator(`[data-language-choice=${language}]`).click();
    assert.equal(await lab.getAttribute('data-case'),name);
    assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${name}/${width}/${language}`);
    assert.ok(await lab.locator('.cp-node').evaluateAll(nodes=>nodes.every(n=>n.scrollHeight<=n.clientHeight+1&&n.scrollWidth<=n.clientWidth+1)),`node fit: ${name}/${width}/${language}`);
    if(language!=='both')assert.equal(await page.locator(`#cp-caption [data-lang=${language==='en'?'zh':'en'}]`).isVisible(),false);
   }
  }
 }
 await page.setViewportSize({width:1440,height:1000});await page.locator('[data-language-choice=both]').click();
 await page.locator('[data-cp-case=first]').click();await page.locator('.cp-controls summary').click();
 await page.evaluate(()=>document.activeElement?.blur());await page.mouse.move(0,0);
 await lab.screenshot({path:'/tmp/feng-compile-1440.png'});
 await page.locator('#compile-fusion').screenshot({path:'/tmp/feng-compile-fusion-1440.png'});
 await page.setViewportSize({width:390,height:1000});await page.locator('[data-language-choice=zh]').click();
 await page.evaluate(()=>document.activeElement?.blur());await lab.screenshot({path:'/tmp/feng-compile-390.png'});
 await page.emulateMedia({reducedMotion:'no-preference'});
 await lab.locator('[data-icon=replay]').click();
 await page.waitForFunction(()=>document.querySelector('#compile-route-lab').dataset.stage==='1');
 assert.equal(await lab.locator('.cp-dot').isVisible(),true);
 await page.evaluate(()=>scrollTo(0,0));await page.waitForTimeout(250);
 const state=await lab.locator('.cp-dot').getAttribute('cx');
 await page.waitForTimeout(250);assert.equal(await lab.locator('.cp-dot').getAttribute('cx'),state);
 await lab.scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelector('#compile-route-lab').dataset.demoFrame==='complete');
 assert.equal(await lab.getAttribute('data-demo-playing'),'false');
 assert.equal(await lab.locator('.cp-dot').isVisible(),false);
 const nojs=await browser.newPage({javaScriptEnabled:false,viewport:{width:390,height:1000}});
 await nojs.goto('http://localhost:4225/blog/'+article);
 assert.equal(await nojs.locator('.cp-node:visible').count(),7);
 assert.equal(await nojs.locator('.cp-edge').count(),7);
 assert.equal(await nojs.locator('.cp-controls').isVisible(),false);
 assert.equal(await nojs.locator('.cp-pane').count(),2);
 assert.ok(await nojs.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
 await nojs.close();assert.deepEqual(errors,[]);
 console.log('Compile: C++ fusion arithmetic/access counts, Python syntax, four cache scenarios, compiler bypass, finite animation, languages, mobile layout and no-JS diagram passed.');
}finally{await browser.close();server.close();}
