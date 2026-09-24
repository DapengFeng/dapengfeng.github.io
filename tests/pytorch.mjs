import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {execFileSync} from 'node:child_process';
import {load} from 'cheerio';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';

const source=load(await fs.readFile('content/posts/pytorch-01-what-is-pytorch.html','utf8'));
const temp=await fs.mkdtemp(path.join(os.tmpdir(),'feng-pytorch-check-'));
try{
 const file=path.join(temp,'model.cpp'),bin=path.join(temp,'model');
 await fs.writeFile(file,source('#pt-cpp-model').text());
 execFileSync(process.env.CXX||'g++',['-std=c++17','-O2','-Wall','-Wextra',file,'-o',bin],{timeout:30000});
 assert.equal(execFileSync(bin,{encoding:'utf8',timeout:10000}).trim(),'58 64\n139 154');
}finally{await fs.rm(temp,{recursive:true,force:true});}
assert.equal(source('code[data-godbolt="c++"]').length,1);
assert.equal(source('code[data-godbolt="python"]').length,0,'PyTorch requires its own local environment');
assert.equal(source('[data-lang] .shared-equation,[data-lang] .pt-demo').length,0,'share math and experiments');
const server=serve(4199),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'}),errors=[];
 page.on('pageerror',error=>errors.push(error.message));
 await page.addInitScript(()=>localStorage.setItem('feng-language','en'));
 await page.goto('http://localhost:4199/blog/pytorch-01-what-is-pytorch.html');
 assert.equal(await page.locator('.pt-demo').count(),2);
 assert.equal(await page.locator('.compiler-check').count(),1);
 assert.equal(await page.locator('[data-pt-stage][aria-current="step"]').getAttribute('data-pt-stage'),'7');
 await page.locator('#pt-route > .reading-explore > summary').click();
 await page.locator('#pt-multiply > .reading-explore > summary').click();
 await page.locator('[data-pt-stage="6"]').click();
 assert.match(await page.locator('#pt-route-info').innerText(),/mm_out_cpu/);
 await page.locator('[data-pt-device="cuda"]').click();
 assert.match(await page.locator('#pt-route-info').innerText(),/mm_out_cuda/);
 await page.locator('[data-pt-grad="off"]').click();await page.locator('[data-pt-stage="7"]').click();
 assert.match(await page.locator('#pt-route-info').innerText(),/grad_fn is None/);
 await page.locator('[data-language-choice="zh"]').click();
 assert.match(await page.locator('#pt-route-info').innerText(),/不含梯度历史/);
 assert.equal(await page.locator('#pt-route-info [data-lang="en"]').first().isVisible(),false);
 assert.equal(await page.locator('#pt-route').getAttribute('data-device'),'cuda');
 await page.locator('[data-language-choice="both"]').click();
 for(let i=0;i<12;i++)await page.locator('#pt-matrix-next').click();
 assert.deepEqual(await page.locator('#pt-matrix-y .pt-cell').allTextContents(),['58','64','139','154']);
 await page.locator('#pt-a00').fill('-2');
 assert.deepEqual(await page.locator('#pt-matrix-y .pt-cell').allTextContents(),['37','40','139','154']);
 for(let i=0;i<12;i++)await page.locator('#pt-matrix-next').click();
 assert.deepEqual(await page.locator('#pt-matrix-y .pt-cell').allTextContents(),['37','40','139','154']);
 // Backward graph uses values independently checked against PyTorch 2.10.0 CPU.
 await page.locator('[data-pt-grad="on"]').click();
 await page.locator('[data-pt-phase="backward"]').click();
 assert.equal(await page.locator('#pt-forward-scroll').isVisible(),false);
 assert.equal(await page.locator('#pt-backward-scroll').isVisible(),true);
 assert.equal(await page.locator('#pt-loss-value').textContent(),'415');
 await page.locator('#pt-reduce').click();
 assert.deepEqual(await page.locator('#pt-g-values .pt-mini-matrix>span').allTextContents(),['1','1','1','1']);
 assert.equal(await page.locator('#pt-ag-values').textContent(),'None');
 await page.locator('#pt-da').click();
 assert.deepEqual(await page.locator('#pt-da-values .pt-mini-matrix>span').allTextContents(),['15','19','23','15','19','23']);
 assert.equal(await page.locator('#pt-ag-values').textContent(),'None','contribution not accumulated yet');
 await page.locator('#pt-acc-a').click();
 assert.deepEqual(JSON.parse(await page.locator('#pt-route').getAttribute('data-grad-a')),[15,19,23,15,19,23]);
 assert.deepEqual(JSON.parse(await page.locator('#pt-route').getAttribute('data-grad-b')),[5,5,7,7,9,9]);
 assert.equal(await page.locator('#pt-route-next').isDisabled(),true);
 await page.locator('[data-pt-loss="mean"]').click();
 assert.equal(await page.locator('#pt-loss-value').textContent(),'103.75');
 assert.equal(await page.locator('#pt-route').getAttribute('data-grad-a'),'null','comparison resets leaf gradients');
 await page.locator('#pt-reduce').click();
 assert.equal(await page.locator('#pt-reduce b').textContent(),'MeanBackward0');
 assert.deepEqual(await page.locator('#pt-g-values .pt-mini-matrix>span').allTextContents(),['0.25','0.25','0.25','0.25']);
 await page.locator('#pt-acc-b').click();
 assert.deepEqual(JSON.parse(await page.locator('#pt-route').getAttribute('data-grad-a')),[3.75,4.75,5.75,3.75,4.75,5.75]);
 assert.deepEqual(JSON.parse(await page.locator('#pt-route').getAttribute('data-grad-b')),[1.25,1.25,1.75,1.75,2.25,2.25]);
 await page.waitForFunction(()=>document.querySelectorAll('#pt-backward-graph .pt-edge').length===8);
 assert.equal(await page.locator('#pt-backward-graph .pt-edge[data-saved=true]').count(),1);
 await page.locator('[data-pt-grad="off"]').click();
 assert.equal(await page.locator('#pt-route-next').isDisabled(),true);
 assert.equal(await page.locator('#pt-route-play').isDisabled(),true);
 assert.equal(await page.locator('#pt-acc-a').isDisabled(),true);
 assert.equal(await page.locator('#pt-route').getAttribute('data-grad-a'),'null');
 assert.match(await page.locator('#pt-route-info').innerText(),/No backward path was recorded/);
 await page.locator('[data-pt-grad="on"]').click();
 // Finite forward-to-backward playback is covered in reading-demos.mjs.
 for(const width of [1440,768,390,320]){
  await page.setViewportSize({width,height:1000});
  for(const language of ['en','zh','both']){
  await page.locator(`[data-language-choice="${language}"]`).click();
  assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${width}/${language}`);
  assert.ok(await page.locator('.site-skip').evaluate(node=>node.getBoundingClientRect().bottom<0),'skip link remains hidden until focused');
  for(const direction of ['forward','backward']){
   await page.locator(`[data-pt-phase="${direction}"]`).click();
   await page.locator('#pt-route-reset').click();
   const count=direction==='forward'?7:5;
   for(let i=0;i<count;i++)await page.locator('#pt-route-next').click();
   await page.waitForFunction(direction=>{
    const region=document.querySelector(`#pt-${direction}-scroll`);
    const n=region.querySelector('[aria-current=step]').getBoundingClientRect(),r=region.getBoundingClientRect();
    return n.top>=r.top&&n.bottom<=r.bottom+1&&n.left>=r.left&&n.right<=r.right+1;
   },direction);
   assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
  }

  }
 }
 await page.locator('.site-skip').focus();
 assert.ok(await page.locator('.site-skip').evaluate(node=>node.getBoundingClientRect().top>=0),'keyboard focus reveals skip link');
 await page.locator('#pt-route-reset').focus();
 await page.setViewportSize({width:1440,height:1000});
 await page.locator('[data-language-choice="both"]').click();
 await page.locator('#pt-route').screenshot({path:'/tmp/feng-pytorch-route.png'});
 await page.setViewportSize({width:390,height:1000});
 await page.locator('#pt-multiply').screenshot({path:'/tmp/feng-pytorch-matrix-mobile.png'});
 if(process.argv.includes('--live')){
  await page.locator('.compiler-run').click();
  await page.waitForFunction(()=>['passed','failed','rejected'].includes(document.querySelector('.compiler-check').dataset.state),{},{timeout:45000});
  assert.equal(await page.locator('.compiler-check').getAttribute('data-state'),'passed',await page.locator('.compiler-check').innerText());
  assert.equal((await page.locator('.compiler-stdout').innerText()).trim(),'58 64\n139 154');
 }
 await page.goto('http://localhost:4199/series/pytorch-internals/');
 assert.ok(await page.locator('.series-episodes a').count()>=2);
 for(const slug of ['pytorch-01-what-is-pytorch','pytorch-02-tensor-strides-storage'])
  assert.equal(await page.locator(`.series-episodes a[href="/blog/${slug}.html"]`).count(),1);
 assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),'series mobile overflow');
 const context=await browser.newContext({javaScriptEnabled:false});
 const staticPage=await context.newPage();await staticPage.goto('http://localhost:4199/blog/pytorch-01-what-is-pytorch.html');
 assert.ok((await staticPage.locator('#article-content').innerText()).includes('mm_out_cpu'));
 assert.equal(await staticPage.locator('#pt-route-next').isDisabled(),true);
 assert.deepEqual(errors,[]);await context.close();
 console.log('PyTorch checks passed: backward gradients, reductions, graph edges, auto-scroll, no-grad blocking, C++ output, route modes, arithmetic, playback, language preservation, responsive layout, series, and no-JS reading.');
}finally{await browser.close();server.close();}
