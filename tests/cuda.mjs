import {chromiumExecutable} from './browser-options.mjs';
import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {execFileSync} from 'node:child_process';
import {load} from 'cheerio';
import {serve} from '../scripts/serve.mjs';
const slug='cuda-rust-two-tracks-blog',live=process.argv.includes('--live');
const source=load(await fs.readFile(`content/posts/${slug}.html`,'utf8'));
const temporary=await fs.mkdtemp(path.join(os.tmpdir(),'feng-cuda-'));
try {
 const file=path.join(temporary,'verify.cpp'),binary=path.join(temporary,'verify');
 await fs.writeFile(file,source('code[data-godbolt="c++"]').text());
 execFileSync(process.env.CXX||'g++',['-std=c++17','-O2','-Wall','-Wextra',file,'-o',binary],{timeout:30000});
 const output=execFileSync(binary,{encoding:'utf8',timeout:10000});
 assert.match(output,/Verified cases: 27/);assert.match(output,/9 missing, 8 multiply-written/);
 console.log(output.trim());
}finally{await fs.rm(temporary,{recursive:true,force:true});}
const server=serve(4190),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox','--no-proxy-server']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));await page.addInitScript(()=>localStorage.setItem('feng-language','en'));
 await page.goto(`http://localhost:4190/blog/${slug}.html`);
 await page.locator('#cuda-map > .reading-explore > summary').click();
 assert.equal(await page.locator('.article-toc nav a').count(),8);
 assert.equal(await page.locator('#cuda-map').count(),1);
 assert.equal(await page.locator('pre.syntax-block').count(),2);
 assert.ok(await page.locator('pre.syntax-block .syntax-keyword').count()>5);
 assert.equal(await page.locator('#cuda-threads .cuda-masked').count(),7);
 await page.locator('#cuda-bug').click();assert.match(await page.locator('#cuda-result').innerText(),/Missing outputs: 9; multiply-written outputs: 8/);
 // Validate both diagrams for partial, complete and one-element groups.
 for(const n of [1,7,8,9,17,32,64])for(const exponent of[1,2,3,4]){
  await page.locator('#cuda-n').fill(String(n));await page.locator('#cuda-b').fill(String(exponent));
  for(const broken of[false,true]){
   if(await page.locator('#cuda-bug').getAttribute('aria-pressed')!==String(broken))await page.locator('#cuda-bug').click();
   const b=2**exponent,g=Math.ceil(n/b),missing=broken?Math.max(0,n-b):0,conflicts=broken&&g>1?Math.min(b,n):0;
   assert.match(await page.locator('#cuda-result').innerText(),new RegExp(`Missing outputs: ${missing}; multiply-written outputs: ${conflicts}`));
   assert.deepEqual(await page.locator('.cuda-interval').allTextContents(),Array.from({length:g},(_,k)=>`[${k*b}, ${Math.min((k+1)*b,n)})`));
  }
 }
 await page.locator('#cuda-n').fill('17');await page.locator('#cuda-b').fill('3');await page.locator('#cuda-bug').click();
 for(const width of[1440,768,390,320]){
  await page.setViewportSize({width,height:1000});
  for(const mode of['en','zh','both']){
   await page.locator(`[data-language-choice=${mode}]`).click();
   assert.equal(await page.locator('#cuda-n').inputValue(),'17');
   assert.equal(await page.locator('#cuda-map').count(),1);
   assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`overflow ${width}/${mode}`);
  }
 }
 await page.setViewportSize({width:1440,height:1000});await page.locator('#cuda-map').screenshot({path:'/tmp/feng-cuda-map-desktop.png'});
 await page.setViewportSize({width:390,height:1000});await page.locator('#cuda-map').screenshot({path:'/tmp/feng-cuda-map-mobile.png'});
 if(live){
  await page.locator('.compiler-run').click();await page.waitForFunction(()=>document.querySelector('.compiler-check').dataset.state!=='pending',null,{timeout:35000});
  assert.equal(await page.locator('.compiler-check').getAttribute('data-state'),'passed',await page.locator('.compiler-status').innerText());
  assert.match(await page.locator('.compiler-stdout').innerText(),/Verified cases: 27/);
 }
 assert.deepEqual(errors,[]);console.log(`CUDA lesson: C++ reference, correct/broken mappings, shared bilingual controls and mobile layouts passed${live?'; live Godbolt output verified':''}.`);
}finally{await browser.close();server.close();}
