import {chromiumExecutable} from './browser-options.mjs';
import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {execFileSync} from 'node:child_process';
import {load} from 'cheerio';
import {serve} from '../scripts/serve.mjs';
const live=process.argv.includes('--live');
const slugs=['band-storage-gaxpy','symmetric-storage-gaxpy','fast-matrix-vector-products','frank-wolfe-algorithm','matrix-multiplication'];
const dates=['2021-07-19','2021-07-19','2019-01-12','2019-01-04','2019-01-07'];
const temp=await fs.mkdtemp(path.join(os.tmpdir(),'feng-lessons-'));
let programs=0;
try{
 for(const [index,slug]of slugs.entries()){
  const $=load(await fs.readFile(`content/posts/${slug}.html`,'utf8'));
  const metadata=JSON.parse($('#article-metadata').text());assert.equal(metadata.date,dates[index]);assert.equal(metadata.updated,'2026-09-23');
  assert.equal($('pre code:not([data-godbolt="c++"])').length,0,slug+' contains unconverted code');
  assert.equal($('[data-lang] [data-math]:not([data-display="inline"])').length,0,'share equations');
  assert.equal($('[data-lang] input,[data-lang] button,[data-lang] canvas').length,0,'share controls');
  for(const code of $('code[data-godbolt="c++"]').toArray()){
   const file=path.join(temp,`${programs}.cpp`),bin=path.join(temp,`${programs}`);await fs.writeFile(file,$(code).text());
   execFileSync(process.env.CXX||'g++',['-std=c++17','-O2','-Wall','-Wextra',file,'-o',bin],{timeout:30000});
   const output=execFileSync(bin,{encoding:'utf8',timeout:10000});assert.ok(output.trim(),slug+' produced no output');
   console.log(slug,output.trim().split('\n').join(' | '));programs++;
  }
 }
 assert.equal(programs,6);
}finally{await fs.rm(temp,{recursive:true,force:true});}
const server=serve(4188),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox','--no-proxy-server']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'}),errors=[];page.on('pageerror',error=>errors.push(error.message));await page.addInitScript(()=>localStorage.setItem('feng-language','both'));
 for(const slug of slugs){
  await page.goto(`http://localhost:4188/blog/${slug}.html`);
  assert.ok(await page.locator('.article-toc nav a').count()>=5);
  assert.equal(await page.locator('.lesson-lab').count(),1);
  for(const width of[1440,768,390,320]){await page.setViewportSize({width,height:1000});for(const language of['en','zh','both']){await page.locator(`[data-language-choice=${language}]`).click();assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${slug}: overflow at ${width}/${language}`);}}
  await page.setViewportSize({width:1440,height:1000});await page.locator('[data-language-choice=en]').click();
  if(await page.locator('.lesson-lab > .reading-explore > summary').count())await page.locator('.lesson-lab > .reading-explore > summary').click();
  if(slug==='band-storage-gaxpy'){
   await page.locator('#band-p').fill('2');await page.locator('#band-q').fill('1');await page.locator('#band-i').fill('3');await page.locator('#band-j').fill('2');
   assert.match(await page.locator('#band-result').innerText(),/offset 10/);assert.equal(await page.locator('#band-packed .selected').innerText(),'-1');
   await page.locator('#band-i').fill('0');await page.locator('#band-j').fill('5');assert.match(await page.locator('#band-result').innerText(),/outside the band/);assert.equal(await page.locator('#band-packed .selected').count(),0);
   await page.locator('#band-p').fill('5');await page.locator('#band-q').fill('5');assert.match(await page.locator('#band-result').innerText(),/visited positions: 36/);
  }else if(slug==='symmetric-storage-gaxpy'){
   assert.match(await page.locator('#sym-result').innerText(),/offset 5/);assert.equal(await page.locator('#sym-dense .selected').count(),2);
   await page.locator('#sym-i').fill('3');await page.locator('#sym-j').fill('3');assert.match(await page.locator('#sym-result').innerText(),/offset 9/);assert.equal(await page.locator('#sym-dense .selected').count(),1);
  }else if(slug==='fast-matrix-vector-products'){
   const values=()=>page.locator('#fft-spectrum strong').allTextContents();assert.deepEqual(await values(),['0.0','4.0','0.0','0.0','0.0','0.0','0.0','4.0']);
   await page.locator('#fft-impulse').click();assert.deepEqual(await values(),Array(8).fill('1.0'));await page.locator('#fft-constant').click();assert.deepEqual(await values(),['8.0',...Array(7).fill('0.0')]);
   await page.locator('#fft-alternating').click();assert.equal((await values())[4],'8.0');await page.locator('#fft-sine').click();await page.locator('#fft-bin').fill('3');assert.equal((await values())[3],'4.0');assert.equal((await values())[5],'4.0');
  }else if(slug==='frank-wolfe-algorithm'){
   assert.equal(await page.locator('#fw-steps').inputValue(),'8');await page.locator('#fw-reset').click();
   await page.locator('#fw-next').click();assert.match(await page.locator('#fw-result').innerText(),/x = \(0.8000, 0.0000\)/);
   await page.locator('#fw-steps').fill('40');const numbers=(await page.locator('#fw-result').innerText()).match(/f = ([0-9.]+).*G = ([0-9.]+)/s);assert.ok(numbers);assert.ok(+numbers[1]>=.04 && +numbers[1]-.04<=+numbers[2]+1e-6);
   await page.locator('#fw-inside').click();assert.equal(await page.locator('#fw-steps').inputValue(),'0');await page.locator('#fw-next').click();assert.match(await page.locator('#fw-result').innerText(),/x = \(0.0000, 0.3500\)/);
  }else{
   await page.locator('#matrix-step').fill('7');assert.match(await page.locator('#matrix-result').innerText(),/trace: 8/);await page.locator('#matrix-row').click();assert.match(await page.locator('#matrix-result').innerText(),/trace: 1/);assert.match(await page.locator('#matrix-result').innerText(),/0 → 1 → 2 → 3 → 4 → 5 → 6 → 7/);
  }
  const inputs=await page.locator('.lesson-lab input').evaluateAll(nodes=>nodes.map(n=>n.value));await page.locator('[data-language-choice=zh]').click();assert.deepEqual(await page.locator('.lesson-lab input').evaluateAll(nodes=>nodes.map(n=>n.value)),inputs);await page.locator('[data-language-choice=both]').click();
  await page.locator('.lesson-lab').screenshot({path:path.join(os.tmpdir(),`feng-lesson-${slug}.png`)});
  if(live){for(const panel of await page.locator('.compiler-check').all()){
   await panel.locator('.compiler-run').click();await page.waitForFunction(()=>![...document.querySelectorAll('.compiler-check')].some(n=>n.dataset.state==='pending'),{},{timeout:35000});
   assert.equal(await panel.getAttribute('data-state'),'passed',await panel.innerText());assert.ok((await panel.locator('.compiler-stdout').innerText()).trim());
  }}
 }
 assert.deepEqual(errors,[]);console.log(`${programs} C++ programs and 5 interactive lessons verified${live?' including live Godbolt runs':''}.`);
}finally{await browser.close();server.close();}
