import {chromiumExecutable} from './browser-options.mjs';
import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import {serve} from '../scripts/serve.mjs';
const server=serve(4194);
const browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));
 await page.addInitScript(()=>{
  localStorage.setItem('feng-language','both');
  Object.defineProperty(navigator,'clipboard',{value:{writeText:async text=>{if(window.failCopy)throw Error('Denied');window.copiedMath=text;}}});
 });
 for(const slug of ['matrix-multiplication','band-storage-gaxpy','spike_notes','fast-matrix-vector-products']){
  await page.goto(`http://localhost:4194/blog/${slug}.html`);
  await page.evaluate(()=>document.querySelectorAll('details').forEach(n=>n.open=true));
  const formula=page.locator('.formula-block').first(), button=formula.locator('.formula-copy');
  const latex=await formula.getAttribute('data-latex');
  assert.equal(await page.locator('.formula-inline .formula-copy').count(),0);
  for(const mode of ['en','zh','both']){
   await page.locator(`[data-language-choice=${mode}]`).click();
   await button.focus();await page.keyboard.press('Enter');
   await page.waitForFunction(()=>window.copiedMath!==undefined);
   assert.equal(await page.evaluate(()=>window.copiedMath),latex);
   assert.equal(await button.getAttribute('title'),mode==='en'?'Copied':mode==='zh'?'已复制':'Copied / 已复制');
   assert.equal(await button.locator('svg').count(),1);
   assert.ok(await button.locator('.copied-symbol').isVisible());
   assert.equal(await button.getAttribute('data-copied'),'true');
  }
  await page.evaluate(()=>window.failCopy=true);await button.click();
  assert.equal(await formula.locator('textarea').inputValue(),latex);
  assert.equal(await formula.locator('textarea').evaluate(n=>n.selectionEnd-n.selectionStart),latex.length);
  await page.evaluate(()=>window.failCopy=false);await button.click();
  assert.equal(await formula.locator('textarea').count(),0);
  for(const width of [1440,768,390,320]){
   await page.setViewportSize({width,height:1000});
   await page.evaluate(()=>new Promise(resolve=>requestAnimationFrame(()=>requestAnimationFrame(resolve))));
   assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${slug}/${width}: page overflow`);
   assert.ok(await page.locator('.formula-scroll').evaluateAll(nodes=>nodes.every(n=>getComputedStyle(n).overflowX==='auto')));
   assert.deepEqual(await page.locator('.formula-block').evaluateAll(nodes=>nodes.filter(n=>n.getBoundingClientRect().height>0).flatMap(block=>{
    const number=block.querySelector('svg[data-labels] g[id]');if(!number)return [{missingTag:true}];
    const n=number.getBoundingClientRect(),f=block.querySelector('g[data-mml-node=mlabeledtr]').getBoundingClientRect(),svg=block.querySelector('mjx-container>svg').getBoundingClientRect(),copy=block.querySelector('.formula-copy').getBoundingClientRect();
    // AMS aligns tags on the mathematical baseline, which need not be the ink-box midpoint.
    const contained=n.top>=svg.top-1&&n.bottom<=svg.bottom+1;
    return contained&&n.left>=f.right-1&&n.top>=copy.bottom&&block.querySelectorAll('.eq-number').length===0?[]:[{number:number.id,contained,right:n.left>=f.right-1,copyClear:n.top>=copy.bottom}];
   })),[],`${slug}/${width}: native AMS tags must stay inside the formula and clear of the copy button`);

   if(slug==='matrix-multiplication'){
    const alignment=await page.locator('.formula-block').evaluateAll(nodes=>nodes.map(block=>{
     const viewport=block.querySelector('.formula-scroll').getBoundingClientRect();
     const math=block.querySelector('mjx-container').getBoundingClientRect();
     const formula=block.querySelector('g[data-mml-node=mlabeledtr]').getBoundingClientRect();
     const number=block.querySelector('svg[data-labels] g[id]').getBoundingClientRect();
     return {id:block.dataset.equationNumber,fits:math.width<=viewport.width+1,
      centerError:Math.abs((formula.left+formula.right-math.left-math.right)/2),
      rightGap:math.right-number.right};
    }));
    for(const item of alignment){
     if(item.fits)assert.ok(item.centerError<2,`${slug}/${width}/${item.id}: center the expression within its block`);
     assert.ok(item.rightGap>=0&&item.rightGap<4,`${slug}/${width}/${item.id}: align the native tag to the right`);
    }
   }

   if([1440,390].includes(width)){
    await formula.scrollIntoViewIfNeeded();
    await page.screenshot({path:`/tmp/feng-math-${slug}-${width}.png`});
   }
  }

 }
 const offline=await browser.newPage({javaScriptEnabled:false});
 await offline.goto('http://localhost:4194/blog/spike_notes.html');
 assert.ok(await offline.locator('.formula-block svg').count()>60);
 assert.ok(await offline.locator('.formula-inline svg').count()>200);
 assert.deepEqual(errors,[]);
 console.log('AMS SVG, bilingual copy feedback, keyboard copying, clipboard fallback, mobile overflow, and no-JavaScript rendering passed.');
}finally{await browser.close();server.close();}
