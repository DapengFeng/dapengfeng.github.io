import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const server=serve(4210),browser=await chromium.launch({executablePath:chromiumExecutable(),args:['--no-sandbox']});
const errors=[];
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'});
 page.on('pageerror',e=>errors.push(e.message));await page.addInitScript(()=>localStorage.setItem('feng-language','en'));
 await page.clock.install();
 // CSS smooth scrolling runs on the compositor, outside Playwright's fake clock.
 // Wait for an actual offscreen observation before advancing animation timers.
 const scrollOffscreen=async selector=>{
  await page.evaluate(selector=>{
   window.testOffscreenObserved=false;
   const observer=new IntersectionObserver(([entry])=>{
    if(!entry.isIntersecting){observer.disconnect();window.testOffscreenObserved=true;}
   });
   observer.observe(document.querySelector(selector));
   scrollTo({top:0,left:0,behavior:'instant'});
  },selector);
  await page.waitForFunction(()=>window.testOffscreenObserved===true);
 };
 for(const [slug,id]of [['matrix-multiplication','matrix-lab'],['frank-wolfe-algorithm','fw-lab'],['pytorch-01-what-is-pytorch','pt-route'],['pytorch-02-tensor-strides-storage','tensor-layout'],['rust-vs-cpp-blog','memory-map']]){
  await page.goto(`http://localhost:4210/blog/${slug}.html`);
  const root=id==='memory-map'?page.locator('.memory-demo'):page.locator('#'+id);
  assert.equal(await root.getAttribute('data-demo-frame'),'complete');
  await root.locator('.reading-playback button').first().scrollIntoViewIfNeeded();
  await page.clock.runFor(20000);assert.equal(await root.getAttribute('data-demo-frame'),'complete','reduced motion retains complete default');
  await root.locator('.reading-playback button').first().click();
  const y=await page.evaluate(()=>scrollY);
  await page.clock.runFor(18000);
  assert.equal(await root.getAttribute('data-demo-frame'),'complete',id+' finishes');
  assert.ok(Math.abs(await page.evaluate(()=>scrollY)-y)<2,id+' must not scroll the page');
  if(id==='pt-route'){
   assert.equal(await root.getAttribute('data-phase'),'backward');
   assert.deepEqual(JSON.parse(await root.getAttribute('data-grad-a')),[15,19,23,15,19,23]);
  }
  if(id==='matrix-lab')assert.deepEqual(await root.locator('#matrix-comparison strong').allTextContents(),['1 distinct cache lines1 个不同缓存行','8 distinct cache lines8 个不同缓存行']);
 }
 await page.goto('http://localhost:4210/blog/pytorch-01-what-is-pytorch.html');
 assert.deepEqual(await page.locator('#pt-matrix-y .pt-cell').allTextContents(),['58','64','139','154']);
 await page.locator('#pt-multiply > .reading-explore > summary').click();
 await page.locator('#pt-a00').fill('-2');assert.deepEqual(await page.locator('#pt-matrix-y .pt-cell').allTextContents(),['37','40','139','154']);
 await page.locator('#pt-matrix-play').click();await page.clock.runFor(1800);
 await page.locator('#pt-a00').fill('3');await page.clock.runFor(7000);
 assert.deepEqual(await page.locator('#pt-matrix-y .pt-cell').allTextContents(),['72','80','139','154'],'reader input stops playback');
 await page.emulateMedia({reducedMotion:'no-preference'});
 await page.goto('http://localhost:4210/blog/matrix-multiplication.html');
 await page.locator('#matrix-comparison').scrollIntoViewIfNeeded();await page.waitForFunction(()=>document.querySelector('#matrix-lab').dataset.demoPlaying==='true');await page.clock.runFor(1500);
 assert.equal(await page.locator('#matrix-lab').getAttribute('data-demo-playing'),'true');
 await scrollOffscreen('#matrix-comparison');
 const paused=await page.locator('#matrix-lab').getAttribute('data-demo-frame');
 assert.notEqual(paused,'complete','pause check must start during playback');
 assert.ok(await page.locator('#matrix-comparison').evaluate(el=>el.getBoundingClientRect().top>=innerHeight),'diagram is below the viewport');
 await page.clock.runFor(5000);
 assert.equal(await page.locator('#matrix-lab').getAttribute('data-demo-frame'),paused,'offscreen pauses');
 await page.locator('#matrix-comparison').scrollIntoViewIfNeeded();await page.waitForFunction(old=>document.querySelector('#matrix-lab').dataset.demoFrame!==old,paused);await page.clock.runFor(8000);
 assert.equal(await page.locator('#matrix-lab').getAttribute('data-demo-frame'),'complete');
 await scrollOffscreen('#matrix-comparison');await page.locator('#matrix-comparison').scrollIntoViewIfNeeded();await page.clock.runFor(3000);
 assert.equal(await page.locator('#matrix-lab').getAttribute('data-demo-frame'),'complete','autoplay only once');
 await page.locator('.reading-playback button').first().click();await page.clock.runFor(800);await page.emulateMedia({reducedMotion:'reduce'});
 await page.waitForFunction(()=>document.querySelector('#matrix-lab').dataset.demoFrame==='complete');
 assert.equal(await page.locator('#matrix-lab').getAttribute('data-demo-frame'),'complete','changing motion preference stops and restores full diagram');
 for(const slug of ['spike_notes','fast-matrix-vector-products','pytorch-01-what-is-pytorch','pytorch-02-tensor-strides-storage']){
  await page.goto(`http://localhost:4210/blog/${slug}.html`);
  for(const width of [1440,390,320]){await page.setViewportSize({width,height:1000});for(const lang of ['en','zh','both']){
   await page.locator(`[data-language-choice=${lang}]`).click();
   assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${slug}/${width}/${lang}`);
   if(lang==='en')assert.doesNotMatch((await page.locator('.reading-comparison,.reading-overview,[data-process]').allInnerTexts()).join(''),/[\u3400-\u9fff]/);
  }}
 }
 await page.setViewportSize({width:1440,height:900});
 for(const slug of ['spike_notes','matrix-multiplication','pytorch-01-what-is-pytorch','pytorch-02-tensor-strides-storage','fast-matrix-vector-products','frank-wolfe-algorithm','band-storage-gaxpy','cuda-rust-two-tracks-blog']){
  await page.goto(`http://localhost:4210/blog/${slug}.html`);await page.locator('[data-language-choice=both]').click();
  const sizes=await page.locator('.experiment,.lesson-lab,.pt-demo,.tl-demo').evaluateAll(nodes=>nodes.map(n=>[n.id,n.getBoundingClientRect().height]));
  for(const [id,height]of sizes)assert.ok(height<=820,`${id}: ${height}px exceeds the compact bilingual reading budget`);
 }
 const nojs=await browser.newPage({javaScriptEnabled:false});
 for(const slug of ['spike_notes','fast-matrix-vector-products','pytorch-01-what-is-pytorch','pytorch-02-tensor-strides-storage']){
  await nojs.goto(`http://localhost:4210/blog/${slug}.html`);assert.ok(await nojs.locator('.reading-comparison,.reading-overview,[data-process]').first().isVisible(),'static explanation needs no script');
 }
 assert.deepEqual(errors,[]);console.log('Reader diagrams: complete defaults, finite playback, backward gradients, input ownership, viewport pause, reduced motion, bilingual comparisons, mobile and no-JS passed.');
}finally{await browser.close();server.close();}
