import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const server=serve(4211),browser=await chromium.launch({executablePath:chromiumExecutable(),args:['--no-sandbox']});
try{
 const page=await browser.newPage({viewport:{width:1440,height:900},reducedMotion:'reduce'}),errors=[];
 page.on('pageerror',e=>errors.push(e.message));await page.addInitScript(()=>localStorage.setItem('feng-language','both'));await page.clock.install();
 const cases=[['frank-wolfe-algorithm','fw'],['pytorch-02-tensor-strides-storage','tensor'],['band-storage-gaxpy','band'],['symmetric-storage-gaxpy','symmetric'],['cuda-rust-two-tracks-blog','cuda'],['spike_notes','lif'],['waves-and-phase','waves']];
 for(const [slug,name]of cases){
  await page.goto(`http://localhost:4211/blog/${slug}.html`);const figure=page.locator(`[data-process=${name}]`),root=page.locator(name==='waves'?'[data-wave-lab]':'#'+({fw:'fw-lab',tensor:'tensor-layout',band:'band-lab',symmetric:'symmetric-lab',cuda:'cuda-map',lif:'lab-lif'})[name]);
  assert.equal(await root.getAttribute('data-demo-frame'),'complete');
  assert.ok((await root.boundingBox()).height<=820,`${name} stays compact in bilingual desktop mode`);
  const step=root.locator('.reading-playback button[data-icon=step]');
  if(name==='fw'){
   for(let i=0;i<4;i++)await step.click();
   assert.equal(await figure.getAttribute('data-iterate'),'0.800000');assert.equal(await figure.getAttribute('data-gamma'),'0.800000');
  }
  if(name==='lif'){
   await step.click();assert.ok(Math.abs(Number(await figure.getAttribute('data-potential'))-1)<1e-12,'integration reaches threshold before reset');
   await step.click();assert.equal(await figure.getAttribute('data-reset'),'false');
   const t=await figure.getAttribute('data-time');await step.click();assert.equal(await figure.getAttribute('data-time'),t,'reset takes no physical time');assert.equal(await figure.getAttribute('data-potential'),'0');
   await root.locator('> .reading-explore > summary').click();await page.locator('#lif-input').fill('0.8');
   assert.equal(await figure.getAttribute('data-spiking'),'false');assert.equal(await figure.getAttribute('data-reset'),'false');assert.ok(Number(await figure.getAttribute('data-potential'))<.8);
   await page.locator('#lif-input').fill('1.6');await root.locator('> .reading-explore > summary').click();
  }
  if(name==='cuda'){
   const good=JSON.parse(await figure.getAttribute('data-correct-writes')),bad=JSON.parse(await figure.getAttribute('data-broken-writes'));
   assert.deepEqual(good,Array(17).fill(1));assert.deepEqual(bad,[...Array(8).fill(3),...Array(9).fill(0)]);
  }
  for(const width of [1440,390,320]){await page.setViewportSize({width,height:900});for(const lang of ['en','zh','both']){
   await page.locator(`[data-language-choice=${lang}]`).click();assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${name}/${width}/${lang}`);
   if(lang==='en')assert.doesNotMatch(await figure.innerText(),/[\u3400-\u9fff]/);
   const controls=await root.locator('.reading-playback button').evaluateAll(buttons=>buttons.map(b=>({text:b.innerText,label:b.getAttribute('aria-label'),title:b.title,en:b.dataset.ariaLabelEn,zh:b.dataset.ariaLabelZh,width:b.getBoundingClientRect().width,height:b.getBoundingClientRect().height,icons:b.querySelectorAll('svg[aria-hidden=true]').length})));
   for(const b of controls){assert.equal(b.text,'');assert.equal(b.icons,1);assert.equal(b.width,44);assert.equal(b.height,44);assert.equal(b.label,lang==='en'?b.en:lang==='zh'?b.zh:`${b.en} / ${b.zh}`);assert.equal(b.title,b.label);}

  }}
  await page.setViewportSize({width:1440,height:900});
 }
 // Continuous interpolation, cancellation, and pause/resume use visible geometry.
 await page.emulateMedia({reducedMotion:'no-preference'});await page.goto('http://localhost:4211/blog/pytorch-02-tensor-strides-storage.html');
 const tensor=page.locator('[data-process=tensor]'),root=page.locator('#tensor-layout');
 await tensor.scrollIntoViewIfNeeded();await page.waitForFunction(()=>document.querySelector('#tensor-layout').dataset.demoPlaying==='true');
 await page.clock.runFor(350);const before=await tensor.locator('[data-token="6"]').evaluate(n=>n.style.transform);
 const storageBefore=await tensor.locator('[data-buffer-a="6"]').boundingBox();await page.clock.runFor(500);
 const after=await tensor.locator('[data-token="6"]').evaluate(n=>n.style.transform);assert.notEqual(after,before,'logical values move continuously');assert.deepEqual(await tensor.locator('[data-buffer-a="6"]').boundingBox(),storageBefore,'backing storage never moves during transpose');
 const play=root.locator('.reading-playback button').first();await play.click();assert.match(await play.getAttribute('aria-label'),/Resume/);const paused=await tensor.locator('[data-token="6"]').evaluate(n=>n.style.transform);await page.clock.runFor(2000);assert.equal(await tensor.locator('[data-token="6"]').evaluate(n=>n.style.transform),paused,'pause freezes the in-flight value');
 await play.click();await page.clock.runFor(7000);assert.equal(await root.getAttribute('data-demo-frame'),'complete');assert.equal(await tensor.locator('.process-token').count(),0);
 assert.equal(await tensor.locator('[data-buffer-a="6"] b').textContent(),'6');assert.equal(await tensor.locator('[data-buffer-b="7"] b').textContent(),'6');
 const staticPage=await browser.newPage({javaScriptEnabled:false});
 for(const [slug,name]of cases){await staticPage.goto(`http://localhost:4211/blog/${slug}.html`);assert.ok(await staticPage.locator(`[data-process=${name}]`).isVisible());}
 assert.deepEqual(errors,[]);console.log('Process diagrams: seven models, exact reset timing, copy addresses, GPU collisions, continuous motion, pause/resume, compact bilingual layouts and static fallback passed.');
}finally{await browser.close();server.close();}
