import {chromium,expect} from '@playwright/test';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import fs from 'node:fs/promises';
import {createRequire} from 'node:module';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
import {addDays,dateInfo,DAY_MS,dailySceneAt} from '../src/scripts/daily-math.js';
import {loadDailyMath} from '../scripts/daily-math.mjs';
const {entries:scenes}=await loadDailyMath();
const server=serve(4236),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
const require=createRequire(import.meta.url),errors=[];
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'});
 page.on('pageerror',e=>errors.push(e.message));
 await page.addInitScript(()=>localStorage.setItem('feng-language','both'));
 const start=Date.parse('2026-10-06T08:00:00Z');
 await page.clock.install({time:new Date(start)});
 // Freeze before navigation so a slow CI page load cannot overtake the pause time.
 await page.clock.pauseAt(new Date(start+3600000));
 await page.goto('http://localhost:4236/');await page.waitForSelector('[data-rendered="true"]');
 await page.clock.runFor(32);
 const root=page.locator('[data-daily-math]');
 const currentCopy=page.locator('[data-math-copy]:not([hidden])');
 assert.equal(await page.locator('.home-hero .daily-math-note').count(),0,'mathematics is separate from the foreground introduction');
 assert.equal(await page.locator('.daily-math-background > canvas').count(),1);
 assert.equal(await page.locator('.daily-math-background > .daily-math-note').count(),1,'animation and explanation belong to the same background layer');
 assert.equal(await root.locator('button').count(),0,'the background has no playback control');
 const canvas=async()=>createHash('sha256').update(await page.locator('#surface-canvas').evaluate(c=>c.toDataURL())).digest('hex');
 async function setReducedMotion(reducedMotion){
  // emulateMedia updates the query before Chromium delivers its change event.
  // Wait on the real browser event before advancing the paused animation clock.
  await page.evaluate(()=>{
   window.dailyMathMotionChanged=false;
   matchMedia('(prefers-reduced-motion: reduce)').addEventListener('change',()=>{window.dailyMathMotionChanged=true;},{once:true});
  });
  await page.emulateMedia({reducedMotion});
  await expect.poll(()=>page.evaluate(()=>window.dailyMathMotionChanged),{
   timeout:5000,message:`Chromium must deliver the ${reducedMotion} media change before checking animation`,
  }).toBe(true);
 }
 const initial=await canvas();await page.clock.runFor(1000);assert.equal(await canvas(),initial,'paused image must stay still');
 await setReducedMotion('no-preference');await page.clock.runFor(800);
 assert.notEqual(await canvas(),initial,'animation plays automatically when reduced motion is disabled');
 await setReducedMotion('reduce');await page.clock.runFor(32);
 const still=await canvas();await page.clock.runFor(800);
 assert.equal(await canvas(),still,'changing the system preference to reduced motion stops the animation');
 const renders=new Set();
 await fs.mkdir('test-results/daily-math',{recursive:true});
 for(let i=0;i<scenes.length;i++){
  const timestamp=Date.parse(scenes[i].date+'T08:00:00Z'),{scene,date}=dailySceneAt(timestamp,scenes);
  await page.clock.setSystemTime(new Date(timestamp));await page.evaluate(()=>window.dispatchEvent(new Event('pageshow')));
  await page.waitForFunction(id=>document.querySelector('[data-daily-math]').dataset.scene===id,scene.id);
  await page.clock.runFor(32);
  assert.equal(await root.getAttribute('data-scene'),scene.id);assert.equal(await root.getAttribute('data-date'),date);
  assert.equal(await page.locator('[data-math-source]').getAttribute('href'),scene.source.url);
  assert.match(await page.locator('[data-math-title]').innerText(),new RegExp(scene.zh));
  assert.equal(await currentCopy.count(),1,'only one topic explanation is visible');
  assert.equal(await currentCopy.getAttribute('data-math-copy'),scene.id);
  for(const [attribute,key] of [['data-math-description','description'],['data-math-reading','reading']]){
   for(const lang of ['en','zh']){
    const paragraph=currentCopy.locator(`[${attribute}] [data-lang=${lang}]`);
    const source=scene[key+(lang==='en'?'En':'Zh')];
    assert.equal(await paragraph.evaluate(node=>{
     const clone=node.cloneNode(true);
     clone.querySelectorAll('mjx-container').forEach(math=>math.replaceWith(document.createTextNode(`\\(${math.getAttribute('aria-label')}\\)`)));
     return clone.textContent;
    }),source,'daily changes preserve both the explanation and inline mathematics');
    assert.equal(await paragraph.locator('mjx-container[display]').count(),0,'prose uses inline typesetting');
   }
  }
  const expression=page.locator('[data-math-expression]:visible');
  assert.equal(await expression.count(),1,'only the current topic has a visible equation');
  assert.equal(await expression.getAttribute('data-math-expression'),scene.id);
  assert.equal(await expression.getAttribute('data-latex'),scene.formula);
  assert.equal(await expression.locator('mjx-container[display="true"]:visible').getAttribute('aria-label'),scene.formula);
  renders.add(await canvas());
  await root.screenshot({path:`test-results/daily-math/${scene.id}.png`});
  await page.setViewportSize({width:320,height:568});await page.clock.runFor(32);
  const formula=page.locator('[data-math-formula]');
  assert.ok(await formula.evaluate(node=>{
   const images=[...node.querySelectorAll('mjx-container > svg')].filter(image=>image.getBoundingClientRect().width>0);
   if(images.length!==1)return false;
   const box=images[0].getBoundingClientRect(),container=node.getBoundingClientRect();
   return box.left>=container.left-1&&box.right<=container.right+1;
  }),`${scene.id}: the complete display equation fits a narrow screen`);
  await formula.screenshot({path:`test-results/daily-math/formula-${scene.id}-320.png`});
  assert.ok(await currentCopy.evaluate(node=>[...node.querySelectorAll('mjx-container')].every(math=>{
   const rect=math.getBoundingClientRect(),bounds=node.getBoundingClientRect(),image=math.querySelector('svg');
   return getComputedStyle(image).display==='inline-block'&&rect.left>=bounds.left-1&&rect.right<=bounds.right+1;
  })),`${scene.id}: inline formulas stay in the text flow without overflowing`);
  await currentCopy.screenshot({path:`test-results/daily-math/inline-${scene.id}-320.png`});
  await page.setViewportSize({width:1440,height:1000});await page.clock.runFor(32);
 }
 assert.equal(renders.size,scenes.length,'each day must render a different construction');
 const boundary=Date.parse('2026-10-13T16:00:00Z');
 await page.clock.setSystemTime(new Date(boundary-1000));await page.evaluate(()=>window.dispatchEvent(new Event('pageshow')));
 await page.waitForFunction(()=>document.querySelector('[data-daily-math]').dataset.scene==='euler');
 await page.clock.runFor(1100);
 await page.waitForFunction(()=>document.querySelector('[data-daily-math]').dataset.scene==='taylor');
 assert.equal(await root.getAttribute('data-date'),'2026-10-14','open page changes at midnight without a reload');
 assert.equal(await root.getAttribute('data-scene'),dailySceneAt(boundary,scenes).scene.id);
 assert.equal(await page.locator('[data-math-expression]:visible').getAttribute('data-math-expression'),dailySceneAt(boundary,scenes).scene.id,'the display equation updates at midnight');
 assert.equal(await currentCopy.getAttribute('data-math-copy'),dailySceneAt(boundary,scenes).scene.id,'inline mathematics changes with the midnight topic');
 // Check actual screen shapes, including short landscape and high-DPI/4K displays.
 assert.equal(await page.locator('.daily-math-space').count(),0,'background must not reserve an illustration panel');
 for(const [width,height] of [[3840,2160],[2560,1440],[1920,1080],[1440,900],[1366,768],[1024,768],[960,720],[768,1024],[390,844],[320,568],[844,390]]){
  await page.setViewportSize({width,height});
  for(const mode of ['en','zh','both']){
   await page.evaluate(()=>window.scrollTo({top:0,behavior:'instant'}));await page.clock.runFor(32);
   await page.locator(`[data-language-choice="${mode}"]`).click();await page.clock.runFor(32);
   assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${width}×${height} ${mode}: page overflow`);
   const title=page.locator('[data-math-title]');
   assert.equal(await title.locator('[data-lang="en"]').isVisible(),mode!=='zh');
   assert.equal(await title.locator('[data-lang="zh"]').isVisible(),mode!=='en');
   assert.equal(await currentCopy.locator('[data-math-reading] [data-lang=en]').isVisible(),mode!=='zh');
   assert.equal(await currentCopy.locator('[data-math-reading] [data-lang=zh]').isVisible(),mode!=='en');
   if(width<=960){
    const mathBox=await page.locator('.daily-math-note').boundingBox(),profileBox=await page.locator('.hero-profile').boundingBox();
    assert.ok(mathBox.y+mathBox.height<=profileBox.y+1,'narrow screens prioritize the mathematical explanation before the author profile');
    const sourceBox=await page.locator('[data-math-source]').boundingBox(),nameBox=await page.locator('.home-intro-name').boundingBox();
    assert.ok(sourceBox.y+sourceBox.height+24<=nameBox.y,'leave readable spacing between the annotations and the profile');
    const principle=await currentCopy.locator('[data-math-description]').boundingBox();
    if(height>=720)assert.ok(principle.y<height-60,'the principle starts within the opening screen');
   }
   const geometry=await page.evaluate(()=>{
    const stage=document.querySelector('[data-daily-math]'),canvas=stage.querySelector('canvas');
    const rect=node=>{const r=node.getBoundingClientRect();return {x:r.x,y:r.y,width:r.width,height:r.height,bottom:r.bottom,right:r.right};};
    return {stage:rect(stage),canvas:rect(canvas),note:rect(stage.querySelector('.daily-math-note')),position:getComputedStyle(canvas).position,pointerEvents:getComputedStyle(canvas).pointerEvents,pixels:canvas.width*canvas.height};
   });
   assert.deepEqual(geometry.canvas,geometry.stage,`${width}×${height} ${mode}: canvas must fill the complete hero`);
   assert.equal(geometry.position,'absolute');assert.equal(geometry.pointerEvents,'none');
   assert.ok(geometry.note.bottom<=geometry.stage.bottom+1,'description must stay inside the hero');
   assert.ok(geometry.note.right<=width+1,'description must not be clipped');
   assert.ok(geometry.pixels<=8010000,'backing canvas has a bounded memory footprint');
   if(width>=1366&&height>=768&&mode!=='both')assert.ok(geometry.stage.height<=height-geometry.stage.y+1,'single-language desktop hero fits the first screen');
  }
  await root.screenshot({path:`test-results/daily-math/layout-${width}x${height}.png`});
  await page.locator('[data-math-source]').click({trial:true});
 }
 // Descriptions and source links remain accessible above the decorative canvas.
 await page.setViewportSize({width:1440,height:1000});await page.clock.runFor(32);
 await page.clock.resume();
 await page.addScriptTag({path:require.resolve('axe-core/axe.min.js')});
 assert.deepEqual(await page.evaluate(async()=>{const result=await axe.run(document.querySelector('[data-daily-math]'),{runOnly:{type:'tag',values:['wcag2a','wcag2aa','wcag21aa']}});return result.violations;}),[]);
 await page.locator('.hero-actions .lime-button').click();assert.ok(page.url().endsWith('#featured'));
 assert.deepEqual(errors,[]);
 // End of stock: clear the old topic rather than loop it or call it today's entry.
 const exhausted=addDays(scenes.at(-1).date,1);
 await page.clock.setSystemTime(new Date(exhausted+'T04:00:00Z'));await page.evaluate(()=>window.dispatchEvent(new Event('pageshow')));
 await page.waitForFunction(date=>document.querySelector('[data-daily-math]').dataset.date===date,exhausted);
 assert.equal(await root.getAttribute('data-scene'),'');assert.equal(await page.locator('[data-math-expression]').count(),0);
 assert.ok(await page.locator('.daily-math-note a[href="/math/"]').isVisible());
 // A failed request must not leave the previous day's content on screen.
 await page.route('**/assets/daily-math/2026-10-25.json',route=>route.abort());
 await page.clock.setSystemTime(new Date('2026-10-25T04:00:00Z'));await page.evaluate(()=>window.dispatchEvent(new Event('pageshow')));
 await page.waitForFunction(()=>document.querySelector('[data-daily-math]').dataset.date==='2026-10-25');
 assert.equal(await root.getAttribute('data-scene'),'');
 await page.unroute('**/assets/daily-math/2026-10-25.json');
 await page.goto('http://localhost:4236/math/2026-10-06/');await page.waitForFunction(()=>document.querySelector('[data-daily-math]').dataset.scene==='fourier');
 await page.clock.setSystemTime(new Date('2028-02-29T04:00:00Z'));await page.evaluate(()=>window.dispatchEvent(new Event('pageshow')));
 assert.equal(await page.locator('[data-daily-math]').getAttribute('data-date'),'2026-10-06','an archived entry never follows the current day');
 await page.clock.setSystemTime(new Date('2026-10-08T04:00:00Z'));await page.goto('http://localhost:4236/math/');
 assert.equal(await page.locator('[data-entry-date]:visible').count(),3,'archive reveals only dates already published');
 assert.equal(await page.locator('[data-entry-date="2026-10-09"]').isVisible(),false);
 await page.goto('http://localhost:4236/');await page.waitForFunction(()=>document.querySelector('[data-daily-math]').dataset.scene==='standing');
 assert.equal(await page.locator('[data-math-copy]').count(),1,'homepage loads one topic, independent of library size');
 // JavaScript-disabled readers retain a useful topic and source.
 const staticPage=await browser.newPage({javaScriptEnabled:false});await staticPage.goto('http://localhost:4236/');
 if(scenes.some(t=>t.date===dateInfo().date)){
  assert.ok(await staticPage.locator('[data-math-copy] [data-math-description]').textContent());
  assert.equal(await staticPage.locator('[data-math-expression]:visible mjx-container[display="true"]:visible').count(),1,'display typesetting works without JavaScript');
 }else assert.ok(await staticPage.locator('.daily-math-note a[href="/math/"]').isVisible());
 console.log('Daily mathematics passed: all scheduled drawings, midnight rollover, autoplay/reduced motion, bilingual responsive layout, accessibility, and no-JS fallback.');
}finally{await browser.close();server.close();}
