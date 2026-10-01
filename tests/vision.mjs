import {chromium} from '@playwright/test';
import sharp from 'sharp';
import assert from 'node:assert/strict';
import {createRequire} from 'node:module';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const require=createRequire(import.meta.url);
const server=serve(4227),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
const url='http://localhost:4227/blog/human-visual-system.html';
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[],requests=[];
 page.on('pageerror',e=>errors.push(e.message));page.on('request',r=>requests.push(r.url()));
 await page.addInitScript(()=>localStorage.setItem('feng-language','both'));
 await page.goto(url);
 assert.equal(await page.locator('#vision-essay h2').count(),8);
 assert.equal(await page.locator('#vision-essay .v-plate').count(),6);
 assert.equal(await page.locator('#vision-essay canvas').count(),1);
 assert.ok(!requests.some(u=>/three\.(core|module)\.min\.js/.test(u)),'3D renderer is deferred until the figure is near the viewport');
 assert.equal(await page.locator('#vision-essay .compiler-editor').count(),2);
 assert.ok(!requests.some(u=>/\/vision(?:-anatomy|-renderer|-models)?\.js/.test(u)),'static plates must not load obsolete 3D modules');
 // A real WebGL model: independent views, keyboard camera and anatomical selection.
 const eye=page.locator('[data-eye-viewer]');await eye.scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelector('[data-eye-viewer]')?.eyeViewer,{},{timeout:120000});
 assert.equal(await eye.getAttribute('data-eye-state'),'ready');
 const cut=await page.locator('.eye-stage canvas').screenshot();
 assert.ok((await sharp(cut).stats()).channels.slice(0,3).some(c=>c.stdev>25),'model must render visible geometry and shading');
 await page.locator('[data-eye-mode=whole]').click();
 await page.waitForFunction(()=>document.querySelector('[data-eye-viewer]').dataset.mode==='whole');
 const whole=await page.locator('.eye-stage canvas').screenshot();assert.notDeepEqual(whole,cut,'intact and cutaway views must differ');
 await page.locator('[data-eye-part=lens]').click();
 await page.waitForFunction(()=>document.querySelector('[data-eye-viewer]').dataset.mode==='cut');
 assert.equal(await page.locator('[data-eye-mode=cut]').getAttribute('aria-pressed'),'true');
 await page.locator('[data-eye-part=nerve]').click();
 await page.waitForFunction(()=>Number(document.querySelector('[data-eye-viewer]').dataset.yaw)<-2);
 assert.ok(await page.locator('.eye-part-title').textContent().then(s=>s.includes('Optic nerve')&&s.includes('视神经')));
 await page.locator('[data-eye-mode=layers]').click();
 await page.waitForFunction(()=>document.querySelector('[data-eye-viewer]').dataset.mode==='layers');
 const layers=await page.locator('.eye-stage canvas').screenshot();assert.notDeepEqual(layers,cut,'separated view must expose the structures');
 const canvas=page.locator('.eye-stage canvas');await canvas.focus();await canvas.press('Home');
 await page.waitForFunction(()=>document.querySelector('[data-eye-viewer]').dataset.mode==='cut');
 const initialYaw=Number(await eye.getAttribute('data-yaw'));await canvas.press('ArrowRight');
 await page.waitForFunction(yaw=>Number(document.querySelector('[data-eye-viewer]').dataset.yaw)>yaw,initialYaw);
 await canvas.press('+');await page.waitForFunction(()=>Number(document.querySelector('[data-eye-viewer]').dataset.zoom)>1);
 await canvas.press('Home');await page.waitForFunction(()=>document.querySelector('[data-eye-viewer]').dataset.zoom==='1.00');
 // A reader who scrolls away should not pay for a perpetual animation loop.
 await page.evaluate(()=>scrollTo(0,0));await page.waitForTimeout(300);
 const frames=await eye.getAttribute('data-frames');await page.waitForTimeout(300);assert.equal(await eye.getAttribute('data-frames'),frames);
 // All figure captions and prose are authored as one English/Chinese pair.
 const unpaired=await page.locator('#vision-essay p.parallel-text,#vision-essay .v-plate figcaption').evaluateAll(nodes=>nodes.filter(n=>n.querySelectorAll(':scope > [data-lang=en]').length!==1||n.querySelectorAll(':scope > [data-lang=zh]').length!==1).map(n=>n.textContent));
 assert.deepEqual(unpaired,[]);
 for(const width of [1440,1024,768,390,320]){
  await page.setViewportSize({width,height:1000});
  for(const language of ['en','zh','both']){
   await page.locator(`[data-language-choice=${language}]`).click();
   await page.evaluate(()=>new Promise(requestAnimationFrame));
   assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`page overflow ${width}/${language}`);
   const figures=await page.locator('#vision-essay .v-plate').evaluateAll(nodes=>nodes.map(n=>({fits:n.scrollWidth<=n.clientWidth+1,hasPlot:[...n.querySelectorAll('.v-svg')].every(s=>s.getBoundingClientRect().width>150)})));
   assert.ok(figures.every(x=>x.fits&&x.hasPlot),`figure overflow ${width}/${language}`);
   const clippedLabels=await page.locator('#vision-essay .v-svg').evaluateAll(nodes=>nodes.flatMap(svg=>{
    const v=svg.viewBox.baseVal;
    return [...svg.querySelectorAll('text')].filter(t=>{const b=t.getBBox();return b.x<v.x-3||b.y<v.y-3||b.x+b.width>v.x+v.width+3||b.y+b.height>v.y+v.height+3;}).map(t=>t.textContent);
   }));
   assert.deepEqual(clippedLabels,[],`figure labels clipped ${width}/${language}`);
   if(language!=='both'){
    const hidden=language==='en'?'zh':'en';
    assert.equal(await page.locator(`#vision-essay .v-plate [data-lang=${hidden}]:visible`).count(),0);
   }
  }
 }
 await page.setViewportSize({width:1440,height:1000});await page.locator('[data-language-choice=zh]').click();
 for(const id of ['vision-anatomy','vision-physics','vision-cell','vision-color','vision-contrast','vision-events'])await page.locator('#'+id).screenshot({path:`/tmp/feng-${id}-desktop.png`});
 await page.setViewportSize({width:390,height:844});await page.locator('[data-language-choice=both]').click();
 await page.locator('.eye-stage').scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>{const r=document.querySelector('[data-eye-viewer]'),b=r.querySelector('.eye-stage').getBoundingClientRect();return Math.abs(r.eyeViewer.camera.aspect-b.width/b.height)<.001;});
 const mobile=await page.locator('.eye-stage').screenshot({path:'/tmp/feng-eye-mobile-stage.png'});
 assert.ok((await sharp(mobile).stats()).channels.slice(0,3).some(c=>c.stdev>25),'mobile resizing must retain rendered anatomy');
 await page.screenshot({path:'/tmp/feng-eye-mobile-viewport.png'});
 await page.evaluate(()=>document.activeElement?.blur());
 await page.locator('#vision-anatomy').screenshot({path:'/tmp/feng-vision-anatomy-mobile.png'});
 await page.locator('#vision-contrast').screenshot({path:'/tmp/feng-vision-contrast-mobile.png'});
 await page.addScriptTag({path:require.resolve('axe-core/axe.min.js')});
 const audit=await page.evaluate(()=>axe.run(document.getElementById('vision-essay'),{runOnly:{type:'tag',values:['wcag2a','wcag2aa','wcag21aa']}}));
 assert.deepEqual(audit.violations.map(v=>({id:v.id,targets:v.nodes.map(n=>n.target)})),[]);
 const noJS=await browser.newPage({javaScriptEnabled:false,viewport:{width:390,height:844}});await noJS.goto(url);
 assert.equal(await noJS.locator('#vision-essay .v-plate').count(),6);
 await noJS.locator('.eye-poster').scrollIntoViewIfNeeded();await noJS.waitForFunction(()=>document.querySelector('.eye-poster').complete);
 assert.equal(await noJS.locator('.eye-poster').evaluate(img=>img.complete&&img.naturalWidth>0),true);
 assert.equal(await noJS.locator('.eye-tools').isVisible(),false);
 for(const figure of await noJS.locator('#vision-essay .v-plate').all())assert.equal(await figure.isVisible(),true);
 assert.equal(await noJS.locator('#vision-essay pre code[data-godbolt]').count(),2);
 const duplicateIds=await page.locator('[id]').evaluateAll(nodes=>{const seen=new Set();return nodes.map(n=>n.id).filter(id=>seen.has(id)||!seen.add(id));});
 assert.deepEqual(duplicateIds,[]);
 assert.deepEqual(errors,[]);
 // A lost GPU context leaves a useful, locally rendered anatomical plate.
 await eye.scrollIntoViewIfNeeded();
 await eye.evaluate(r=>r.eyeViewer.renderer.getContext().getExtension('WEBGL_lose_context').loseContext());
 await page.waitForFunction(()=>document.querySelector('[data-eye-viewer]').dataset.eyeState==='fallback');
 assert.equal(await page.locator('.eye-poster').isVisible(),true);assert.equal(await page.locator('.eye-tools').isVisible(),false);
 const fallback=await browser.newPage({viewport:{width:390,height:844}});
 await fallback.addInitScript(()=>{const get=HTMLCanvasElement.prototype.getContext;HTMLCanvasElement.prototype.getContext=function(type,...args){return /^webgl|experimental-webgl/.test(type)?null:get.call(this,type,...args);};});
 await fallback.goto(url);await fallback.locator('[data-eye-viewer]').scrollIntoViewIfNeeded();
 await fallback.waitForFunction(()=>document.querySelector('[data-eye-viewer]').dataset.eyeState==='fallback');
 assert.equal(await fallback.locator('.eye-poster').isVisible(),true);await fallback.close();
 console.log('Vision article passed: 3D anatomy/camera/separation/GPU fallback, six static plates, bilingual reading, 5 widths, accessible labels/contrast, no-JavaScript figures, and shared code editors.');
}finally{await browser.close();server.close();}
