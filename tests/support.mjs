import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import {createRequire} from 'node:module';
import {readFile} from 'node:fs/promises';
import {load} from 'cheerio';
import {serve} from '../scripts/serve.mjs';
import {loadSupport,supportDialog,supportTrigger} from '../scripts/support.mjs';
import {chromiumExecutable} from './browser-options.mjs';

const require=createRequire(import.meta.url),server=serve(4222);
const browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
const url='http://localhost:4222/blog/jiuzhaigou-water-and-mountains.html';
try{
 // Exercise the UI even when the build has no payment variable (for example, a fork).
 // Only the intercepted test page gets this fixture; dist and the deployment stay untouched.
 const profile='https://paypal.me/example';
 const fixture=load(await readFile('dist/blog/jiuzhaigou-water-and-mountains.html','utf8'),{scriptingEnabled:false});
 fixture('#article-support-dialog,[data-support-open],noscript:has(.support-fallback)').remove();
 fixture('[data-share-open]').after(supportTrigger());
 fixture('.article-reading-main').append(supportDialog(loadSupport({PAYPAL_ME_URL:profile})));
 const fixtureHTML=fixture.html();
 const page=await browser.newPage({viewport:{width:1440,height:1000}}),errors=[],requests=[];
 await page.route(url,r=>r.fulfill({contentType:'text/html',body:fixtureHTML}));
 page.on('pageerror',e=>errors.push(e.message));
 page.on('request',r=>requests.push(r.url()));
 await page.route('https://giscus.app/**',r=>r.fulfill({contentType:'text/html',body:'<!doctype html><title>Mock discussion</title>'}));
 await page.addInitScript(()=>localStorage.setItem('feng-language','zh'));
 await page.goto(url);
 const dialog=page.locator('#article-support-dialog'),trigger=page.locator('[data-support-open]').first();
 const dock=page.locator('.article-action-dock');
 assert.equal(await page.locator('[data-support-open]').count(),1,'one floating support entry');
 assert.equal(await page.locator('[data-share-open]').count(),1,'one floating sharing entry');
 const beforeScroll=await dock.boundingBox();
 await page.evaluate(()=>scrollTo({top:1200,behavior:'instant'}));
 const afterScroll=await dock.boundingBox();
 assert.ok(Math.abs(beforeScroll.y-afterScroll.y)<1,'toolbar stays fixed while the article scrolls');
 for(const width of [1440,768,390,320]){
  await page.setViewportSize({width,height:800});
  const box=await dock.boundingBox();
  assert.ok(box.x>=0&&box.x+box.width<=width&&box.y+box.height<=800,'toolbar fits viewport');
  const boxes=await dock.locator('button:visible').evaluateAll(nodes=>nodes.map(n=>{const r=n.getBoundingClientRect();return {x:r.x,width:r.width,height:r.height};}));
  assert.equal(boxes.length,3);
  for(let i=0;i<boxes.length;i++){
   assert.equal(boxes[i].width,44);assert.equal(boxes[i].height,44);
   if(i)assert.ok(boxes[i].x>=boxes[i-1].x+boxes[i-1].width,'targets do not overlap');
  }
 }
 await page.screenshot({path:'/tmp/feng-action-dock-mobile.png'});
 await page.setViewportSize({width:1440,height:1000});
 await page.screenshot({path:'/tmp/feng-action-dock-desktop.png'});
 await trigger.click();
 assert.ok(await dialog.isVisible());
 assert.ok(!(await dock.isVisible()),'toolbar yields to the open panel');
 assert.equal(await dialog.locator('[role=tablist],.support-unavailable,[data-support-image]').count(),0);
 assert.doesNotMatch(await dialog.innerText(),/WeChat|Alipay|微信|支付宝|暂未开通/);
 const paymentLink=dialog.locator('.support-paypal-link'),custom=dialog.locator('[data-support-custom]');
 assert.equal(await paymentLink.getAttribute('href'),'https://paypal.me/example/3USD');
 assert.doesNotMatch(await paymentLink.getAttribute('title'),/new window|新窗口/i);
 assert.ok(await custom.isVisible(),'custom amount is editable immediately');
 assert.ok(await custom.isEnabled());
 assert.equal(await custom.inputValue(),'10');
 const choose=async value=>dialog.locator(`.support-amounts label:has(input[value="${value}"])`).click();
 for(const amount of [1,3,5]){
  await choose(amount);
  assert.equal(await paymentLink.getAttribute('href'),`https://paypal.me/example/${amount}USD`);
  assert.ok(await custom.isVisible(),'custom amount stays visible with a preset selected');
 }
 await custom.click();
 assert.equal(await custom.inputValue(),'10');
 assert.equal(await paymentLink.getAttribute('href'),'https://paypal.me/example/10USD');
 assert.deepEqual(await custom.evaluate(n=>[n.selectionStart,n.selectionEnd]),[0,2],'initial 10 is selected for replacement');
 await page.keyboard.type('20');
 assert.equal(await custom.inputValue(),'20','typing replaces the initial amount instead of appending to it');
 assert.equal(await paymentLink.getAttribute('href'),'https://paypal.me/example/20USD');
 await custom.fill('12.50');
 assert.equal(await paymentLink.getAttribute('href'),'https://paypal.me/example/12.5USD');
 for(const invalid of ['','0','-1','1.001','abc','1e2']){
  await custom.fill(invalid);
  assert.equal(await paymentLink.getAttribute('href'),null,'invalid amount cannot link to checkout');
  assert.equal(await paymentLink.getAttribute('aria-disabled'),'true');
 }
 await paymentLink.click();
 assert.ok(await custom.evaluate(n=>n===document.activeElement),'invalid amount returns focus to input');
 await choose(1);
 assert.equal(await paymentLink.getAttribute('href'),'https://paypal.me/example/1USD','preset recovers from invalid custom input');
 await custom.fill('0.01');
 assert.equal(await dialog.locator('[name="support-amount"]:checked').count(),0,'typing directly deselects the preset');
 assert.equal(await paymentLink.getAttribute('href'),'https://paypal.me/example/0.01USD');
 await custom.fill('10');
 await choose(3);
 await dialog.locator('input[value="3"]').focus();
 await page.keyboard.press('ArrowRight');
 assert.equal(await paymentLink.getAttribute('href'),'https://paypal.me/example/5USD','keyboard changes selected amount');
 await page.keyboard.press('Tab');
 assert.ok(await custom.evaluate(n=>n===document.activeElement),'custom amount is in keyboard tab order');
 assert.equal(await paymentLink.getAttribute('href'),'https://paypal.me/example/5USD','tabbing past the custom input must not change the payment amount');
 assert.ok(!requests.some(u=>u.includes('paypal.me')),'opening the panel must not contact PayPal');
 for(let i=0;i<8;i++){
  await page.keyboard.press('Tab');
  assert.ok(await dialog.evaluate(n=>n.contains(document.activeElement)),'native modal keeps keyboard focus inside');
 }
 await page.keyboard.press('Escape');
 assert.ok(!(await dialog.isVisible()));
 assert.ok(await dock.isVisible());
 assert.ok(await trigger.evaluate(n=>n===document.activeElement));
 await page.addScriptTag({path:require.resolve('axe-core/axe.min.js')});
 for(const mode of ['en','zh','both']){
  await page.locator(`[data-language-choice=${mode}]`).click();
  await trigger.click();
  await custom.click();
   const violations=await page.evaluate(async()=>(await axe.run(document.getElementById('article-support-dialog'),{runOnly:{type:'tag',values:['wcag2a','wcag2aa','wcag21aa']}})).violations.map(v=>({id:v.id,nodes:v.nodes.map(n=>n.target)})));
   assert.deepEqual(violations,[],`${mode} accessibility`);
   for(const [width,height] of [[1440,1000],[768,800],[390,844],[320,568],[667,320]]){
    await page.setViewportSize({width,height});
    const rect=await dialog.boundingBox();
    assert.ok(rect.x>=0&&rect.x+rect.width<=width+1&&rect.y>=0&&rect.y+rect.height<=height+1,`${mode} ${width}: fits viewport`);
    assert.ok(await dialog.evaluate(n=>n.scrollWidth<=n.clientWidth+1),'no horizontal overflow');
   }
  await page.setViewportSize({width:1440,height:1000});
  await page.keyboard.press('Escape');
 }
 await page.locator('[data-support-open]').last().click();
 await page.screenshot({path:'/tmp/feng-support-desktop.png'});
 await page.setViewportSize({width:390,height:844});
 await page.screenshot({path:'/tmp/feng-support-mobile.png'});
 await page.locator('[data-support-close]').click();
 assert.ok(await page.locator('[data-support-open]').last().evaluate(n=>n===document.activeElement));
 await page.locator('[data-support-open]').last().click();
 await page.mouse.click(2,2);
 assert.ok(!(await dialog.isVisible()),'clicking the backdrop closes the panel');
 await page.locator('.discussion-launcher').click();
 assert.match(await page.locator('.discussion-shell').getAttribute('class'),/is-floating/);
 assert.ok(!(await dock.isVisible()),'toolbar does not cover the floating discussion');
 await page.locator('.discussion-close').click();
 assert.ok(await dock.isVisible());
 await page.locator('#article-discussions').scrollIntoViewIfNeeded();
 await page.waitForFunction(()=>document.querySelector('.discussion-launcher').hidden);
 assert.ok(await trigger.isVisible(),'support remains available at the end of the article');
 assert.ok(await dock.locator('[data-share-open]').isVisible(),'sharing remains available at the end of the article');
 assert.deepEqual(errors,[]);
 const nojs=await browser.newPage({javaScriptEnabled:false});
 await nojs.route(url,r=>r.fulfill({contentType:'text/html',body:fixtureHTML}));
 await nojs.goto(url);
 assert.equal(await nojs.locator('[data-support-open]:visible').count(),0);
 assert.ok(!(await nojs.locator('.article-action-dock').isVisible()),'no empty toolbar without JavaScript');
 assert.equal(await nojs.locator('.support-fallback a').getAttribute('href'),'https://paypal.me/example');
 assert.ok(await nojs.locator('.support-fallback a').isVisible());
 console.log('Support: payment targets, keyboard/focus, languages, responsive layout, accessibility and no-JS fallback passed.');
}finally{await browser.close();server.close();}
