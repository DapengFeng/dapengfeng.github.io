import {chromiumExecutable} from './browser-options.mjs';
import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import {serve} from '../scripts/serve.mjs';
const server=serve(4182),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox','--no-proxy-server']});
const endpoint='https://api.country.is/';
async function scenario({locale='en-US',country='CN',saved,fail=false,invalid=false,slow=false,blocked=false}={}){
 const context=await browser.newContext({locale});let requests=0;
 await context.addInitScript(({saved,blocked})=>{if(blocked){for(const key of ['localStorage','sessionStorage'])Object.defineProperty(window,key,{get(){throw Error('Storage blocked');}});}else if(saved)localStorage.setItem('feng-language',saved);},{saved,blocked});
 let release;
 await context.route(endpoint,async route=>{requests++;if(slow)await new Promise(resolve=>release=resolve);try{if(fail)await route.abort();else await route.fulfill({json:invalid?{country:null}:{country}});}catch{}});
 const page=await context.newPage();await page.goto('http://localhost:4182/');
 return {context,page,requests:()=>requests,release:()=>release?.()};
}
try{
 for(const country of ['CN','HK','MO','TW','US','DE']){
  const s=await scenario({country,locale:country==='US'?'zh-CN':'en-US'});
  await s.page.waitForFunction(()=>document.documentElement.dataset.languageSource==='auto');
  const expected=['CN','HK','MO','TW'].includes(country)?'zh':'en';
  assert.equal(await s.page.locator('html').getAttribute('data-language'),expected);
  assert.equal(await s.page.evaluate(()=>localStorage.getItem('feng-language')),null);
  await s.page.reload();assert.equal(s.requests(),1);assert.equal(await s.page.locator('html').getAttribute('data-language'),expected);
  await s.context.close();
 }
 for(const saved of ['en','zh','both']){
  const s=await scenario({saved});assert.equal(s.requests(),0);assert.equal(await s.page.locator('html').getAttribute('data-language'),saved);await s.context.close();
 }
 for(const options of [{fail:true,locale:'zh-CN'},{invalid:true},{blocked:true,country:'CN'}]){
  const s=await scenario(options);await s.page.waitForFunction(()=>document.documentElement.dataset.languageSource==='auto');assert.equal(await s.page.locator('html').getAttribute('data-language'),options.locale==='zh-CN'||options.blocked?'zh':'en');await s.context.close();
 }
 const race=await scenario({slow:true});await race.page.locator('[data-language-choice=both]').click();race.release();await race.page.reload();assert.equal(await race.page.locator('html').getAttribute('data-language'),'both');assert.equal(race.requests(),1);await race.context.close();
 const timeout=await scenario({slow:true,locale:'zh-CN'});await timeout.page.waitForFunction(()=>document.documentElement.dataset.languageSource==='auto');assert.equal(await timeout.page.locator('html').getAttribute('data-language'),'zh');timeout.release();await timeout.context.close();
 // Exercise a language switch while an article experiment is already running.
 const s=await scenario({saved:'both'});await s.page.goto('http://localhost:4182/blog/waves-and-phase.html');await s.page.locator('[data-wave-phase]').fill('2.5');await s.page.locator('[data-language-choice=zh]').click();assert.equal(await s.page.locator('[data-wave-phase]').inputValue(),'2.5');await s.context.close();
 console.log('Language checks passed: country mapping, saved preferences, session cache, browser fallback, timeout, blocked storage, manual-choice race, experiment state.');
}finally{await browser.close();server.close();}
