import {chromiumExecutable} from './browser-options.mjs';
import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import {serve} from '../scripts/serve.mjs';
const server=serve(4182),browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox','--no-proxy-server']});

async function scenario({locale='en-US',saved,automatic,blocked=false}={}){
 const context=await browser.newContext({locale});let locationRequests=0;
 await context.addInitScript(({saved,automatic,blocked})=>{
  if(blocked){for(const key of ['localStorage','sessionStorage'])Object.defineProperty(window,key,{get(){throw Error('Storage blocked');}});}
  else{if(saved)localStorage.setItem('feng-language',saved);if(automatic)sessionStorage.setItem('feng-auto-language',automatic);}
 },{saved,automatic,blocked});
 await context.route('https://api.country.is/**',async route=>{locationRequests++;await route.fulfill({json:{country:'CN'}});});
 const page=await context.newPage();await page.goto('http://localhost:4182/');
 return {context,page,locationRequests:()=>locationRequests};
}

try{
 for(const [locale,expected] of [['en-US','en'],['en-GB','en'],['zh-CN','zh'],['zh-TW','zh'],['fr-FR','en']]){
  // A session value left by an old site version must not override the browser.
  const s=await scenario({locale,automatic:expected==='en'?'zh':'en'});
  assert.equal(await s.page.locator('html').getAttribute('data-language'),expected);
  assert.equal(await s.page.locator('html').getAttribute('data-language-source'),'browser');
  assert.equal(await s.page.evaluate(()=>localStorage.getItem('feng-language')),null);
  await s.page.reload();assert.equal(await s.page.locator('html').getAttribute('data-language'),expected);
  assert.equal(s.locationRequests(),0);await s.context.close();
 }
 for(const saved of ['en','zh','both']){
  const s=await scenario({saved,locale:saved==='zh'?'en-US':'zh-CN'});
  assert.equal(await s.page.locator('html').getAttribute('data-language'),saved);
  assert.equal(await s.page.locator('html').getAttribute('data-language-source'),'manual');
  assert.equal(s.locationRequests(),0);await s.context.close();
 }
 for(const options of [{saved:'invalid',locale:'zh-CN'},{blocked:true,locale:'zh-CN'}]){
  const s=await scenario(options);assert.equal(await s.page.locator('html').getAttribute('data-language'),'zh');
  await s.page.locator('[data-language-choice=both]').click();assert.equal(await s.page.locator('html').getAttribute('data-language'),'both');
  assert.equal(s.locationRequests(),0);await s.context.close();
 }
 const s=await scenario();
 await s.page.locator('[data-language-choice=both]').click();await s.page.reload();
 assert.equal(await s.page.locator('html').getAttribute('data-language'),'both');
 assert.equal(await s.page.evaluate(()=>sessionStorage.getItem('feng-auto-language')),null);
 await s.page.goto('http://localhost:4182/blog/waves-and-phase.html');
 const reading=s.page.locator('.article-byline [data-reading-minutes]').first();
 assert.equal(await reading.count(),1,'article metadata carries per-language reading estimates');
 for(const language of ['en','zh','both']){
  await s.page.locator(`[data-language-choice=${language}]`).click();
  const expected=await reading.getAttribute(`data-minutes-${language}`);
  assert.ok(Number(expected)>=1);
  assert.match(await reading.innerText(),new RegExp(`^${expected} ${language==='zh'?'分钟阅读':'min read'}`));
 }
 // Preserve a running experiment when changing the reading language.
 await s.page.locator('[data-wave-lab] > .reading-explore > summary').click();
 await s.page.locator('[data-wave-phase]').fill('2.5');await s.page.locator('[data-language-choice=zh]').click();
 assert.equal(await s.page.locator('[data-wave-phase]').inputValue(),'2.5');
 assert.equal(s.locationRequests(),0);await s.context.close();

 // A mixed index ranks titles above tags and prose and uses a language-specific
 // excerpt around the match, rather than concatenating two article editions.
 const search=await scenario({saved:'en'});
 const records=[
  {kind:'article',url:'/blog/example.html',title:'信号笔记',titleEn:'Signal notes',description:'一篇信号笔记',descriptionEn:'A short signal note',searchTextEn:'Earlier context. '.repeat(30)+'sampling preserves the signal. <img src=x onerror=alert(1)>',searchTextZh:'前面的解释。'.repeat(30)+'采样保留了信号。'},
  {kind:'series',url:'/series/example/',title:'信号处理系列',titleEn:'Signals series',description:'从采样开始',descriptionEn:'A series about signals',tags:['sampling','采样']},
  {kind:'math',url:'/math/example/',title:'采样定理',titleEn:'Sampling theorem',description:'样本之间的关系',descriptionEn:'How samples relate',tags:[]},
 ];
 await search.context.route('**/search-index.json',route=>route.fulfill({json:records}));
 await search.page.locator('.search-trigger').first().click();
 await search.page.locator('#global-search').fill('sampling');
 await search.page.waitForFunction(()=>document.querySelector('.search-result')?.getAttribute('href')==='/math/example/');
 assert.deepEqual(await search.page.locator('.search-result').evaluateAll(nodes=>nodes.map(node=>node.getAttribute('href'))),['/math/example/','/series/example/','/blog/example.html']);
 assert.match(await search.page.locator('.search-result').last().locator(':scope > span [data-lang=en]').innerText(),/sampling preserves the signal/);
 assert.equal(await search.page.locator('.search-result img').count(),0,'excerpts are escaped');
 await search.page.locator('[data-close-search]').click();
 await search.page.locator('[data-language-choice=zh]').click();
 await search.page.locator('.search-trigger').first().click();
 await search.page.locator('#global-search').fill('采样');
 await search.page.waitForFunction(()=>document.querySelector('.search-result:last-child > span [data-lang=zh]')?.textContent.includes('采样保留'));
 const excerpt=await search.page.locator('.search-result').last().locator(':scope > span').innerText();
 assert.match(excerpt,/采样保留了信号/);assert.doesNotMatch(excerpt,/Earlier context|sampling preserves/);
 await search.context.close();
 console.log('Language/search checks passed: browser locale, persisted choice, blocked storage, no IP lookup, reading estimates, native language, experiment state, mixed search ranking and localized excerpts.');
}finally{await browser.close();server.close();}
