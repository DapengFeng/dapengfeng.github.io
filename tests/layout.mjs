import {chromiumExecutable} from './browser-options.mjs';
import {chromium} from '@playwright/test';
import fs from 'node:fs/promises';
import assert from 'node:assert/strict';
import {serve} from '../scripts/serve.mjs';
const before=process.argv.includes('--audit'),server=serve(4191);
const browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
const posts=(await fs.readdir('content/posts')).filter(f=>f.endsWith('.html'));
const urls=['/','/blog/','/categories/','/archive/','/about/','/404.html',...posts.map(p=>'/blog/'+p)];
try{
 const page=await browser.newPage(),errors=[],tight=[];page.on('pageerror',e=>errors.push(e.message));await page.addInitScript(()=>localStorage.setItem('feng-language','both'));
 for(const url of urls){
  await page.goto('http://localhost:4191'+url);
  for(const width of [1440,768,390,320]){
   await page.setViewportSize({width,height:1000});
   for(const language of ['both','en','zh']){
    await page.locator(`[data-language-choice="${language}"]`).click();
    await page.evaluate(()=>new Promise(requestAnimationFrame));
    const fits=await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1);
    if(before&&!fits)console.log('OVERFLOW',url,width,language,await page.evaluate(()=>[...document.querySelectorAll('body *')].filter(n=>n.getBoundingClientRect().right>innerWidth+1&&n.getBoundingClientRect().width>0&&!n.closest('pre,svg,.katex,.math-display,.math-inline')).slice(-20).map(n=>({tag:n.tagName,cls:n.className,id:n.id,right:n.getBoundingClientRect().right}))));
    else assert.ok(fits,`${url} overflow: ${width}/${language}`);
    if(language==='both'&&[1440,390].includes(width)){
     const gaps=await page.evaluate(()=>{
      const hits=[];
      document.querySelectorAll('.lesson,.article-section,#article>section,.card-body,.controls,.compiler-check,.card-copy,.prose').forEach(parent=>{
       const nodes=[...parent.children].filter(n=>n.getBoundingClientRect().height>0&&!['SCRIPT','STYLE','SPAN'].includes(n.tagName));
       for(let i=1;i<nodes.length;i++){
        const a=nodes[i-1],b=nodes[i],ar=a.getBoundingClientRect(),br=b.getBoundingClientRect(),gap=br.top-ar.bottom;
        if(br.top>=ar.top&&br.left<ar.right&&br.right>ar.left&&gap<12&&!a.matches('h2,h3,h4')&&!b.matches('h2,h3,h4')&&b.textContent.trim())hits.push({parent:parent.className||parent.id,a:a.tagName+'.'+a.className,b:b.tagName+'.'+b.className,gap:Math.round(gap),text:b.textContent.trim().slice(0,45)});
       }
      });return hits.slice(0,15);
     });
     if(gaps.length)tight.push({url,width,gaps});
    }
   }
  }
 }
 if(before)console.log(JSON.stringify(tight,null,2));else assert.deepEqual(tight,[],'content blocks need at least 12px clearance');assert.deepEqual(errors,[]);
 for(const [url,target,name]of [['/','.category-section','home'],['/blog/rust-vs-cpp-blog.html','#tradeoffs','rust'],['/blog/spike_notes.html','#article','spike'],['/blog/waves-and-phase.html','#topic-2','waves'],['/blog/band-storage-gaxpy.html','.lesson','band'],['/about/','.about-grid','about']]){
  await page.goto('http://localhost:4191'+url);await page.locator('[data-language-choice=both]').click();
  for(const width of [1440,390]){await page.setViewportSize({width,height:1000});await page.locator(target).evaluate(n=>window.scrollTo(0,Math.max(0,n.getBoundingClientRect().top+scrollY-28)));await page.screenshot({path:`/tmp/feng-layout-${before?'before':'after'}-${name}-${width}.png`});}
 }
 // Open long derivations and check the standalone maintenance guide too.
 for(const url of ['/blog/spike_notes.html','/blog/rust-vs-cpp-blog.html',new URL('../README.html',import.meta.url).href]){
  await page.goto(url.startsWith('file:')?url:'http://localhost:4191'+url);
  await page.evaluate(()=>document.querySelectorAll('details').forEach(n=>n.open=true));
  for(const width of [1440,768,390,320]){
   await page.setViewportSize({width,height:1000});
   for(const mode of ['both','en','zh']){
    await page.locator(`[data-language-choice=${mode}]`).click();
    await page.evaluate(()=>new Promise(requestAnimationFrame));
    assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`expanded content overflow: ${url}/${width}/${mode}`);
   }
  }
 }
 console.log(`Checked ${urls.length} pages × 4 widths × 3 languages; no overflow or script errors.`);
}finally{await browser.close();server.close();}
