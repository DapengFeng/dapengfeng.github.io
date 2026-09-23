import {chromium} from '@playwright/test';
import {createRequire} from 'node:module';
import fs from 'node:fs/promises';
import assert from 'node:assert/strict';
import {serve} from '../scripts/serve.mjs';
import {chromiumExecutable} from './browser-options.mjs';
const require=createRequire(import.meta.url),server=serve(4198);
const browser=await chromium.launch({executablePath:chromiumExecutable(),headless:true,args:['--no-sandbox']});
const failures=[];
await fs.rm('test-results/accessibility.json',{force:true});
try{
 const page=await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:'reduce'});
 await page.addInitScript(()=>localStorage.setItem('feng-language','en'));
 const posts=(await fs.readdir('content/posts')).filter(f=>f.endsWith('.html'));
 const series=(await fs.readdir('dist/series').catch(()=>[])).map(id=>'/series/'+id+'/');
 const urls=['/','/blog/','/categories/','/archive/','/about/','/404.html',...series,...posts.map(f=>'/blog/'+f)];
 async function audit(url,mode){
  const violations=await page.evaluate(async()=>{
   const result=await axe.run(document,{runOnly:{type:'tag',values:['wcag2a','wcag2aa','wcag21aa']}});
   return result.violations.map(v=>({id:v.id,impact:v.impact,help:v.help,nodes:v.nodes.map(n=>({target:n.target,failure:n.failureSummary}))}));
  });
  if(violations.length)failures.push({url,mode,violations});
 }
 for(const url of urls){
  await page.goto('http://localhost:4198'+url);
  await page.evaluate(()=>document.querySelectorAll('details').forEach(n=>n.open=true));
  await page.addScriptTag({path:require.resolve('axe-core/axe.min.js')});
  for(const mode of ['en','zh','both']){
   await page.locator(`[data-language-choice=${mode}]`).click();
   await audit(url,mode);
  }
  if(url.includes('pytorch-01')){
   await page.locator('[data-pt-phase=backward]').click();
   await page.locator('#pt-acc-a').click();
   for(const mode of ['en','zh','both']){await page.locator(`[data-language-choice=${mode}]`).click();await audit(url,'backward-'+mode);}
   await page.locator('[data-pt-grad=off]').click();await audit(url,'no-grad-backward');
  }
  if(url.includes('rust-vs-cpp'))for(const button of await page.locator('[data-code-mode]').all()){
   await button.click();await audit(url,'code-'+await button.getAttribute('data-code-mode'));
  }
  if(url.includes('spike_notes'))for(const button of await page.locator('[data-ultra-mode]').all()){
   await button.click();await audit(url,'ultra-'+await button.getAttribute('data-ultra-mode'));
  }
 }
 if(failures.length){await fs.mkdir('test-results',{recursive:true});await fs.writeFile('test-results/accessibility.json',JSON.stringify(failures,null,2));}
 assert.deepEqual(failures,[],'Accessibility regressions; see test-results/accessibility.json');
 console.log(`Accessibility checks passed: ${urls.length} pages × 3 languages, expanded content, and every UltraLIF tab.`);
}finally{await browser.close();server.close();}
