import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import * as cheerio from 'cheerio';
import {loadContent,parseContent} from '../scripts/content.mjs';
const posts=await loadContent();
test('all articles contain dated bilingual content and valid category metadata',()=>{
 assert.ok(posts.length>=10);assert.equal(new Set(posts.map(p=>p.url)).size,posts.length);
 for(const p of posts){if(['spike_notes','rust-vs-cpp-blog','cuda-rust-two-tracks-blog'].includes(p.slug))assert.ok(p.hasTranslation,p.slug);assert.ok(p.html.includes('data-lang="en"'));assert.ok(p.html.includes('data-lang="zh"'));assert.ok(p.title,p.slug);}
 for(const slug of ['spike_notes','rust-vs-cpp-blog','cuda-rust-two-tracks-blog']){const p=posts.find(p=>p.slug===slug);assert.equal(p.date,'2026-09-23');assert.equal(p.dateKind,'published');}
});
test('new standalone HTML is indexed from metadata and its scripts and styles survive',async()=>{
 const dir=await fs.mkdtemp(path.join(os.tmpdir(),'feng-content-'));
 try{
  const file=path.join(dir,'example.html');
  await fs.writeFile(file,'<!doctype html><html><head><title>A new experiment</title><meta name="article:published_time" content="2026-09-23"><meta name="category" content="math"><style>body {color:red} h2 {color:blue}</style></head><body><h2 class="parallel-text"><span data-lang="en">Experiment</span><span data-lang="zh">实验</span></h2><h2 class="parallel-text"><span data-lang="en">Experiment</span><span data-lang="zh">实验</span></h2><div data-math="y=Ax"></div><button id="demo">Run</button><script>window.demo = true;</script></body></html>');
  const p=await parseContent(file);assert.equal(p.category,'math');assert.equal(p.headings.length,2);assert.equal(p.headings[0].titleEn,'Experiment');assert.equal(p.headings[0].titleZh,'实验');assert.notEqual(p.headings[0].id,p.headings[1].id);assert.match(p.styles,/\.legacy-content h2/);assert.match(p.html,/window.demo = true/);assert.match(p.html,/<mjx-container/);
  await fs.writeFile(file,'<title>Missing date</title><meta name="category" content="math"><p>Test</p>');await assert.rejects(()=>parseContent(file),/date are required/);
  await fs.writeFile(file,'<script type="application/json" id="article-metadata">{"draft":true}</script><p>Not published</p>');assert.equal(await parseContent(file),null);
 }finally{await fs.rm(dir,{recursive:true,force:true});}
});
test('generated pages have unique IDs, resolvable internal links and heading targets',async()=>{
 const files=await fs.readdir('dist',{recursive:true});let checked=0;
 for(const file of files.filter(f=>f.endsWith('.html'))){
  const $=cheerio.load(await fs.readFile(path.join('dist',file),'utf8')),ids=$('[id]').map((_,el)=>$(el).attr('id')).get();assert.equal(new Set(ids).size,ids.length,`duplicate IDs in ${file}`);
  for(const el of $('a[href],img[src],script[src],link[href]').toArray()){
   const value=$(el).attr('href')||$(el).attr('src');if(!value||!value.startsWith('/')||value.startsWith('//'))continue;
   const [pathname,hash]=value.split('#'),local=pathname.split('?')[0];let target=path.join('dist',local);if(local.endsWith('/'))target=path.join(target,'index.html');await assert.doesNotReject(()=>fs.access(target),`${file} links to missing ${value}`);checked++;
  }
  for(const el of $('.article-toc a[href^="#"]').toArray()){const id=$(el).attr('href').slice(1);assert.ok(ids.includes(id),`missing TOC target ${id}`);}
 }
 assert.ok(checked>100);
});
test('RSS, sitemap and search index contain every published note',async()=>{
 const search=JSON.parse(await fs.readFile('dist/search-index.json','utf8')),feed=await fs.readFile('dist/feed.xml','utf8'),sitemap=await fs.readFile('dist/sitemap.xml','utf8');
 assert.equal(search.length,posts.length);for(const p of posts){assert.ok(feed.includes(p.url));assert.ok(sitemap.includes(p.url));assert.ok(search.find(x=>x.slug===p.slug));}
});
test('English precedes its Chinese counterpart without separate full-article editions',()=>{
 const spike=posts.find(p=>p.slug==='spike_notes'),$=cheerio.load(spike.html);
 assert.ok($('.parallel-text').length>0,'article contains paired bilingual text');
 $('.parallel-text').each((_,node)=>{assert.equal($(node).children().eq(0).attr('data-lang'),'en');assert.equal($(node).children().eq(1).attr('data-lang'),'zh');});
 assert.equal(spike.headings.filter(h=>h.level===2).length,8);
 assert.ok(spike.headings.filter(h=>h.level===2).every(h=>h.titleEn&&h.titleZh));
 for(const slug of ['waves-and-phase','benchmark-with-evidence']){
  const doc=cheerio.load(posts.find(p=>p.slug===slug).html);
  assert.equal(doc('.bilingual-chapter').length,5);
  doc('.paragraph-pair').each((_,node)=>assert.deepEqual(doc(node).children().map((_,x)=>doc(x).attr('data-lang')).get(),['en','zh']));
 }
});

test('one bilingual HTML is sufficient for discovery, metadata, contents and search',async()=>{
 const file='content/posts/single-file-workflow-check.html';
 try {
  await fs.copyFile('examples/article.html',file);
  const discovered=(await loadContent()).find(p=>p.slug==='single-file-workflow-check');
  assert.ok(discovered);
  assert.equal(discovered.titleEn,'A visual note: vectors, formulas, and code');
  assert.equal(discovered.title,'可视化笔记：向量、公式与代码');
  assert.equal(discovered.url,'/blog/single-file-workflow-check.html');
  assert.ok(discovered.hasTranslation);
  assert.equal(discovered.headings[0].titleEn,'One article, two languages');
  assert.equal(discovered.headings[0].titleZh,'一篇文章，两种语言');
  assert.match(discovered.searchText,/Code with syntax colors/);
  assert.match(discovered.searchText,/带语法配色的代码/);
  assert.match(discovered.html,/<mjx-container/);
  const raw=await fs.readFile(file,'utf8');
  await fs.writeFile(file,raw.replaceAll('data-lang="zh"','data-unused="zh"'));
  await assert.rejects(()=>parseContent(file),/include both English and Chinese/);
 } finally { await fs.rm(file,{force:true}); }
});

test('full explanations are paired and preserve the same inline formulas',()=>{
 for(const p of posts){
  const $=cheerio.load(p.html);
  assert.doesNotMatch(p.html,/中文阅读版|English reading edition|original Chinese edition contains/);
  $('.parallel-text').each((_,node)=>{
   const en=$(node).children('[data-lang="en"]'),zh=$(node).children('[data-lang="zh"]');
   // A wrapper can contain a separately paired heading or shared controls.
   if(!en.length&&!zh.length)return;
   assert.equal(en.length,1,p.slug);assert.equal(zh.length,1,p.slug);
   assert.equal($(node).children('[data-lang]').first().attr('data-lang'),'en',p.slug);
   const equations=part=>part.find('[data-math],.math-inline').map((_,el)=>$(el).attr('data-math')||$(el).attr('aria-label')).get().sort();
   assert.deepEqual(equations(en),equations(zh),`${p.slug}: inline formulas differ`);
  });
 }
 const rust=cheerio.load(posts.find(p=>p.slug==='rust-vs-cpp-blog').html);
 assert.ok(rust('.parallel-text').length>175);
 assert.ok(rust('#article-terms').length);
 const matrix=cheerio.load(posts.find(p=>p.slug==='matrix-multiplication').html);
 assert.equal(matrix('h3.parallel-text').length,4);
});

test('every published display equation has a number',()=>{
 for(const post of posts){
  const $=cheerio.load(post.html);
  $('.formula-block').each((_,el)=>{
   assert.ok($(el).find('.eq-number').text().trim()||$(el).find('[data-mml-node="mlabeledtr"]').length,`${post.slug}: unnumbered display formula`);
  });
 }
});

test('Spike display numbers are contiguous within each numbered chapter',()=>{
 const spike=posts.find(p=>p.slug==='spike_notes'),$=cheerio.load(spike.html);
 const chapters=['overview-title','lif-title','ultralif-title','bptt-title','stdp-title','eventprop-title','comparison-title','references-title'];
 let section=0,sequence=0,total=0;
 $('h2,.formula-block').each((_,el)=>{
  const node=$(el);
  if(node.is('h2')){section=chapters.indexOf(node.attr('id'))+1;sequence=0;return;}
  assert.ok(section>0,'a display equation belongs to a numbered chapter');
  assert.equal(node.attr('data-equation-number'),`${section}.${++sequence}`);
  total++;
 });
 assert.equal(total,72);
});

test('Spike contents, body chapters, and equation prefixes share one-based numbering',async()=>{
 const $=cheerio.load(await fs.readFile('dist/blog/spike_notes.html','utf8'));
 const chapters=$('#article-content section.chapter');
 const toc=$('.article-toc nav a');
 assert.equal(chapters.length,8);
 assert.equal(toc.length,8);
 chapters.each((i,el)=>{
  const number=String(i+1).padStart(2,'0'),node=$(el);
  assert.equal(node.find('.chapter-no').first().text(),number);
  assert.equal(toc.eq(i).children('span').first().text(),number);
  assert.equal(toc.eq(i).attr('href'),'#'+node.find('h2').first().attr('id'));
  node.find('.formula-block[data-equation-number]').each((j,formula)=>assert.equal($(formula).attr('data-equation-number'),`${i+1}.${j+1}`));
 });
});

test('shared math glyphs resolve locally on every generated page',async()=>{
 for(const file of (await fs.readdir('dist',{recursive:true})).filter(f=>f.endsWith('.html'))){
  const $=cheerio.load(await fs.readFile('dist/'+file,'utf8'));
  const idCounts=new Map();
  $('[id]').each((_,el)=>{const id=$(el).attr('id');idCounts.set(id,(idCounts.get(id)||0)+1);});
  $('mjx-container use').each((_,el)=>{
   const href=$(el).attr('href');
   assert.ok(href?.startsWith('#'),`${file}: glyph must stay local`);
   assert.equal(idCounts.get(href.slice(1)),1,`${file}: missing or duplicate glyph ${href}`);
  });
  assert.equal($('mjx-container svg:not([aria-hidden="true"])').length,0,`${file}: decorative glyphs should not repeat accessible math labels`);
  $('mjx-container').each((_,el)=>assert.ok($(el).attr('aria-label'),`${file}: formula needs an accessible name`));
 }
});
