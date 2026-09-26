import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {load} from 'cheerio';
import {collectSeries,renderSeriesPreviews} from '../scripts/series.mjs';
import {seriesNavigation} from '../scripts/templates.mjs';

const first={slug:'first',url:'/blog/first.html',series:{id:'example-series',title:'示例专题',titleEn:'Example series',part:1}};
const third={slug:'third',url:'/blog/third.html',series:{...first.series,part:3}};
test('series order and neighboring links follow published installment numbers',()=>{
 const groups=collectSeries([third,{slug:'unrelated'},first]);
 assert.deepEqual(groups[0].posts.map(p=>p.slug),['first','third']);
 const nav=load(seriesNavigation(first,[third,first]));
 assert.equal(nav('a[rel="next"]').attr('href'),third.url);
 assert.equal(nav('a[rel="prev"]').length,0);
 assert.equal(load(seriesNavigation(third,[third,first]))('a[rel="prev"]').attr('href'),first.url);
 assert.throws(()=>collectSeries([first,{...third,series:first.series}]),/Duplicate series part/);
 assert.throws(()=>collectSeries([{...first,series:{...first.series,id:'../bad'}}]),/series requires/);
 assert.throws(()=>collectSeries([{...first,series:{...first.series,part:0}}]),/positive integer/);
 assert.throws(()=>collectSeries([first,{...third,series:{...third.series,title:'不一致'}}]),/Inconsistent/);
});
test('PyTorch series is discoverable and contains only published episodes',async()=>{
 const $=load(await fs.readFile('dist/series/pytorch-internals/index.html','utf8'));
 const published=[];
 for(const file of await fs.readdir('content/posts')){
  if(!file.endsWith('.html'))continue;
  const doc=load(await fs.readFile('content/posts/'+file,'utf8')),raw=doc('#article-metadata').text();
  if(!raw)continue;
  const metadata=JSON.parse(raw);
  if(!metadata.draft&&metadata.series?.id==='pytorch-internals')published.push({part:metadata.series.part,url:'/blog/'+file});
 }
 published.sort((a,b)=>a.part-b.part);
 assert.deepEqual($('.series-episodes a').map((_,el)=>$(el).attr('href')).get(),published.map(p=>p.url));
 assert.equal($('.series-episodes li').eq(1).attr('value'),'2');
 assert.equal($('.series-episodes a').eq(1).attr('href'),'/blog/pytorch-02-tensor-strides-storage.html');
 assert.equal($('.series-episodes li').attr('value'),'1');
 assert.equal($('.series-episodes a').attr('href'),'/blog/pytorch-01-what-is-pytorch.html');
 const second=load(await fs.readFile('dist/blog/pytorch-02-tensor-strides-storage.html','utf8'));
 assert.equal(second('a[rel=prev]').attr('href'),'/blog/pytorch-01-what-is-pytorch.html');
 const firstPage=load(await fs.readFile('dist/blog/pytorch-01-what-is-pytorch.html','utf8'));
 assert.equal(firstPage('a[rel=next]').attr('href'),'/blog/pytorch-02-tensor-strides-storage.html');
 assert.equal(firstPage('.series-preview [data-series-part=2]').attr('href'),'/blog/pytorch-02-tensor-strides-storage.html');
 assert.equal(firstPage('.series-preview [data-series-part=2] time').attr('datetime'),'2026-09-24');
 const thirdPage=load(await fs.readFile('dist/blog/pytorch-03-operator-dispatch.html','utf8'));
 assert.equal(firstPage('.series-preview [data-series-part=3]').attr('href'),'/blog/pytorch-03-operator-dispatch.html');
 assert.equal(second('a[rel=next]').attr('href'),'/blog/pytorch-03-operator-dispatch.html');
 assert.equal(thirdPage('a[rel=prev]').attr('href'),'/blog/pytorch-02-tensor-strides-storage.html');
 const afterThird=published.find(p=>p.part>3);
 assert.equal(thirdPage('a[rel=next]').attr('href'),afterThird?.url);
 const fourthPage=load(await fs.readFile('dist/blog/pytorch-04-autograd-engine.html','utf8'));
 assert.equal(thirdPage('a[rel=next]').attr('href'),'/blog/pytorch-04-autograd-engine.html');
 assert.equal(fourthPage('a[rel=prev]').attr('href'),'/blog/pytorch-03-operator-dispatch.html');
 assert.equal(fourthPage('.series-directory [data-series-part=4] time').attr('datetime'),'2026-09-26');
 const fifthPage=load(await fs.readFile('dist/blog/pytorch-05-cuda-streams-timing.html','utf8'));
 assert.equal(fourthPage('a[rel=next]').attr('href'),'/blog/pytorch-05-cuda-streams-timing.html');
 assert.equal(fifthPage('a[rel=prev]').attr('href'),'/blog/pytorch-04-autograd-engine.html');
 assert.equal(firstPage('.series-directory [data-series-part=5]').attr('href'),'/blog/pytorch-05-cuda-streams-timing.html');
 assert.equal(fifthPage('.series-directory [data-series-part=5] time').attr('datetime'),'2026-09-26');
 for(const [page,part]of [[firstPage,1],[second,2],[thirdPage,3],[fourthPage,4],[fifthPage,5]]){
  assert.equal(page('.series-directory').length,1);
  assert.equal(page('.series-directory [data-series-part]').length,6);
  assert.equal(page('.series-directory [aria-current=page]').attr('data-series-part'),String(part));
  assert.equal(page('#article-content .series-preview').length,0,'no duplicated contents in article body');
  assert.equal(page('.series-navigation .series-directory').length,0,'keep full directory out of the reading entrance');
  assert.equal(page('.reading-layout').next().hasClass('series-footer'),true,'directory follows the complete article');
  assert.equal(page('.series-footer').next().hasClass('related-section'),true,'series comes before broader related reading');
 }
 assert.equal(firstPage('.series-directory').text(),second('.series-directory').text());
 assert.ok((await fs.readFile('dist/sitemap.xml','utf8')).includes('/series/pytorch-internals/'));
 for(const file of ['index.html','blog/index.html','blog/pytorch-01-what-is-pytorch.html']){
  const page=load(await fs.readFile('dist/'+file,'utf8'));
  assert.ok(page('a[href="/series/pytorch-internals/"]').length,file);
 }
});

test('a roadmap declared once is inherited by other installments and replaced on publication',()=>{
 const roadmap=[{part:3,titleEn:'Planned third part',title:'第三期规划'}];
 const owner={...first,series:{...first.series,roadmap}};
 const reader={slug:'second',series:{...first.series,part:2},html:'<div data-series-preview></div>'};
 let page=load(renderSeriesPreviews(reader,[owner,reader]));
 assert.equal(page('[data-series-part=3]').attr('data-state'),'planned');
 assert.equal(page('[data-series-part=3] [data-lang=en]').first().text(),'Planned third part');
 const published={...third,titleEn:'Actual title',title:'实际标题',date:'2026-09-24'};
 page=load(renderSeriesPreviews(reader,[owner,reader,published]));
 assert.equal(page('[data-series-part=3]').attr('href'),third.url);
 assert.equal(page('[data-series-part=3] [data-lang=en]').first().text(),'Actual title');
 assert.throws(()=>collectSeries([{...owner,series:{...owner.series,roadmap:[{part:0}]}}]),/roadmap requires/);
});

test('roadmap links, titles and dates follow publication and removal automatically',()=>{
 const post={...first,titleEn:'First',title:'第一期',date:'2026-09-23',html:'<div data-series-preview><div data-series-part="3"><span data-lang="en">Planned topic</span><span data-lang="zh">规划主题</span></div></div><script>const n = 1 < 2;</script>'};
 const published={...third,titleEn:'Published & revised',title:'实际标题',date:'2026-09-24'};
 let page=load(renderSeriesPreviews(post,[post]));
 assert.equal(page('[data-series-part=3]').is('div'),true);
 assert.equal(page('[data-series-part=3] [data-lang=en]').first().text(),'Planned topic');
 page=load(renderSeriesPreviews(post,[post,published]));
 assert.equal(page('[data-series-part=3]').attr('href'),third.url);
 assert.equal(page('[data-series-part=3] [data-lang=en]').first().text(),'Published & revised');
 assert.equal(page('[data-series-part=3] time').text(),'2026-09-24');
 assert.equal(page('[data-series-part=1]').attr('aria-current'),'page');
 assert.equal(page('script').text(),'const n = 1 < 2;');
 const later={...published,slug:'later',url:'/blog/later.html',series:{...published.series,part:8}};
 page=load(renderSeriesPreviews(post,[later,published,post]));
 assert.deepEqual(page('[data-series-part]').map((_,el)=>Number(page(el).attr('data-series-part'))).get(),[1,3,8]);
 page=load(renderSeriesPreviews(post,[post]));
 assert.equal(page('[data-series-part=3]').attr('href'),undefined);
 assert.equal(page('[data-series-part=3] time').length,0);
 assert.throws(()=>renderSeriesPreviews({...post,html:'<div data-series-preview="missing"></div>'},[post]),/unknown preview series/);
});
