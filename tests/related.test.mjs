import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {load} from 'cheerio';
import {createRelatedRecommender} from '../scripts/related.mjs';
import {loadContent} from '../scripts/content.mjs';
const post=(slug,fields={})=>({slug,url:`/blog/${slug}.html`,titleEn:'',title:'',descriptionEn:'',description:'',tags:[],headings:[],category:'math',...fields});
test('content relevance beats category and publication order in both languages',()=>{
 for(const text of ['matrix multiplication blocked cache locality','矩阵乘法 分块计算 缓存局部性']){
  const a=post('source',{titleEn:text,relatedText:text}),close=post('close',{category:'systems',titleEn:text,relatedText:text}),far=post('far',{titleEn:'ocean tides waves phase',date:'2099-01-01'});
  const recommendations=createRelatedRecommender([far,a,close])(a);
  assert.equal(recommendations[0].slug,'close');assert.ok(!recommendations.some(p=>p.slug==='far'));
 }
});
test('body-only overlap is useful; unavailable articles and self are excluded',()=>{
 const a=post('a',{relatedText:'tensors strides storage transposition'}),b=post('b',{relatedText:'tensors strides storage transposition'});
 const draft=post('draft',{...b,slug:'draft',draft:true}),archived=post('archive',{...b,slug:'archive',archiveOnly:true});
 assert.deepEqual(createRelatedRecommender([a,b,draft,archived])(a).map(p=>p.slug),['b']);
 assert.deepEqual(createRelatedRecommender([a])(a),[]);
});
test('same-series recommendations favor neighboring installments with stable ties',()=>{
 const series={id:'tensor',part:1},a=post('a',{series}),near=post('near',{series:{...series,part:2}}),distant=post('distant',{series:{...series,part:8}});
 assert.deepEqual(createRelatedRecommender([distant,a,near])(a).map(p=>p.slug),['near','distant']);
 const b=post('b',{tags:['FFT']}),c=post('c',{tags:['fft']}),d=post('d',{tags:['fft']});
 assert.deepEqual(createRelatedRecommender([d,c,b])(b).map(p=>p.slug),createRelatedRecommender([b,c,d])(b).map(p=>p.slug));
 assert.equal(createRelatedRecommender([b,c,d])(b,1).length,1);
});
test('newly added content is considered on the next build',()=>{
 const a=post('a',{titleEn:'tensor strides storage'}),old=post('old',{titleEn:'tensor kernel operators'}),fresh=post('fresh',{titleEn:'tensor strides storage'});
 assert.notEqual(createRelatedRecommender([a,old])(a)[0]?.slug,'fresh');
 assert.equal(createRelatedRecommender([a,old,fresh])(a)[0].slug,'fresh');
});
test('published recommendation cards prioritize explained editorial connections and retain relevant fallbacks',async()=>{
 const posts=await loadContent(),recommend=createRelatedRecommender(posts);
 const editorial=JSON.parse(await fs.readFile('content/editorial.json','utf8'));
 for(const p of posts){
  const $=load(await fs.readFile('dist/'+p.url.slice(1),'utf8'));
  const links=$('.related-section .card-link').map((_,el)=>$(el).attr('href')).get();
  const curated=editorial.related.filter(link=>link.source===p.slug).slice(0,3);
  const curatedUrls=curated.map(link=>posts.find(post=>post.slug===link.target)?.url);
  assert.ok(curatedUrls.every(Boolean),`${p.slug}: editorial targets are published`);
  assert.deepEqual(links.slice(0,curated.length),curatedUrls,`${p.slug}: editorial connections come first`);
  assert.equal(new Set(links).size,links.length,`${p.slug}: no duplicate reading suggestions`);
  assert.ok(!links.includes(p.url)&&links.length<=3,`${p.slug}: no self-link or excessive suggestions`);
  const remaining=recommend(p).filter(post=>!curatedUrls.includes(post.url)).slice(0,3-curated.length);
  assert.deepEqual(links.slice(curated.length),remaining.map(post=>post.url),`${p.slug}: automatic relevance fills remaining places`);
  for(const [i,connection]of curated.entries()){
   const reason=$('.related-section .card-connection').eq(i);
   assert.equal(reason.find('[data-lang="en"]').text(),connection.reason.en,`${p.slug}: English connection reason`);
   assert.equal(reason.find('[data-lang="zh"]').text(),connection.reason.zh,`${p.slug}: Chinese connection reason`);
  }
 }
 const matrix=posts.find(p=>p.slug==='matrix-multiplication');
 assert.ok(recommend(matrix).some(p=>p.slug==='band-storage-gaxpy'));
 const pytorch=posts.find(p=>p.slug==='pytorch-01-what-is-pytorch');
 assert.equal(recommend(pytorch)[0].series?.id,pytorch.series.id);
 assert.ok(recommend(matrix).some(p=>p.category!==matrix.category));
});
