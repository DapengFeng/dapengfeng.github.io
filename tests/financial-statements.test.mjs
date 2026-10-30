import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {load} from 'cheerio';
import {statements} from '../content/assets/financial-statements/model.js';
import {collectSeries} from '../scripts/series.mjs';
import {seriesNavigation} from '../scripts/templates.mjs';

test('the fictional shop reconciles its income, cash ledger and closing balance sheet',()=>{
 const v=statements();
 assert.equal(v.profit,30);
 assert.equal(v.operatingCash,11);
 assert.equal(v.cash,100+40-60-35+80-20-4-10);
 assert.equal(v.cash,91);
 assert.equal(v.assets,185);
 assert.equal(v.assets,v.liabilitiesAndEquity);
 assert.equal(v.equity,130);
 assert.equal(v.operatingCash,v.indirectCash);
 assert.equal(v.investingCash,-60);
 assert.equal(v.financingCash,140);
});
test('collection changes cash and receivables without inventing profit or assets',()=>{
 for(let amount=0;amount<=100;amount++){
  const v=statements(amount);
  assert.equal(v.profit,30);
  assert.equal(v.assets,185);
  assert.equal(v.equity,130);
  assert.equal(v.receivables+v.collected,100);
  assert.equal(v.cash,amount+11);
  assert.equal(v.operatingCash,amount-69);
  assert.equal(v.indirectCash,v.operatingCash);
  assert.equal(v.cash,v.operatingCash+v.investingCash+v.financingCash);
  assert.equal(v.assets,v.liabilitiesAndEquity);
 }
 for(const amount of [-1,101,NaN,Infinity,'80'])assert.throws(()=>statements(amount),RangeError);
});
test('the no-JavaScript figure matches the reconciled default example',async()=>{
 const $=load(await fs.readFile('content/posts/financial-statements-01-three-statements.html','utf8'));
 const figure=$('[data-collection-example]'),v=statements(Number(figure.find('input').attr('value')));
 assert.equal(figure.length,1);
 for(const el of figure.find('[data-account]').toArray())assert.equal(Number($(el).text()),v[$(el).attr('data-account')]);
 assert.equal(figure.find('input').attr('disabled'),'disabled');
});

test('the learning path and filenames follow the seven-part order',async()=>{
 const order=[
  'financial-statements-01-three-statements',
  'financial-statements-02-owner-earnings',
  'financial-statements-03-returns-and-capital-allocation',
  'financial-statements-04-management-incentives',
  'financial-statements-05-liquidity-risk',
  'financial-statements-06-value-and-price',
  'financial-statements-07-annual-report-practice'
 ];
 const posts=[];
 for(const file of (await fs.readdir('content/posts')).filter(f=>f.startsWith('financial-statements-'))){
  const $=load(await fs.readFile('content/posts/'+file,'utf8'));
  const meta=JSON.parse($('#article-metadata').text());
  assert.equal($('.fs-article').attr('data-fs-part'),String(meta.series.part),file);
  assert.equal(Number(/^financial-statements-(\d+)-/.exec(file)[1]),meta.series.part,file);
  posts.push({...meta,slug:file.slice(0,-5),url:'/blog/'+file});
 }
 const group=collectSeries(posts)[0];
 assert.deepEqual(group.posts.map(p=>p.slug),order);
 for(const [i,post] of group.posts.entries()){
  const nav=load(seriesNavigation(post,posts));
  assert.equal(nav('a[rel=prev]').attr('href'),group.posts[i-1]?.url);
  assert.equal(nav('a[rel=next]').attr('href'),group.posts[i+1]?.url);
 }
 const guide=JSON.parse(await fs.readFile('content/editorial.json','utf8')).series.find(s=>s.id===group.id);
 assert.deepEqual(guide.paths.flatMap(p=>p.parts),[1,2,3,4,5,6,7]);
});
