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

test('the learning path connects the foundations and insurer cases in order',async()=>{
 const order=[
  'financial-statements-01-three-statements',
  'financial-statements-02-owner-earnings',
  'financial-statements-03-returns-and-capital-allocation',
  'financial-statements-04-management-incentives',
  'financial-statements-05-liquidity-risk',
  'financial-statements-06-value-and-price',
  'financial-statements-07-annual-report-practice',
  'financial-statements-08-cpic-profit-quality',
  'financial-statements-09-nci-new-business-value',
  'financial-statements-10-pingan-capital-allocation',
  'financial-statements-11-china-life-asset-liability',
  'financial-statements-12-picc-underwriting-discipline',
  'financial-statements-13-insurer-comparison'
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
 assert.deepEqual(guide.paths.flatMap(p=>p.parts),[1,2,3,4,5,6,7,8,9,10,11,12,13]);
});

test('insurer comparisons retain equal reporting periods, units and traceable source values',async()=>{
 const data=JSON.parse(await fs.readFile('content/assets/financial-statements/insurance-2024-2026.json','utf8'));
 const $=load(await fs.readFile('content/posts/financial-statements-13-insurer-comparison.html','utf8'));
 assert.equal(data.unit,'RMB million');
 assert.equal(data.periods.m92026.end,'2026-09-30');
 assert.equal(Object.hasOwn(data.periods,'fy2026'),false,'partial-year results must not become annual actuals');
 assert.equal(data.companies.length,5);
 assert.equal($('tr[data-insurer]').length,10);
 for(const company of data.companies){
  for(const key of ['annualSource','quarterSource']){
   assert.match(company[key].url,/^https:\/\//);
   assert.ok(company[key].pages);
   assert.ok($('a[href="'+company[key].url+'"]').length,company.id+': missing evidence link');
  }
  for(const [base,current] of [['fy2024','fy2025'],['m92025','m92026']]){
   assert.equal(data.periods[base].months,data.periods[current].months);
   const row=$('tr[data-insurer="'+company.id+'"][data-current="'+current+'"]');
   assert.equal(row.attr('data-base'),base);
   for(const period of [base,current]){
    assert.ok(company.profit[period]>0);
    const displayed=Number(row.find('[data-value="'+period+'"]').text().replaceAll(',',''));
    assert.equal(displayed,company.profit[period]/100,company.id+': 百万元 to 亿元');
   }
   const rate=(company.profit[current]/company.profit[base]-1)*100;
   const chart=$('[data-chart-insurer="'+company.id+'"] [data-period="'+current+'"]');
   const expected=(rate>=0?'+':'−')+Math.abs(rate).toFixed(1)+'%';
   assert.equal(chart.find('strong').text(),expected);
   assert.equal(row.find('[data-growth]').text(),expected);
   const style=chart.find('i').attr('style');
   const width=Number(/--bar-width:([\d.]+)%/.exec(style)[1]);
   const left=Number(/--bar-left:([\d.]+)%/.exec(style)[1]);
   assert.ok(Math.abs(width-Math.abs(rate))<0.00001,'same −50%/+50% chart scale');
   assert.ok(Math.abs(left-(50+Math.min(0,rate)))<0.00001,'bar starts at zero');
  }
 }
});
