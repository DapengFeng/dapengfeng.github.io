import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {load} from 'cheerio';
import {collectSeries} from '../scripts/series.mjs';
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
 assert.equal($('.series-episodes li').length,1);
 assert.equal($('.series-episodes li').attr('value'),'1');
 assert.equal($('.series-episodes a').attr('href'),'/blog/pytorch-01-what-is-pytorch.html');
 assert.ok((await fs.readFile('dist/sitemap.xml','utf8')).includes('/series/pytorch-internals/'));
 for(const file of ['index.html','blog/index.html','blog/pytorch-01-what-is-pytorch.html']){
  const page=load(await fs.readFile('dist/'+file,'utf8'));
  assert.ok(page('a[href="/series/pytorch-internals/"]').length,file);
 }
});
