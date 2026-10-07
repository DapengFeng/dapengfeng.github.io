import {test} from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {DAY_MS,addDays,dateInfo,dailySceneAt,fourier,lissajous,standing,heatKernel,mobius,phyllotaxis} from '../src/scripts/daily-math.js';
import {loadDailyMath,validateCatalogue,validateHistory,planPublications,inventory,similarityWarnings} from '../scripts/daily-math.mjs';
import {bezier,convolution,taylorSine,uniformSumDensity,binomialPMF,poissonPMF,betaPDF,randomWalk} from '../src/scripts/daily-math-models.js';
const {topics,publications,entries}=await loadDailyMath();
const near=(a,b,tolerance=1e-10)=>assert.ok(Math.abs(a-b)<tolerance,`${a} ≠ ${b}`);

test('further reading uses existing articles and explicit bilingual labels',async()=>{
 const linked=topics.filter(topic=>topic.related?.length);
 assert.ok(linked.length>0);
 for(const topic of linked)for(const link of topic.related){
  assert.match(link.url,/^\/blog\/[a-z0-9-]+\.html$/);
  assert.ok(link.en?.trim()&&link.zh?.trim(),`${topic.id}: reading labels describe the destination in both languages`);
  await fs.access(new URL('../content/posts/'+link.url.slice('/blog/'.length),import.meta.url));
 }
});

test('exact date selection never loops, including at Shanghai midnight and after exhaustion',()=>{
 const midnight=Date.parse('2026-10-06T16:00:00Z');
 const before=dailySceneAt(midnight-1,entries),after=dailySceneAt(midnight,entries);
 assert.equal(before.date,'2026-10-06');assert.equal(after.date,'2026-10-07');
 assert.equal(before.nextChangeAt,midnight);assert.equal(after.nextChangeAt,midnight+DAY_MS);
 assert.notEqual(before.scene.id,after.scene.id);
 assert.equal(dailySceneAt(midnight-DAY_MS,entries).scene.id,before.scene.id);
 for(const date of [addDays(entries[0].date,-1),addDays(entries.at(-1).date,1)])assert.equal(dailySceneAt(Date.parse(date+'T04:00:00Z'),entries).scene,null);
 assert.throws(()=>dateInfo(NaN));assert.equal(addDays('2028-02-28',1),'2028-02-29');assert.equal(addDays('2028-02-29',1),'2028-03-01');assert.throws(()=>addDays('2027-02-29',1));
 // A decade-sized ledger is selected by its dates, never day-of-year or a cycle length.
 const decade=Array.from({length:3653},(_,i)=>({date:addDays('2026-10-06',i),id:`topic-${i}`}));
 for(const i of [0,365,730,1096,3652])assert.equal(dailySceneAt(Date.parse(decade[i].date+'T04:00:00Z'),decade).scene.id,`topic-${i}`);
 assert.equal(dailySceneAt(Date.parse(addDays(decade.at(-1).date,1)+'T04:00:00Z'),decade).scene,null);
});
test('publication validation rejects duplicates, missing review and mismatched concepts',()=>{
 for(const field of ['id','concept','en','zh','descriptionEn','descriptionZh']){
  const copy=structuredClone(topics);copy[1][field]=copy[0][field];assert.throws(()=>validateCatalogue(copy,publications),/Duplicate/);
 }
 const duplicate=[...publications,publications[0]];assert.throws(()=>validateCatalogue(topics,duplicate),/Duplicate/);
 const draft=structuredClone(topics);draft[0].status='draft';assert.throws(()=>validateCatalogue(draft,publications),/reviewed/);
 const wrong=structuredClone(publications);wrong[0].concept='another';assert.throws(()=>validateCatalogue(topics,wrong),/matching/);
 const missing=structuredClone(topics);delete missing[0].readingZh;assert.throws(()=>validateCatalogue(missing,publications),/Missing/);
 assert.deepEqual(similarityWarnings(topics),[]);
});
test('history is append-only for published dates, and planning never backfills missed days',()=>{
 const prior=publications.slice(0,3),today=prior[1].date;
 assert.doesNotThrow(()=>validateHistory(publications,prior,today));
 const removed=publications.slice(1);assert.throws(()=>validateHistory(removed,prior,today),/immutable/);
 const altered=structuredClone(publications);altered[0].id='replacement';assert.throws(()=>validateHistory(altered,prior,today),/immutable/);
 assert.throws(()=>validateHistory(publications,prior.slice(1),today),/backdate/);
 const next=planPublications(topics,prior,{today:'2026-10-07',count:2});assert.deepEqual(next.slice(0,3),prior);assert.equal(next.length,5);assert.equal(next[3].date,'2026-10-09');
 const resumed=planPublications(topics,prior,{today:'2026-11-10',count:1});assert.equal(resumed.at(-1).date,'2026-11-10');assert.doesNotThrow(()=>validateHistory(resumed,prior,'2026-11-10'));
 assert.equal(inventory(entries,entries[0].date).remaining,entries.length-1);assert.equal(inventory(entries,addDays(entries.at(-1).date,1)).days,0);
});
test('new diagram models preserve their defining identities and probability mass',()=>{
 near(convolution(0),0);near(convolution(1),1);near(convolution(1.7),.3);
 near(taylorSine(.6,5),Math.sin(.6),1e-9);
 const points=[[0,0],[1,2],[3,-1],[4,0]];assert.deepEqual(bezier(points,0),points[0]);assert.deepEqual(bezier(points,1),points[3]);
 for(const n of [1,2,4,8]){let mass=0,mean=0,variance=0;const dx=.002;for(let x=-5+dx/2;x<5;x+=dx){const v=uniformSumDensity(x,n)*dx;mass+=v;mean+=x*v;variance+=x*x*v;}near(mass,1,.002);near(mean,0,.001);near(variance,1,.004);}
 for(const n of [8,16,32,64]){let mass=0;for(let k=0;k<=n;k++)mass+=binomialPMF(k,n,3/n);near(mass,1,1e-10);}
 near(Array.from({length:30},(_,k)=>poissonPMF(k)).reduce((a,b)=>a+b),1,1e-10);
 for(const [a,b] of [[1,1],[7,3]]){let mass=0;for(let x=.0005;x<1;x+=.001)mass+=betaPDF(x,a,b)*.001;near(mass,1,1e-5);}
 const walk=randomWalk(7);assert.deepEqual(walk,randomWalk(7));for(let i=1;i<walk.length;i++)near(Math.hypot(walk[i][0]-walk[i-1][0],walk[i][1]-walk[i-1][1]),1);
 for(let n=0;n<20;n++)near(Math.sin(2*Math.PI*7*n/6),Math.sin(2*Math.PI*n/6));
});

test('the Fourier sum has the displayed odd spectrum and a closed period',()=>{
 near(fourier(.42),fourier(.42+2*Math.PI));near(fourier(-.42),-fourier(.42));
 for(const n of [1,2,3,4,13,15]){
  let coefficient=0;const samples=4096;
  for(let i=0;i<samples;i++){const t=2*Math.PI*i/samples;coefficient+=fourier(t)*Math.sin(n*t)*2/samples;}
  near(coefficient,n%2&&n<=13?4/(Math.PI*n):0,1e-12);
 }
});

test('Lissajous curve closes; membrane edges and nodal lines remain still',()=>{
 lissajous(0).forEach((x,i)=>near(x,lissajous(2*Math.PI)[i]));
 for(const t of [0,.3,4.7])for(const a of [.2,.8,2]){
  for(const x of [0,Math.PI/2,Math.PI])near(standing(x,a,t),0);
  for(const y of [0,Math.PI/3,2*Math.PI/3,Math.PI])near(standing(a,y,t),0);
 }
 near(standing(.7,.4,Math.PI/(2*Math.sqrt(13))),0);
});

test('heat diffusion conserves mass while its variance increases',()=>{
 for(const t of [.28,1,2.48]){
  let mass=0,secondMoment=0;const dx=.1;
  for(let x=-12+dx/2;x<12;x+=dx)for(let y=-12+dx/2;y<12;y+=dx){const density=heatKernel(x,y,t)*dx*dx;mass+=density;secondMoment+=(x*x+y*y)*density;}
  near(mass,1,1e-5);near(secondMoment,4*t,.001);
 }
 assert.ok(heatKernel(0,0,1)>heatKernel(0,0,2));
});

test('Möbius boundary closes after two turns, and golden-angle points occupy equal-area increments',()=>{
 const start=mobius(0,.7),half=mobius(2*Math.PI,.7),end=mobius(4*Math.PI,.7);
 start.forEach((x,i)=>near(x,end[i]));half.forEach((x,i)=>near(x,mobius(0,-.7)[i]));
 assert.ok(Math.hypot(...half.map((x,i)=>x-start[i]))>1);
 for(const n of [1,2,55,610]){const [x,y]=phyllotaxis(n);near(x*x+y*y,n,1e-9);}
});
