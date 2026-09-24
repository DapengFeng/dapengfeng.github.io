import {test} from 'node:test';
import assert from 'node:assert/strict';
import {frankWolfeStep,cudaWrites,lifCycle,bandAddress,symmetricAddress,transposeOrder,superposedWave} from '../src/scripts/process-models.js';
test('Frank–Wolfe oracle and exact line search stay feasible and minimize the segment',()=>{
 for(const x of [[0,0],[.8,0],[.3,.4]])for(const c of [[.8,.6],[.25,.35]]){
  const s=frankWolfeStep(x,c);assert.ok(s.next.every(v=>v>=-1e-12));assert.ok(s.next[0]+s.next[1]<=1+1e-12);
  assert.ok(s.objective(s.gamma)<=s.objective(0)+1e-12);
  for(let k=0;k<=100;k++)assert.ok(s.objective(s.gamma)<=s.objective(k/100)+1e-12);
 }
 const first=frankWolfeStep([0,0],[.8,.6]);assert.deepEqual(first.s,[1,0]);assert.equal(first.gamma,.8);assert.deepEqual(first.next,[.8,0]);assert.ok(Math.abs(first.objective(.8)-.18)<1e-12);
});
test('CUDA address-count model distinguishes guarded tails from missing offsets',()=>{
 for(let n=1;n<=64;n++)for(const width of [2,4,8,16]){const m=cudaWrites(n,width);assert.ok(m.writes.every(v=>v===1));assert.equal(m.masked,Math.ceil(n/width)*width-n);}
 const bad=cudaWrites(17,8,true);assert.equal(bad.missing,9);assert.equal(bad.collisions,8);assert.equal(bad.writes.reduce((a,b)=>a+b),24);
});
test('LIF event time satisfies the threshold, then resets and holds; silent inputs never fire',()=>{
 const m=lifCycle(1.6,20,1);assert.ok(Math.abs(1.6*(1-Math.exp(-m.tint/20))-1)<1e-12);assert.ok(m.potential(m.tint-1e-5)>.99999);
 assert.equal(m.potential(m.tint),0);assert.equal(m.potential(m.tint+1),0);assert.equal(m.horizon-m.tint,2);
 for(const J of [.8,1]){const silent=lifCycle(J,20,1);assert.equal(silent.spike,false);assert.equal(silent.rate,0);assert.ok(silent.potential(silent.horizon)<1);}
});
test('Packed examples and transpose preserve values and logical addresses',()=>{
 assert.deepEqual(bandAddress(2,1,6,1,1),{row:2,column:1,offset:5});assert.equal(bandAddress(0,3,6,1,1),null);
 const slots=[];for(let i=0;i<6;i++)for(let j=0;j<6;j++){const a=bandAddress(i,j,6,1,1);if(a)slots.push(a.offset);}assert.equal(new Set(slots).size,16);
 assert.equal(symmetricAddress(1,2,4).offset,5);assert.deepEqual(symmetricAddress(1,2,4),symmetricAddress(2,1,4));
 const order=transposeOrder(3,4);assert.deepEqual(order,[0,4,8,1,5,9,2,6,10,3,7,11]);assert.equal(order.indexOf(6),7);
});
test('Same-phase waves double amplitude and opposite-phase waves cancel at every sampled time',()=>{
 for(const t of [0,.7,2*Math.PI])for(let x=0;x<8;x+=.1){assert.ok(Math.abs(superposedWave(x,t,Math.PI))<1e-12);assert.ok(Math.abs(superposedWave(x,t,0)-2*Math.sin(1.5*x-t))<1e-12);}
});
