import {fourier} from './daily-math.js';
import {taylorSine,convolution,lerp,bezier,normal,uniformSumDensity,poissonPMF,binomialPMF,betaPDF,randomWalk} from './daily-math-models.js';
const TAU=Math.PI*2,GREEN='#c5ec96',GOLD='#e5bf82',CYAN='#83bec8',FAINT='#9ab79d45';
const walks=Array.from({length:8},(_,i)=>randomWalk(41+i*97));
// Draw in a shared mathematical coordinate system; equal units stay equal on both axes.
export function drawDailyConstruction(ctx,w,h,id,time){
 const scale=Math.min(w/13,h/7),cx=w/2,cy=h/2;
 const p=([x,y])=>[cx+x*scale,cy-y*scale];
 function line(points,color=GREEN,width=1.8,fill){
  ctx.beginPath();points.forEach((v,i)=>{const [x,y]=p(v);i?ctx.lineTo(x,y):ctx.moveTo(x,y);});
  if(fill){ctx.closePath();ctx.fillStyle=fill;ctx.fill();}ctx.strokeStyle=color;ctx.lineWidth=width;ctx.stroke();
 }
 function dot(v,color=GREEN,r=3){const [x,y]=p(v);ctx.beginPath();ctx.arc(x,y,r,0,TAU);ctx.fillStyle=color;ctx.fill();}
 function label(v,text,color=GOLD){const [x,y]=p(v);ctx.font='16px ui-monospace,monospace';ctx.fillStyle=color;ctx.textAlign='left';ctx.fillText(text,x,y);}
 function curve(fn,a,b,color=GREEN,n=240){line(Array.from({length:n+1},(_,i)=>fn(a+(b-a)*i/n)),color);}
 function circle(center,r,color=FAINT){curve(t=>[center[0]+r*Math.cos(t),center[1]+r*Math.sin(t)],0,TAU,color,120);}
 function axes(x=5,y=2.4){line([[-x,0],[x,0]],FAINT,1);line([[0,-y],[0,y]],FAINT,1);}
 function plot(fn,a,b,color=GREEN){curve(x=>[x,fn(x)],a,b,color);}
 const phase=(time*.22)%1,progress=(1-Math.cos(time*.35))/2;
 ctx.save();
 if(id==='euler'){
  const center=[-3,0],r=1.7,a=time*.65,end=[center[0]+r*Math.cos(a),r*Math.sin(a)];
  circle(center,r);line([center,end],GOLD);line([end,[end[0],0]],CYAN);line([end,[0,end[1]]],FAINT,1);dot(end);axes(5,2.4);
  curve(x=>[x,r*Math.sin(a-x)],0,5);dot([0,end[1]]);label([-3.2,-2.25],'eⁱθ');label([2,2.35],'Im');
 }else if(id==='taylor'){
  axes();plot(Math.sin,-4.5,4.5,CYAN);const n=1+Math.floor(time*.5)%5;
  const points=[];for(let x=-4.5;x<=4.5;x+=.025){const y=taylorSine(x,n);if(Math.abs(y)<=2.7)points.push([x,y]);else if(points.length){line(points);points.length=0;}}if(points.length)line(points);
  label([-4.7,2.7],`m = ${n}`);
 }else if(id==='beats'){
  axes(5.5);const t=time*.2;
  plot(x=>Math.sin(10*(x+t))+Math.sin(11*(x+t)),-5.5,5.5);
  plot(x=>2*Math.abs(Math.cos((x+t)/2)),-5.5,5.5,GOLD);plot(x=>-2*Math.abs(Math.cos((x+t)/2)),-5.5,5.5,GOLD);
 }else if(id==='gibbs'){
  axes(5.2);line([[-5,-1],[0,-1],[0,1],[5,1]],FAINT);const n=[3,11,31][Math.floor(time*.4)%3];
  plot(x=>fourier(x*.55,n),-5,5);label([-4.7,2.35],`N = ${n}`);
 }else if(id==='convolution'){
  const t=-.2+phase*2.4,x=v=>-4+v*3,lo=Math.max(0,t-1),hi=Math.min(1,t);
  line([[x(-.5),.4],[x(2.5),.4]],FAINT);line([[x(0),.4],[x(0),1.8],[x(1),1.8],[x(1),.4]],CYAN);
  line([[x(t-1),.4],[x(t-1),2.1],[x(t),2.1],[x(t),.4]],GOLD);
  if(hi>lo)line([[x(lo),.4],[x(hi),.4],[x(hi),1.8],[x(lo),1.8]],GREEN,1,'#c5ec9630');
  curve(v=>[x(v),-2.2+convolution(v)*1.5],-.2,2.2);dot([x(t),-2.2+convolution(t)*1.5]);label([-4.6,2.7],`t = ${t.toFixed(2)}`);
 }else if(id==='gaussian'){
  const sigma=.5+progress*1.2;
  for(const [center,s,color,title] of [[-3,sigma,GREEN,'x'],[3,1/sigma,GOLD,'ω']]){
   line([[center-2.5,-1.5],[center+2.5,-1.5]],FAINT);curve(x=>[center+x,-1.5+2.7*Math.exp(-x*x/(2*s*s))],-2.5,2.5,color);label([center,1.8],title,color);
  }label([-1,2.8],`σ = ${sigma.toFixed(2)}`);
 }else if(id==='aliasing'){
  axes(5.4);const x=t=>-5+t*5;curve(t=>[x(t),1.6*Math.sin(TAU*t)],0,2,GREEN,500);curve(t=>[x(t),1.6*Math.sin(7*TAU*t)],0,2,CYAN,700);
  for(let n=0;n<=12;n++)dot([x(n/6),1.6*Math.sin(TAU*n/6)],GOLD,3.2);
  const n=Math.floor(time*2)%13;dot([x(n/6),1.6*Math.sin(TAU*n/6)],GOLD,6);label([-4.7,2.5],'1 Hz / 7 Hz · fₛ = 6 Hz');
 }else if(id==='bezier'){
  const points=[[-4,-1.8],[-2,2.6],[1,-2.2],[4,1.7]],t=progress;
  line(points,FAINT);points.forEach(v=>dot(v,CYAN));curve(t=>bezier(points,t),0,1);
  const a=points.slice(1).map((v,i)=>lerp(points[i],v,t)),b=a.slice(1).map((v,i)=>lerp(a[i],v,t));line(a,GOLD);line(b,CYAN);a.forEach(v=>dot(v,GOLD));b.forEach(v=>dot(v,CYAN));dot(bezier(points,t),GREEN,5);
 }else if(id==='cycloid'){
  const t=phase*TAU,at=t=>[t-Math.sin(t)-Math.PI,1-Math.cos(t)-1];
  line([[-5,-1],[5,-1]],FAINT);curve(at,0,TAU,FAINT);curve(at,0,t);circle([t-Math.PI,0],1,CYAN);line([[t-Math.PI,0],at(t)],GOLD);dot(at(t),GREEN,5);
 }else if(id==='spiral'){
  const at=a=>[.1*Math.exp(.2*a)*Math.cos(a),.1*Math.exp(.2*a)*Math.sin(a)],a=phase*16;
  curve(at,0,16,FAINT);curve(at,0,a);line([[0,0],at(a)],GOLD);dot(at(a));
 }else if(id==='ellipse'){
  const a=3.6,b=1.8,c=Math.sqrt(a*a-b*b),v=[a*Math.cos(time*.5),b*Math.sin(time*.5)],d=Math.hypot(v[0]+c,v[1]);
  curve(t=>[a*Math.cos(t),b*Math.sin(t)],0,TAU);line([[-c,0],v],GREEN);line([[c,0],v],GOLD);dot([-c,0]);dot([c,0],GOLD);dot(v,CYAN,5);
  line([[-a,-2.5],[-a+d,-2.5]],GREEN,5);line([[-a+d,-2.5],[a,-2.5]],GOLD,5);label([-.5,2.6],'2a');
 }else if(id==='involute'){
  const at=t=>[Math.cos(t)+t*Math.sin(t),Math.sin(t)-t*Math.cos(t)],t=phase*3.6;
  circle([0,0],1,CYAN);curve(at,0,3.6,FAINT);curve(at,0,t);const v=[Math.cos(t),Math.sin(t)];line([v,at(t)],GOLD);dot(v,GOLD);dot(at(t));
 }else if(['determinant','eigenvectors','svd'].includes(id)){
  axes();const rotate=(v,a)=>[v[0]*Math.cos(a)-v[1]*Math.sin(a),v[0]*Math.sin(a)+v[1]*Math.cos(a)];
  let map;
  if(id==='determinant'){
   const q=progress;map=([x,y])=>[(1+.3*q)*x+.6*q*y,.4*q*x+(1+.1*q)*y];
   for(let k=-3;k<=3;k++){line([map([-3,k]),map([3,k])],FAINT,1);line([map([k,-2]),map([k,2])],FAINT,1);}
   line([[0,0],[1,0],[1,1],[0,1],[0,0]].map(map),GREEN,2,'#c5ec9620');label([-4,2.65],`det = ${(1+.4*q-.21*q*q).toFixed(2)}`);
  }else{
   if(id==='eigenvectors'){const q=progress;map=([x,y])=>[(1+q)*x+q*y,q*x+(1+q)*y];}
   else{const stage=(time*.25)%3;map=v=>{v=rotate(v,-.65*Math.min(1,stage));if(stage>1){const q=Math.min(1,stage-1);v=[v[0]*(1+q),v[1]*(1-.4*q)];}return stage>2?rotate(v,.9*(stage-2)):v;};label([-4,2.6],['Vᵀ','Σ','U'][Math.floor(stage)]);}
   circle([0,0],1,FAINT);for(let i=0;i<16;i++){const a=i/16*TAU,v=[Math.cos(a),Math.sin(a)];line([[0,0],v],FAINT,1);line([[0,0],map(v)],i%8===2?GREEN:i%8===6?GOLD:FAINT,i%4===2?2.5:1);}
   curve(a=>map([Math.cos(a),Math.sin(a)]),0,TAU,CYAN);dot(map([1,0]),GOLD,5);
  }
 }else if(id==='gradient'){
  for(const radius of [.4,.8,1.3,2,3,4])curve(a=>[radius*Math.cos(a),radius/2*Math.sin(a)],0,TAU,FAINT);
  const n=Math.floor(time*1.5)%17,pts=Array.from({length:n+1},(_,k)=>[2.6*.8**k,1.5*.2**k]);line(pts);pts.forEach(v=>dot(v));dot(pts.at(-1),GOLD,5);label([-4,2.6],`k = ${n}`);
 }else if(id==='newton'){
  axes(4);const f=x=>x*x-2;plot(f,-2.1,2.1,CYAN);let x=2.4;
  const n=Math.floor(time*.7)%5;
  for(let i=0;i<n;i++)x=(x+2/x)/2;
  const next=(x+2/x)/2,y=f(x);line([[x,y],[next,0]],GOLD,2);line([[x,0],[x,y]],FAINT);dot([x,y]);dot([next,0],GOLD,4);label([-4,2.6],`x = ${x.toFixed(5)}`);
 }else if(id==='bisection'){
  const f=x=>x*x*x-x-2,x=v=>(v-1.5)*7,y=v=>v*.5;let lo=1,hi=2;
  const n=Math.floor(time*.8)%7;for(let i=0;i<n;i++){const mid=(lo+hi)/2;if(f(mid)<0)lo=mid;else hi=mid;}
  line([[x(lo),-2],[x(hi),-2],[x(hi),2.3],[x(lo),2.3]],FAINT,1,'#c5ec9620');line([[-4,0],[4,0]],FAINT);curve(v=>[x(v),y(f(v))],1,2);const mid=(lo+hi)/2;dot([x(mid),y(f(mid))],GOLD,5);label([-4,2.7],`[${lo.toFixed(3)}, ${hi.toFixed(3)}]`);
 }else if(id==='riemann'){
  const n=[4,8,16,32][Math.floor(time*.5)%4],x=v=>-4+v*8,y=v=>-1.5+v*3,f=v=>1/(1+v*v);
  for(let k=0;k<n;k++){const a=k/n,b=(k+1)/n;line([[x(a),y(0)],[x(b),y(0)],[x(b),y(f(a))],[x(a),y(f(a))]],GREEN,1,'#c5ec9618');}
  curve(v=>[x(v),y(f(v))],0,1,GOLD);label([-4,2.6],`n = ${n}`);
 }else if(id==='divergence'||id==='curl'){
  circle([0,0],1.7,GOLD);
  for(let j=0;j<18;j++)for(let k=1;k<=4;k++){
   const radius=id==='curl'?k*.65:.4*Math.exp(((time*.3+k*.5)%2)),angle=j/18*TAU+(id==='curl'?time*.45:0);
   const at=t=>[radius*Math.exp(id==='divergence'?t:0)*Math.cos(angle+(id==='curl'?t:0)),radius*Math.exp(id==='divergence'?t:0)*Math.sin(angle+(id==='curl'?t:0))];
   line([at(-.14),at(0)],FAINT,1.3);dot(at(0),k%2?GREEN:CYAN,2);
  }
 }else if(id==='walk'){
  const n=1+Math.floor(time*7)%120,factor=.16;circle([0,0],Math.sqrt(n)*factor,GOLD);
  walks.forEach((path,i)=>{const color=i%2?CYAN:GREEN,pts=path.slice(0,n+1).map(([x,y])=>[x*factor,y*factor]);line(pts,color,1);dot(pts.at(-1),color);});label([-4,2.65],`n = ${n}`);
 }else if(id==='clt'){
  const n=1+Math.floor(time*.55)%8,at=(x,d)=>[x*1.25,-1.4+d*6];line([[-5,-1.4],[5,-1.4]],FAINT);curve(x=>at(x,normal(x)),-4,4,GOLD);
  if(n===1){const r=Math.sqrt(3);line([at(-4,0),at(-r,0),at(-r,1/(2*r)),at(r,1/(2*r)),at(r,0),at(4,0)]);}
  else curve(x=>at(x,uniformSumDensity(x,n)),-4,4);
  label([-4.7,2.7],`n = ${n}`);
 }else if(id==='poisson'){
  const n=[8,16,32,64][Math.floor(time*.5)%4];for(let k=0;k<=10;k++){const x=-4.5+k*.9,v=binomialPMF(k,n,3/n)*13;line([[x,-1.5],[x,-1.5+v]],GREEN,9);dot([x,-1.5+poissonPMF(k)*13],GOLD,3.5);if(k%2===0)label([x-.12,-1.95],String(k),CYAN);}label([-4.7,2.6],`n = ${n}, np = 3`);
 }else if(id==='bayes'){
  const sequence=[1,1,0,1,0,1,1,1],n=Math.floor(time*.65)%9,heads=sequence.slice(0,n).reduce((a,b)=>a+b,0),a=1+heads,b=1+n-heads;
  const at=(x,y)=>[-4+x*8,-1.5+y];line([at(0,1),at(1,1)],GOLD);curve(x=>at(x,betaPDF(x,a,b)),0,1);
  for(let k=0;k<8;k++)dot([-3+k*.8,2.8],k<n?(sequence[k]?GREEN:CYAN):FAINT,4);
  label([-4,-2.2],'0');label([4,-2.2],'1');
 }
 ctx.restore();
}
