import {dateInfo,fourier,lissajous,standing,heatKernel,mobius,phyllotaxis} from './daily-math.js';

import {drawDailyConstruction} from './daily-math-drawings.js';
import {renderers} from './daily-math-models.js';

const root=document.querySelector('[data-daily-math]');
if(root)mount(root);
function mount(root){
 const canvas=root.querySelector('#surface-canvas'),ctx=canvas.getContext('2d');
 if(!ctx)return;
 const reduced=matchMedia('(prefers-reduced-motion: reduce)');
 let current,visible=true,time=1.4,last=0,frame=0,dateTimer=0,width=0,height=0;
 const note=root.querySelector('.daily-math-note'),fixed=root.dataset.mathFixedDate;
 const initial=JSON.parse(root.querySelector('[data-math-entry]').textContent),initialHTML=note.innerHTML;
 const cache=new Map(initial.id?[[initial.date,{...initial,html:initialHTML}]]:[]);
 let request=0,controller,pendingDate;
 const empty='<div class="daily-math-identity"><h2 id="daily-math-title">'+window.LabI18n.bi('Mathematics archive','数学往期')+'</h2><a class="daily-math-source" href="/math/">'+window.LabI18n.bi('Browse the collection','浏览往期内容')+'</a></div>';
 function show(entry,date){
  current={date,scene:entry};time=1.4;root.dataset.date=date;root.dataset.scene=entry?.id||'';
  note.innerHTML=entry?.html||empty;draw();
 }
 async function selectDate(){
  clearTimeout(dateTimer);const calendar=dateInfo(),date=fixed||calendar.date;
  if(!fixed)dateTimer=setTimeout(selectDate,Math.min(3600000,Math.max(100,calendar.nextChangeAt-Date.now()+30)));
  if(current?.date===date&&current.scene||pendingDate===date)return;
  controller?.abort();const token=++request;controller=new AbortController();pendingDate=date;
  show(null,date); // Never display yesterday's construction under today's date.
  try{
   let entry=cache.get(date);
   if(!entry){
    const response=await fetch(`/assets/daily-math/${date}.json`,{signal:controller.signal,cache:'no-cache'});
    if(!response.ok)throw Error('Daily entry unavailable');
    entry=await response.json();
    if(entry.date!==date||!entry.id||!renderers.includes(entry.renderer)||typeof entry.html!=='string')throw Error('Invalid daily entry');
   }
   if(token!==request)return;
   cache.set(date,entry);if(cache.size>4)cache.delete(cache.keys().next().value);
   show(entry,date);
  }catch(error){if(token===request&&error.name!=='AbortError')show(null,date);}
  finally{if(token===request)pendingDate=null;}
 }
 function line(points,color='#bcec8d',lineWidth=1,close=false){ctx.beginPath();points.forEach(([x,y],i)=>i?ctx.lineTo(x,y):ctx.moveTo(x,y));if(close)ctx.closePath();ctx.strokeStyle=color;ctx.lineWidth=lineWidth;ctx.stroke();}
 function dot(x,y,r=3,color='#e0ffc0'){ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.fillStyle=color;ctx.shadowColor=color;ctx.shadowBlur=12;ctx.fill();ctx.shadowBlur=0;}
 function circle(x,y,r,color){ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.strokeStyle=color;ctx.lineWidth=1;ctx.stroke();}
 function drawFourier(w,h){
  const scale=Math.min(w*.14,h*.25),cx=w*.24,cy=h*.43,start=w*.56,length=w*.41;
  let x=cx,y=cy;
  for(let k=0;k<7;k++){const n=2*k+1,r=4/Math.PI/n*scale,nx=x+r*Math.cos(n*time),ny=y-r*Math.sin(n*time);circle(x,y,r,k?'#8fb77b70':'#dcb47fbd');line([[x,y],[nx,ny]],k?'#aacd8e':'#ecd1a7',1.4);x=nx;y=ny;}
  ctx.setLineDash([3,6]);line([[x,y],[start,y]],'#badb9a7f');ctx.setLineDash([]);dot(x,y,3.5);
  line([[start,cy],[start+length,cy]],'#abc99e30');
  const wave=[];for(let i=0;i<=220;i++)wave.push([start+length*i/220,cy-scale*fourier(time-i/220*Math.PI*2)]);line(wave,'#bfea99',2);dot(start,y,3);
  const base=h*.86,barW=length/11;
  for(let k=0;k<7;k++){const n=2*k+1,bx=start+k*length/7;ctx.fillStyle=k?'#8db88c88':'#dec28e';ctx.fillRect(bx,base-h*.17/n,barW,h*.17/n);ctx.fillStyle='#bed0ab';ctx.font='16px monospace';ctx.textAlign='center';if(w>420||k%2===0)ctx.fillText(String(n),bx+barW/2,base+21);}
 }
 function drawLissajous(w,h){
  const at=t=>{const [x,y]=lissajous(t);return [w*.5+x*w*.42,h*.48-y*h*.38];};
  line([[w*.05,h*.48],[w*.95,h*.48]],'#91b79728');line([[w*.5,h*.06],[w*.5,h*.9]],'#91b79728');
  line(Array.from({length:501},(_,i)=>at(i/500*Math.PI*2)),'#84b8b173',1.4);
  for(let j=0;j<24;j++)line(Array.from({length:8},(_,i)=>at(time*.6-(24-j)*.035+i*.005)),`rgba(207,242,146,${.2+j/30})`,2.5);
  dot(...at(time*.6),4);
 }
 function drawSurface(w,h,id){
  const yaw=-.6,scale=Math.min(w/11,h/6.9);
  const project=(x,y,z)=>[w*.5+(x*Math.cos(yaw)-y*Math.sin(yaw))*scale,h*.54+(x*Math.sin(yaw)+y*Math.cos(yaw))*scale*.46-z*scale];
  const extent=id==='standing'?Math.PI:id==='heat'?4.5:2.5,n=42;
  const heatTime=.28+((time+1)%18)/18*2.2;
  const value=(x,y)=>id==='standing'?standing(x,y,time*.5)*1.15:id==='heat'?heatKernel(x,y,heatTime)*13:(x*x-y*y)*.25;
  const at=(x,y)=>project(id==='standing'?(x-Math.PI/2)*2.5:x,id==='standing'?(y-Math.PI/2)*2.5:y,value(x,y));
  const low=id==='standing'?0:-extent,high=extent;
  for(let i=0;i<=n;i++){
   const v=low+(high-low)*i/n,a=[],b=[];
   for(let j=0;j<=n;j++){const u=low+(high-low)*j/n;a.push(at(u,v));b.push(at(v,u));}
   const alpha=.3+.5*Math.sin(i/n*Math.PI);line(a,`rgba(180,225,139,${alpha})`,.9);line(b,'#91bba87a',.7);
  }
  if(id==='standing'){
   for(const x of [Math.PI/2])line([at(x,0),at(x,Math.PI)],'#e5c995',1.8);
   for(const y of [Math.PI/3,2*Math.PI/3])line([at(0,y),at(Math.PI,y)],'#9bcaca',1.8);
  }
  if(id==='saddle'){
   for(const axis of [0,1]){const points=[];for(let i=0;i<=80;i++){const a=-2.5+i/16;points.push(at(axis?0:a,axis?a:0));}line(points,axis?'#86c6d5':'#e9c78b',2);const a=Math.sin(time*.6)*2.5;dot(...at(axis?0:a,axis?a:0),3,axis?'#86c6d5':'#e9c78b');}
  }
 }
 function drawMobius(w,h){
  const scale=Math.min(w/6.6,h/4.8),yaw=.35+Math.sin(time*.11)*.18;
  const at=(u,v)=>{const [x,y,z]=mobius(u,v),a=x*Math.cos(yaw)-y*Math.sin(yaw),b=x*Math.sin(yaw)+y*Math.cos(yaw);return [w*.5+a*scale,h*.52+b*scale*.45-z*scale,b];};
  const tiles=[];
  for(let i=0;i<100;i++)for(let j=0;j<8;j++){const u=i/100*Math.PI*2,v=-.7+j*.175;const points=[at(u,v),at(u+Math.PI/50,v),at(u+Math.PI/50,v+.175),at(u,v+.175)];tiles.push({points,depth:points.reduce((s,p)=>s+p[2],0)/4,i,j});}
  tiles.sort((a,b)=>a.depth-b.depth);
  for(const {points,i,j} of tiles){ctx.beginPath();points.forEach(([x,y],k)=>k?ctx.lineTo(x,y):ctx.moveTo(x,y));ctx.closePath();ctx.fillStyle=`hsla(${88+j*7},30%,${18+Math.sin(i/100*Math.PI)*18}%,.85)`;ctx.fill();ctx.strokeStyle='#b2d39635';ctx.lineWidth=.65;ctx.stroke();}
  const edge=[];for(let i=0;i<=400;i++)edge.push(at(i/400*Math.PI*4,.7));line(edge,'#cbdca187',1.2);
  const t=time*.7;line(Array.from({length:80},(_,i)=>at(t-i*.012,.7)),'#ebcc96',2.2);dot(...at(t,.7).slice(0,2),4,'#ffe2ad');
 }
 function drawPhyllotaxis(w,h){
  const scale=Math.min(w,h)*.44/Math.sqrt(610),head=Math.floor(time*35)%610;
  for(let n=1;n<=610;n++){const [x,y]=phyllotaxis(n),bright=n<=head,alpha=bright?.88:.2;ctx.beginPath();ctx.arc(w*.5+x*scale,h*.48+y*scale,Math.max(1.1,scale*.55),0,Math.PI*2);ctx.fillStyle=`hsla(${75+n/610*85},50%,${bright?68:45}%,${alpha})`;ctx.fill();}
  const p=phyllotaxis(Math.max(1,head));dot(w*.5+p[0]*scale,h*.48+p[1]*scale,3,'#f0d5a5');
 }
 function draw(){
  if(!current)return;
  const box=canvas.getBoundingClientRect();if(!box.width||!box.height)return;
  width=box.width;height=box.height;const dpr=Math.min(devicePixelRatio||1,1.75,Math.sqrt(8000000/(width*height)));
  if(canvas.width!==Math.round(width*dpr)||canvas.height!==Math.round(height*dpr)){canvas.width=Math.round(width*dpr);canvas.height=Math.round(height*dpr);}
  ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,width,height);
  if(!current.scene){root.dataset.rendered="true";return;}
  // The mathematical construction belongs to the entire hero background, independent of text flow.
  for(let x=20;x<width;x+=48)for(let y=24;y<height;y+=48){ctx.fillStyle='#a5bf8a18';ctx.fillRect(x,y,1.2,1.2);}
  const archived=Boolean(fixed),compact=width<=960;
  // Keep geometry at a stable aspect ratio. Portrait screens crop the atmosphere,
  // rather than stretching curves or reserving an illustration panel in the layout.
  const sceneWidth=archived?width*.96:compact?Math.max(620,width*1.35):Math.max(width*.94,Math.min(height*1.65,width*1.14));
  const sceneHeight=sceneWidth/1.8;
  const centerX=width*(archived?.5:compact?.56:.61),centerY=archived?height*.5:compact?Math.min(height*.32,260):height*.40;
  ctx.save();ctx.translate(centerX-sceneWidth/2,centerY-sceneHeight/2);
  if(current.scene.renderer==='fourier')drawFourier(sceneWidth,sceneHeight);
  else if(current.scene.renderer==='lissajous')drawLissajous(sceneWidth,sceneHeight);
  else if(['standing','heat','saddle'].includes(current.scene.renderer))drawSurface(sceneWidth,sceneHeight,current.scene.renderer);
  else if(current.scene.renderer==='mobius')drawMobius(sceneWidth,sceneHeight);
  else if(current.scene.renderer==='phyllotaxis')drawPhyllotaxis(sceneWidth,sceneHeight);
  else drawDailyConstruction(ctx,sceneWidth,sceneHeight,current.scene.renderer,time);
  ctx.restore();root.dataset.rendered='true';
 }
 function tick(now){frame=0;if(reduced.matches||!visible||document.hidden)return;if(!last||now-last>=1000/30){if(last)time+=Math.min(now-last,100)*.00045;last=now;draw();}frame=requestAnimationFrame(tick);}
 function resume(){cancelAnimationFrame(frame);frame=0;last=0;draw();if(!reduced.matches&&visible&&!document.hidden)frame=requestAnimationFrame(tick);}
 const header=document.querySelector('.lab-header');
 const resize=new ResizeObserver(()=>{
  const headerHeight=header?.getBoundingClientRect().height||0;
  if(root.style.getPropertyValue('--home-header-height')!==`${headerHeight}px`)root.style.setProperty('--home-header-height',`${headerHeight}px`);
  draw();
 });
 resize.observe(root);if(header)resize.observe(header);
 new IntersectionObserver(([entry])=>{visible=entry.isIntersecting;resume();}).observe(root);
 document.addEventListener('visibilitychange',()=>{if(!document.hidden)selectDate();resume();});
 window.addEventListener('languagechange',draw);
 reduced.addEventListener('change',resume);
 window.addEventListener('pagehide',()=>{controller?.abort();pendingDate=null;request++;clearTimeout(dateTimer);cancelAnimationFrame(frame);frame=0;});
 window.addEventListener('pageshow',()=>{selectDate();resume();});
 selectDate();resume();
}
