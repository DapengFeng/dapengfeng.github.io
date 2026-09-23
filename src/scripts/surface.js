(() => {
 const canvas=document.getElementById('surface-canvas');if(!canvas)return;
 const ctx=canvas.getContext('2d'),slider=document.getElementById('surface-frequency'),out=document.getElementById('frequency-value'),button=document.getElementById('surface-toggle');
 const reduced=matchMedia('(prefers-reduced-motion: reduce)');
 let paused=reduced.matches,time=0,last=0,frame=0,visible=true,yaw=-.52,drag=null;
 function project(x,y,z,w,h){const a=x*Math.cos(yaw)-y*Math.sin(yaw),b=x*Math.sin(yaw)+y*Math.cos(yaw);return [w*.5+a*w*.066,h*.48+b*h*.063-z*h*.21];}
 function draw(){
  const rect=canvas.getBoundingClientRect(),w=rect.width,h=rect.height,dpr=Math.min(devicePixelRatio||1,2);if(!w||!h)return;
  if(canvas.width!==Math.round(w*dpr)||canvas.height!==Math.round(h*dpr)){canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);}
  ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
  function path(points,color,width=1){ctx.beginPath();points.forEach((p,i)=>{const q=project(...p,w,h);i?ctx.lineTo(...q):ctx.moveTo(...q);});ctx.strokeStyle=color;ctx.lineWidth=width;ctx.stroke();}
  for(let k=-5;k<=5;k++) {path([[k,-5,-.85],[k,5,-.85]],'#34472b55',.5);path([[-5,k,-.85],[5,k,-.85]],'#34472b55',.5);}
  path([[-5,-5,-.85],[5,-5,-.85]],'#72856388',.8);path([[5,-5,-.85],[5,5,-.85]],'#72856388',.8);
  const freq=Number(slider.value),n=42;
  for(let i=0;i<=n;i++){
    const v=-4.4+i*8.8/n,points=[],other=[];
    for(let j=0;j<=n;j++){const u=-4.4+j*8.8/n,r=Math.hypot(u,v),z=Math.sin(freq*r-time)*Math.exp(-.12*r*r);points.push([u,v,z]);other.push([v,u,z]);}
    const brightness=Math.round(90+90*(1-Math.abs(v)/4.4));path(points,`rgba(${brightness*.82},${brightness+42},${brightness*.5},.85)`,.8);path(other,`rgba(163,213,113,${.2+(1-Math.abs(v)/4.4)*.3})`,.6);
  }
  for(const [label,x,y,z]of [['x',5.25,-5,-.85],['y',5,5.3,-.85],['z',0,0,1.15]]){ctx.fillStyle='#9cac85';ctx.font='16px monospace';ctx.fillText(label,...project(x,y,z,w,h));}
 }
 function tick(now){frame=0;if(!visible||document.hidden||paused)return;if(last)time+=Math.min(40,now-last)*.0007;last=now;draw();frame=requestAnimationFrame(tick);}
 function start(){cancelAnimationFrame(frame);frame=0;last=0;draw();if(!paused&&visible&&!document.hidden)frame=requestAnimationFrame(tick);}
 function updateButton(){button.textContent=paused?'▶':'Ⅱ';button.setAttribute('aria-pressed',String(paused));const zh=document.documentElement.dataset.language==='zh';button.setAttribute('aria-label',paused?(zh?'播放曲面动画':'Play surface animation'):(zh?'暂停曲面动画':'Pause surface animation'));}
 button.addEventListener('click',()=>{paused=!paused;updateButton();start();});
 slider.addEventListener('input',()=>{out.value=Number(slider.value).toFixed(1);draw();});
 canvas.addEventListener('pointerdown',e=>{if(e.pointerType==='touch')return;drag=e.clientX;canvas.setPointerCapture(e.pointerId);});canvas.addEventListener('pointermove',e=>{if(drag===null)return;yaw+=(e.clientX-drag)*.007;drag=e.clientX;draw();});canvas.addEventListener('pointerup',()=>drag=null);canvas.addEventListener('pointercancel',()=>drag=null);
 new ResizeObserver(draw).observe(canvas);new IntersectionObserver(entries=>{visible=entries[0].isIntersecting;start();}).observe(canvas);document.addEventListener('visibilitychange',start);window.addEventListener('languagechange',()=>{updateButton();draw();});reduced.addEventListener('change',()=>{paused=reduced.matches;updateButton();start();});updateButton();start();
})();
