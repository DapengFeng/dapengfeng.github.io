import {mountProcesses} from './process-demos.js';
// Finite explanatory sequences. The complete diagram is the default; controls
// only inspect it. No measurements, compiler calls, or page scrolling happen here.
(() => {
 mountProcesses();
 const pair=(en,zh)=>`<span data-lang="en" lang="en">${en}</span><span data-lang="zh" lang="zh-CN">${zh}</span>`;
 const icons={
  play:'<path d="m8 5 11 7-11 7Z" fill="currentColor" stroke="none"/>',
  pause:'<path d="M8 5v14M16 5v14" stroke-width="4"/>',
  replay:'<path d="M4 10a8 8 0 1 1 1 7M4 4v6h6"/>',
  step:'<path d="m5 5 10 7-10 7Z" fill="currentColor" stroke="none"/><path d="M19 5v14"/>',
  complete:'<rect x="3" y="3" width="6" height="6" rx="1"/><rect x="15" y="3" width="6" height="6" rx="1"/><rect x="3" y="15" width="6" height="6" rx="1"/><path d="m14 18 3 3 5-6"/>'
 };
 function iconButton(button,icon,en,zh){
  button.dataset.icon=icon;
  button.innerHTML=`<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false">${icons[icon]}</svg>`;
  for(const attr of ['title','aria-label']){
   button.setAttribute(`data-${attr}-en`,en);button.setAttribute(`data-${attr}-zh`,zh);
   const lang=document.documentElement.dataset.language;
   button.setAttribute(attr,lang==='en'?en:lang==='zh'?zh:`${en} / ${zh}`);
  }
 }
 function explore(root,nodes,en,zh){
  if(!root||!nodes.length)return;const details=document.createElement('details');details.className='reading-explore';
  const summary=document.createElement('summary');summary.innerHTML=pair(en,zh);details.append(summary);
  nodes.forEach(node=>details.append(node));root.append(details);
 }
 const route=document.getElementById('pt-route');
 if(route)explore(route,[...route.children].filter(n=>!n.matches('.reading-overview')),'Inspect source nodes and compare execution settings','查看源码节点与执行条件');
 const tensor=document.getElementById('tensor-layout');
 if(tensor)explore(tensor,[...tensor.children].filter(n=>!n.matches('.reading-overview,[data-process]')),'Inspect any element, slice, or zero-stride view','查看任意元素、切片与零步幅视图');
 const fft=document.getElementById('fft-lab');
 if(fft)explore(fft,[...fft.children].filter(n=>!n.matches('#fft-comparison')&&n!==document.getElementById('fft-comparison').nextElementSibling),'Inspect samples and change the sine frequency','查看采样值与调整正弦频率');
 for(const id of ['lab-ultra','lab-bptt','lab-stdp']){
  const lab=document.getElementById(id);if(!lab)continue;
  const heading=lab.querySelector('.experiment-header');if(heading)lab.prepend(heading);
  const body=lab.querySelector('.experiment-body');if(body)explore(lab,[body,...lab.querySelectorAll(':scope > .experiment-caption')],'Explore other parameter values','探索其他参数值');
 }
 for(const id of ['lab-lif','lab-event']){
  const lab=document.getElementById(id);if(!lab)continue;
  const controls=lab.querySelector('.controls');
  if(controls)explore(lab,[...controls.querySelectorAll(':scope > .control-group'),...lab.querySelectorAll('#event-check,.process-reference')],'Adjust model parameters and inspect validation','调整模型参数与查看校验');
 }
 for(const id of ['band-lab','symmetric-lab','cuda-map']){
  const lab=document.getElementById(id);if(!lab)continue;
  explore(lab,[...lab.children].filter(el=>!el.matches('[data-process]')),'Change parameters and inspect individual elements','修改参数与查看具体元素');
 }
 const wave=document.querySelector('[data-wave-lab]');
 if(wave)explore(wave,[...wave.children].filter(el=>!el.matches('[data-process]')),'Explore other phase differences','探索其他相位差');
 const multiply=document.getElementById('pt-multiply');
 if(multiply){
  multiply.querySelector('.pt-label').textContent='Y = AB · (2 × 3) @ (3 × 2)';
  explore(multiply,[multiply.querySelector('.lesson-controls'),multiply.querySelector('.pt-toolbar'),multiply.querySelector('#pt-matrix-position')],'Change inputs or inspect every multiply–add','修改输入或查看每次乘加');
 }
 const fw=document.getElementById('fw-lab');
 if(fw){
  explore(fw,[...fw.children].filter(el=>!el.matches('[data-process]')),'Change the target or inspect iterations','更换目标或查看迭代步骤');
 }
 // Connections follow the actual layout, including the one-column mobile view.
 for(const flow of document.querySelectorAll('.reading-overview .reading-flow')){
  const cards=[...flow.children],svg=document.createElementNS('http://www.w3.org/2000/svg','svg');
  svg.classList.add('reading-connections');svg.setAttribute('aria-hidden','true');flow.prepend(svg);
  const marker='reading-arrow-'+flow.closest('[id]').id;
  const draw=()=>{const box=flow.getBoundingClientRect();if(!box.width)return;
   svg.setAttribute('viewBox',`0 0 ${box.width} ${box.height}`);
   svg.innerHTML=`<defs><marker id="${marker}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0 0L10 5L0 10Z" fill="#a1bf91"/></marker></defs>`+cards.slice(1).map((card,i)=>{
    const a=cards[i].getBoundingClientRect(),b=card.getBoundingClientRect();let path;
    if(Math.abs(a.top-b.top)<5)path=`M${a.right-box.left+2},${a.top+a.height/2-box.top}H${b.left-box.left-3}`;
    else{const x=a.left+a.width/2-box.left,y=a.bottom-box.top+2,X=b.left+b.width/2-box.left,Y=b.top-box.top-3,mid=(y+Y)/2;path=`M${x},${y}V${mid}H${X}V${Y}`;}
    return `<path class="reading-connection" data-active="${card.dataset.active==='true'}" data-to="${i+1}" d="${path}" marker-end="url(#${marker})"/>`;
   }).join('');
  };new ResizeObserver(draw).observe(flow);window.addEventListener('languagechange',draw);draw();
 }
 const motion=matchMedia('(prefers-reduced-motion: reduce)');
 for(const root of document.querySelectorAll('#matrix-lab,#fw-lab,#pt-route,#pt-multiply,#tensor-layout,.memory-demo,#band-lab,#symmetric-lab,#cuda-map,#lab-lif,[data-wave-lab]')){
  const demo=root.readingDemo;if(!demo)continue;
  const overview=root.querySelector('.reading-overview');
  if(overview&&!demo.tick){demo.target=overview;const original=demo.render;demo.render=n=>{original(n);overview.querySelectorAll('[data-scene]').forEach(el=>el.dataset.active=String(Number(el.dataset.scene)===n));overview.querySelectorAll('.reading-connection').forEach(el=>el.dataset.active=String(Number(el.dataset.to)===n));};}
  const bar=document.createElement('div');bar.className='reading-playback';
  const play=demo.button?demo.button.cloneNode(false):document.createElement('button');play.type='button';play.className='outline-button';play.disabled=false;
  if(demo.button)demo.button.replaceWith(play);
  const end=document.createElement('button');end.type='button';end.className='outline-button';iconButton(end,'complete','Full diagram','完整示意');
  bar.append(play,end);
  if(demo.barAfter)demo.barAfter.after(bar);else if(overview)overview.after(bar);else if(root.matches('.memory-demo'))root.append(bar);else demo.target.after(bar);
  const position=document.createElement('span');position.className='reading-position';bar.append(position);
  let frame=0,timer=null,raf=null,elapsed=0,lastTime=null,finished=true,playing=false,visible=false,used=false;
  const live=[...root.querySelectorAll('[aria-live]')].map(el=>[el,el.getAttribute('aria-live')]);
  const clear=()=>{clearTimeout(timer);cancelAnimationFrame(raf);timer=null;raf=null;lastTime=null;demo.activity?.(false);};
  function label(){live.forEach(([el,value])=>el.setAttribute('aria-live',playing?'off':value));iconButton(play,playing?'pause':finished?'replay':'play',playing?'Pause':finished?'Replay':'Resume',playing?'暂停':finished?'重播':'继续');play.setAttribute('aria-pressed',String(playing));root.dataset.demoPlaying=String(playing);}
  function render(){elapsed=0;demo.render(frame);demo.tick?.(motion.matches?1:0,frame);position.textContent=`${frame+1} / ${demo.frames}`;root.dataset.demoFrame=String(frame);}
  function stop(){clear();playing=false;label();}
  function finish(){finished=true;elapsed=0;frame=demo.frames-1;stop();demo.complete?demo.complete():demo.render(demo.frames-1);overview?.querySelectorAll('[data-scene],.reading-connection').forEach(el=>delete el.dataset.active);position.innerHTML=pair('Complete','完整');root.dataset.demoFrame='complete';label();}
  function schedule(){
   clear();if(!playing||!visible||document.hidden)return;demo.activity?.(true);
   if(!demo.tick){timer=setTimeout(()=>{if(frame>=demo.frames-1){finish();return;}frame++;render();label();schedule();},demo.interval||1200);return;}
   const tick=now=>{if(!playing||!visible||document.hidden)return;
    if(lastTime!==null)elapsed+=Math.max(0,now-lastTime);lastTime=now;
    const duration=demo.interval||1600,p=Math.min(1,elapsed/duration);demo.tick(motion.matches?1:p,frame);
    if(p>=1){if(frame>=demo.frames-1){finish();return;}frame++;render();label();}
    raf=requestAnimationFrame(tick);
   };raf=requestAnimationFrame(tick);
  }
  function start(){if(demo.canPlay&&!demo.canPlay())return;used=true;finished=false;playing=true;frame=0;render();label();schedule();}
  play.addEventListener('click',()=>{used=true;if(playing){stop();return;}if(!finished){playing=true;label();schedule();}else start();});
  end.addEventListener('click',()=>{used=true;finish();});
  if(demo.tick){const step=document.createElement('button');step.type='button';step.className='outline-button';iconButton(step,'step','Step','单步');end.before(step);step.addEventListener('click',()=>{used=true;stop();frame=finished?0:Math.min(demo.frames-1,frame+1);finished=false;render();elapsed=demo.interval||1600;demo.tick(1,frame);label();});}
  root.addEventListener('processchange',()=>{used=true;finish();});
  // Reader input takes ownership; automatic playback never overwrites edits.
  for(const event of ['pointerdown','keydown','input'])root.addEventListener(event,e=>{if(bar.contains(e.target))return;used=true;finished=true;frame=demo.frames-1;stop();demo.cancel?.();},true);
  root.addEventListener('click',e=>{if(!bar.contains(e.target)){used=true;finished=true;frame=demo.frames-1;stop();demo.cancel?.();}},true);
  const observer=new IntersectionObserver(entries=>{const entry=entries[0];visible=entry.isIntersecting&&entry.intersectionRect.height>=Math.min(160,entry.boundingClientRect.height*.2);if(!visible){clear();return;}if(!used&&!motion.matches&&!document.hidden)start();else schedule();},{threshold:[0,.01,.05,.1,.2,.35,.5,1]});
  observer.observe(demo.target);
  document.addEventListener('visibilitychange',()=>{if(document.hidden)clear();else if(!used&&visible&&!motion.matches)start();else schedule();});
  motion.addEventListener('change',()=>{if(motion.matches){used=true;finish();}});
  window.addEventListener('pagehide',stop);
  window.addEventListener('beforeprint',finish);
  // Do not trigger a network request or change a user's input to make a diagram.
  finish();
 }
})();
