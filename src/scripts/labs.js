(() => {
 const reduced=matchMedia('(prefers-reduced-motion: reduce)').matches;
 function plot(canvas,series,{xmin=0,xmax=100,ymin=-1,ymax=1}={}){
  const rect=canvas.getBoundingClientRect();if(!rect.width)return;
  const w=rect.width,h=rect.height,dpr=Math.min(devicePixelRatio||1,2);canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);
  const c=canvas.getContext('2d');c.setTransform(dpr,0,0,dpr,0,0);const X=x=>42+(x-xmin)/(xmax-xmin)*(w-60),Y=y=>h-30-(y-ymin)/(ymax-ymin)*(h-50);
  c.font='10px monospace';for(let i=0;i<=4;i++){const y=ymin+(ymax-ymin)*i/4;c.strokeStyle='#34402c';c.lineWidth=.6;c.beginPath();c.moveTo(42,Y(y));c.lineTo(w-18,Y(y));c.stroke();c.fillStyle='#95a28a';c.textAlign='right';c.fillText(Math.abs(y)<.01?'0':y.toFixed(1),35,Y(y)+3);const x=xmin+(xmax-xmin)*i/4;c.textAlign='center';c.fillText(x.toFixed(0),X(x),h-10);}
  for(const s of series){c.beginPath();let started=false;for(const [x,y]of s.points){if(!Number.isFinite(x)||!Number.isFinite(y)){started=false;continue;}if(!started){c.moveTo(X(x),Y(y));started=true;}else c.lineTo(X(x),Y(y));}c.strokeStyle=s.color||'#c0f47b';c.lineWidth=s.width||2;c.setLineDash(s.dash||[]);c.stroke();c.setLineDash([]);}
 }
 function bilingual(en,zh){
  const pair=document.createElement('span');pair.className='i18n';
  for(const [language,text]of [['en',en],['zh',zh]]){const span=document.createElement('span');span.dataset.lang=language;span.lang=language==='zh'?'zh-CN':'en';span.textContent=text;pair.append(span);}return pair;
 }
 document.querySelectorAll('[data-wave-lab]').forEach(lab=>{
  const canvas=lab.querySelector('canvas'),phase=lab.querySelector('[data-wave-phase]'),wave=lab.querySelector('[data-wave-k]'),button=lab.querySelector('[data-wave-pause]');
  let paused=reduced,time=0,last=0,frame=0,visible=false;
  function draw(){const k=Number(wave.value),p=Number(phase.value),one=[],two=[],sum=[];for(let x=0;x<=12;x+=.03){const a=Math.sin(k*x-time),b=Math.sin(k*x-time+p);one.push([x,a]);two.push([x,b]);sum.push([x,a+b]);}plot(canvas,[{points:one,color:'#7eabdd',width:1},{points:two,color:'#78976a',width:1},{points:sum}],{xmax:12,ymin:-2.2,ymax:2.2});lab.querySelector('output').value=p.toFixed(2)+' rad';}
  function tick(now){frame=0;if(paused||!visible||document.hidden)return;if(last)time+=Math.min(40,now-last)*.001;last=now;draw();frame=requestAnimationFrame(tick);}
  function start(){cancelAnimationFrame(frame);last=0;draw();if(!paused&&visible&&!document.hidden)frame=requestAnimationFrame(tick);}
  function label(){button.replaceChildren(bilingual(paused?'Play':'Pause',paused?'播放':'暂停'));button.setAttribute('aria-pressed',String(paused));}
  function canvasLabel(){const lang=document.documentElement.dataset.language;canvas.setAttribute('aria-label',lang==='en'?canvas.dataset.labelEn:lang==='zh'?canvas.dataset.labelZh:`${canvas.dataset.labelEn} / ${canvas.dataset.labelZh}`);}
  phase.addEventListener('input',draw);wave.addEventListener('input',draw);
  button.addEventListener('click',()=>{paused=!paused;label();start();});
  new ResizeObserver(draw).observe(canvas);new IntersectionObserver(e=>{visible=e[0].isIntersecting;start();}).observe(lab);
  document.addEventListener('visibilitychange',start);window.addEventListener('languagechange',canvasLabel);label();canvasLabel();start();
 });
 document.querySelectorAll('[data-benchmark-lab]').forEach(lab=>{
 const view={run:lab.querySelector('[data-benchmark-run]'),exp:lab.querySelector('[data-benchmark-export]'),size:lab.querySelector('select'),result:lab.querySelector('[data-benchmark-result]')};
 let saved=null,busy=false,status='idle',errorMessage='';
 const messages={idle:['Not run yet. Results will be measured on this device.','尚未运行。数据来自当前设备的实际测量。'],running:['Warming up and measuring…','预热与测量中…'],timeout:['Measurement timed out. Please retry.','测量超时，请重试。'],failed:['Experiment failed. Refresh and try again.','实验运行失败，请刷新后重试。'],unsupported:['This browser cannot run this experiment.','当前浏览器不支持此实验。']};
 const workerErrors={'Invalid input size':'输入规模无效','Output verification failed':'输出校验失败','Output verification failed during timing':'计时期间输出校验失败'};
 function renderBenchmark(){
  const {run,exp,size,result}=view;
  run.disabled=busy;size.disabled=busy;exp.disabled=!saved||busy;
  if(status==='error'){result.replaceChildren(bilingual(errorMessage,workerErrors[errorMessage]||messages.failed[1]));return;}
  if(!saved){result.replaceChildren(bilingual(...messages[status]));return;}
  const a=saved.summary.plain,b=saved.summary.unrolled,ratio=b.median>0?a.median/b.median:null;
  result.replaceChildren();
  function row(en,zh,value){const line=document.createElement('div');line.className='benchmark-row';line.append(bilingual(en,zh),typeof value==='string'?document.createTextNode(value):value);result.append(line);}
  row('Plain loop / ms per call','普通循环 / 毫秒每次',a.median.toFixed(4));
  row('Unrolled loop / ms per call','四路展开 / 毫秒每次',b.median.toFixed(4));
  row('Plain IQR / ms','普通循环四分位区间 / 毫秒',`[${a.q1.toFixed(4)}, ${a.q3.toFixed(4)}]`);
  row('Unrolled IQR / ms','四路展开四分位区间 / 毫秒',`[${b.q1.toFixed(4)}, ${b.q3.toFixed(4)}]`);
  row('Time ratio (plain / unrolled)','耗时比（普通 / 展开）',ratio?ratio.toFixed(2)+'×':bilingual('Below timer resolution','计时精度不足'));
  row('Samples × calls per sample','采样轮数 × 每轮次数','15 × 16');
  result.append(bilingual('Verified · measured on this device','校验通过 · 当前设备测量'));
 }
 view.size.addEventListener('change',()=>{saved=null;status='idle';renderBenchmark();});
  view.run.addEventListener('click',()=>{
   if(busy)return;saved=null;busy=true;status='running';renderBenchmark();let worker,timer;
   const finish=next=>{clearTimeout(timer);worker?.terminate();busy=false;status=next;renderBenchmark();};
   try{
    worker=new Worker('/assets/benchmark-worker.js');timer=setTimeout(()=>finish('timeout'),30000);
    worker.onmessage=e=>{if(e.data.error){errorMessage=e.data.error;finish('error');return;}saved={...e.data,userAgent:navigator.userAgent,logicalProcessors:navigator.hardwareConcurrency??null};finish('idle');};
    worker.onerror=()=>finish('failed');worker.postMessage({size:Number(view.size.value)});
   }catch{finish('unsupported');}
  });
  view.exp.addEventListener('click',()=>{if(!saved)return;const url=URL.createObjectURL(new Blob([JSON.stringify(saved,null,2)],{type:'application/json'})),a=document.createElement('a');a.href=url;a.download='feng-benchmark-'+saved.timestamp.slice(0,10)+'.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
 renderBenchmark();
 });
})();
