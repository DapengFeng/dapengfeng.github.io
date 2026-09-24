import {frankWolfeStep,cudaWrites,lifCycle,bandAddress,symmetricAddress,transposeOrder,superposedWave} from './process-models.js';
const pair=(en,zh)=>`<span data-lang="en" lang="en">${en}</span><span data-lang="zh" lang="zh-CN">${zh}</span>`;
const q=(root,s)=>root.querySelector(s),all=(root,s)=>[...root.querySelectorAll(s)];
const set=(node,text)=>{if(node.textContent!==text)node.textContent=text;};
const label=(figure,en,zh)=>q(figure,'.process-stage-label').innerHTML=pair(en,zh);
function particles(figure){
 const layer=q(figure,'.process-particles'),positions=new Map();
 const refresh=()=>positions.forEach(args=>api.move(...args));
 new ResizeObserver(refresh).observe(figure);
 window.addEventListener('languagechange',refresh);
 const api={
  clear(){positions.clear();layer.replaceChildren();},
  move(key,value,from,to,p){
   if(!from||!to)return;positions.set(key,[key,value,from,to,p]);let token=q(layer,`[data-token="${key}"]`);
   if(!token){token=document.createElement('span');token.className='process-token';token.dataset.token=key;token.textContent=value;layer.append(token);}
   const r=figure.getBoundingClientRect(),a=from.getBoundingClientRect(),b=to.getBoundingClientRect(),t=p*p*(3-2*p);
   const x=a.left+a.width/2+(b.left+b.width/2-a.left-a.width/2)*t-r.left;
   const y=a.top+a.height/2+(b.top+b.height/2-a.top-a.height/2)*t-r.top-18*Math.sin(Math.PI*p);
   token.style.transform=`translate(${x}px,${y}px) translate(-50%,-50%)`;token.style.opacity=p>=1?'0':'1';
  }
 };
 return api;
}
function highlight(figure,nodes){all(figure,'[data-active]').forEach(n=>n.removeAttribute('data-active'));nodes.filter(Boolean).forEach(n=>n.dataset.active='true');}
function chart(svg,height,xmax,ymin,ymax){
 const width=Math.max(240,Math.round(svg.clientWidth||480)),L=42,R=16,T=22,B=32;
 svg.setAttribute('viewBox',`0 0 ${width} ${height}`);
 const X=x=>L+x/xmax*(width-L-R),Y=y=>height-B-(y-ymin)/(ymax-ymin)*(height-T-B);
 return {width,height,X,Y,axes:`<path class="process-axis" d="M${L} ${T}V${height-B}H${width-R}"/>`,
  point:(x,y,color='#eff4e8',r=6)=>`<circle cx="${X(x)}" cy="${Y(y)}" r="${r}" fill="${color}"/>`,
  text:(x,y,t)=>`<text x="${X(x)}" y="${Y(y)}" text-anchor="middle">${t}</text>`,
  path:(points,color='#c0f47b',extra='')=>`<path d="M${points.map(([x,y])=>`${X(x)},${Y(y)}`).join('L')}" fill="none" stroke="${color}" stroke-width="3" ${extra}/>`};
}
function resize(figure,draw){new ResizeObserver(()=>draw()).observe(figure);window.addEventListener('languagechange',draw);}
function mountFW(root,figure){
 const state=frankWolfeStep([0,0],[.8,.6]);let phase=4,progress=1;
 const names=[['1 · Choose the minimizing vertex','1 · 选择使线性目标最小的顶点'],['2 · Form the feasible search segment','2 · 确定可行搜索线段'],['3 · Minimize the objective along the segment','3 · 在线段上最小化目标'],['4 · Update the iterate','4 · 更新迭代点']];
 function draw(){
  const a=chart(q(figure,'[data-plot="fw-geometry"]'),280,1.15,0,1.1),b=chart(q(figure,'[data-plot="fw-search"]'),170,1,.15,.55);
  const searching=phase===2,gamma=phase<2?0:phase===2?(progress<.65?progress/.65:1-(progress-.65)/.35*(1-state.gamma)):state.gamma;
  const moved=phase===3?progress:phase===4?1:0,x=.8*moved;
  q(figure,'[data-plot="fw-geometry"]').innerHTML=a.axes+
   `<path d="M${a.X(0)},${a.Y(0)}L${a.X(1)},${a.Y(0)}L${a.X(0)},${a.Y(1)}Z" fill="#c0f47b15" stroke="#8daa7f"/>`+
   a.path([[0,0],[.8*(phase===0?progress:1),.6*(phase===0?progress:1)]],'#8dbbf8','stroke-dasharray="5 5"')+a.point(.8,.6,'#8dbbf8')+
   a.path([[0,0],[phase===1?progress:1,0]],phase===0?'#61765c':'#e9b77c')+a.point(1,0,'#e9b77c',phase===0?4+4*progress:8)+
   (searching?a.point(gamma,0,'#bdabff',7):'')+a.point(0,0,'#81937f',4)+a.point(x,0)+
   a.text(.8,.72,'c')+a.text(.9,.09,'s')+a.text(.42,.4,'−∇f(x)')+a.text(1,-.1,'1')+a.text(0,1.04,'1')+a.text(.04,-.1,'0')+a.text(1.08,.09,'x₁')+a.text(.18,1.04,'x₂');
  const points=Array.from({length:61},(_,i)=>[i/60,state.objective(i/60)]);
  q(figure,'[data-plot="fw-search"]').innerHTML=b.axes+b.path(points)+b.path([[.8,.15],[.8,.5]],'#e9b77c','stroke-dasharray="3 5"')+b.point(gamma,state.objective(gamma),'#bdabff')+
   b.text(0,.105,'0')+b.text(.8,.105,'0.8')+b.text(1,.105,'γ')+b.text(.11,.51,'0.5')+b.text(.8,.24,'min');
  figure.dataset.gamma=gamma.toFixed(6);figure.dataset.iterate=x.toFixed(6);
 }
 root.readingDemo={target:figure,barAfter:figure,frames:4,interval:1700,render(n){phase=n;label(figure,...names[n]);draw();},tick(p){progress=p;draw();},complete(){phase=4;progress=1;label(figure,'One iteration: vertex → segment → line search → update','一次迭代：顶点 → 线段 → 线搜索 → 更新');draw();}};
 resize(figure,draw);
}
function mountTensor(root,figure){
 const fly=particles(figure),order=transposeOrder(3,4),names=[['Transpose logical positions; buffer A does not move','转置逻辑位置；存储 A 不搬动'],['Follow the view address: a.T[2,1] → A[6]','沿视图寻址：a.T[2,1] → A[6]'],['Materialize contiguous storage: copy A into a new B','物化连续存储：从 A 复制到新的 B']];
 root.readingDemo={target:figure,barAfter:figure,frames:3,interval:1900,
  render(n){fly.clear();figure.dataset.stage=n;label(figure,...names[n]);highlight(figure,[q(figure,'[data-logical-a="6"]'),q(figure,'[data-logical-t="7"]'),q(figure,'[data-buffer-a="6"]'),...(n===2?[q(figure,'[data-buffer-b="7"]')]:[])]);},
  tick(p,n){if(n===0)order.forEach((v,i)=>fly.move(i,v,q(figure,`[data-logical-a="${v}"]`),q(figure,`[data-logical-t="${i}"]`),p));
   else if(n===1)fly.move('lookup',6,q(figure,'[data-logical-t="7"]'),q(figure,'[data-buffer-a="6"]'),p);
   else order.forEach((v,i)=>fly.move(i,v,q(figure,`[data-buffer-a="${v}"]`),q(figure,`[data-buffer-b="${i}"]`),p));},
  cancel(){fly.clear();},complete(){fly.clear();figure.dataset.stage='complete';label(figure,'A stays fixed; only contiguous() allocates and copies B','A 保持不动；只有 contiguous() 分配并复制 B');highlight(figure,all(figure,'[data-value="6"]'));}};
}
function mountBand(root,figure){
 const fly=particles(figure),names=[['Upper diagonal → packed row 0','上对角线 → 紧凑存储第 0 行'],['Main diagonal → packed row 1','主对角线 → 紧凑存储第 1 行'],['Lower diagonal → packed row 2','下对角线 → 紧凑存储第 2 行']];
 function entries(n){return Array.from({length:6},(_,j)=>[j+n-1,j]).filter(([i])=>i>=0&&i<6);}
 root.readingDemo={target:figure,barAfter:figure,frames:3,interval:1700,render(n){fly.clear();label(figure,...names[n]);highlight(figure,entries(n).flatMap(([i,j])=>[q(figure,`[data-band-a="${i*6+j}"]`),q(figure,`[data-band-b="${n*6+j}"]`)]));},
  tick(p,n){for(const [i,j]of entries(n)){const packed=bandAddress(i,j,6,1,1);fly.move(j,i===j?2:-1,q(figure,`[data-band-a="${i*6+j}"]`),q(figure,`[data-band-b="${packed.row*6+j}"]`),p);}},
  cancel(){fly.clear();},complete(){fly.clear();label(figure,'Columns stay fixed: r = 1 + i − j, offset = 3j + r','列号不变：r = 1 + i − j，offset = 3j + r');highlight(figure,[q(figure,'[data-band-a="13"]'),q(figure,'[data-band-b="13"]')]);}};
}
function mountSymmetric(root,figure){
 const fly=particles(figure),slot=q(figure,'[data-sym-slot]'),names=[['Two equal entries refer to one stored 22; no addition','两个相等元素对应一个存储值 22，不相加'],['Read 22 once; multiply by x₂ and x₁','读取一次 22，分别乘以 x₂ 与 x₁'],['Accumulate both off-diagonal contributions','累加两路非对角贡献']];
 figure.dataset.offset=symmetricAddress(1,2,4).offset;
 root.readingDemo={target:figure,barAfter:figure,frames:3,interval:1700,render(n){fly.clear();label(figure,...names[n]);highlight(figure,n===0?[q(figure,'[data-sym-a="6"]'),q(figure,'[data-sym-a="9"]'),slot]:all(figure,n===1?'[data-product]':'[data-output]'));},
  tick(p,n){if(n===0)[6,9].forEach(i=>fly.move(i,22,q(figure,`[data-sym-a="${i}"]`),slot,p));else [0,1].forEach(i=>fly.move(i,n===1?22:[66,44][i],n===1?slot:q(figure,`[data-product="${i}"]`),q(figure,n===1?`[data-product="${i}"]`:`[data-output="${i}"] b`),p));},
  cancel(){fly.clear();},complete(){fly.clear();label(figure,'One stored value, two multiply–adds: Δy₁ = 66, Δy₂ = 44','一个存储值，两次乘加：Δy₁ = 66，Δy₂ = 44');highlight(figure,[slot,q(figure,'[data-sym-a="6"]'),q(figure,'[data-sym-a="9"]'),...all(figure,'[data-output]')]);}};
}
function mountCUDA(root,figure){
 const fly=particles(figure);let previous=-1;
 function counts(blocks){if(previous===blocks)return;previous=blocks;
  for(const mode of ['correct','broken']){const model=cudaWrites(17,8,mode==='broken',blocks);figure.dataset[mode+'Writes']=JSON.stringify(model.writes);
   all(figure,`[data-gpu="${mode}"] [data-index]`).forEach((cell,i)=>{cell.dataset.writes=model.writes[i];set(q(cell,'b'),model.writes[i]+'×');});}
 }
 root.readingDemo={target:figure,barAfter:figure,frames:3,interval:1800,
  render(n){fly.clear();previous=-1;counts(n);label(figure,n===2?'Tail block: one valid address, seven guarded threads':'Map block '+n+' to its output addresses',n===2?'尾块：一个有效地址，七个线程被屏蔽':'将块 '+n+' 映射到输出地址');
   for(const mode of ['correct','broken'])q(figure,`[data-launch="${mode}"]`).textContent=`b = ${n} · t = 0…7${n===2&&mode==='correct'?' · × × × × × × ×':''}`;
   all(figure,'.process-tile-strip span').forEach((el,i)=>el.dataset.active=String(i===n));},
  tick(p,n){for(const mode of ['correct','broken'])for(let t=0;t<8;t++){const i=mode==='broken'?t:n*8+t;if(i<17)fly.move(mode+t,i,q(figure,`[data-launch="${mode}"]`),q(figure,`[data-gpu="${mode}"] [data-index="${i}"]`),p);}counts(n+(p>=.98?1:0));},
  cancel(){fly.clear();},complete(){fly.clear();previous=-1;counts(3);label(figure,'Correct: 17 unique outputs · broken: 9 missing, 8 multiply written','正确：17 个独立输出 · 错误：9 个漏写，8 个重复写入');for(const mode of ['correct','broken'])q(figure,`[data-launch="${mode}"]`).textContent='b = 0, 1, 2 · t = 0…7';highlight(figure,[]);}};
}
function mountLIF(root,figure){
 const plotArea=q(root,'.plot-area'),reference=document.createElement('div');reference.className='process-reference';
 [...plotArea.children].forEach(el=>reference.append(el));plotArea.append(figure,reference);
 let phase=4,progress=1;
 const inputs=()=>[Number(q(root,'#lif-input').value),Number(q(root,'#lif-tau').value),Number(q(root,'#lif-theta').value)];
 function draw(){
  const [J,tau,theta]=inputs(),m=lifCycle(J,tau,theta),svg=q(figure,'[data-plot="lif-cycle"]'),a=chart(svg,250,m.horizon,0,Math.max(1.2,theta*1.22));
  let t=m.spike?(phase===0?m.tint*progress:phase===1||phase===2?m.tint:phase===3?m.tint+m.ref*progress:m.horizon):m.horizon*(phase===0?progress:1);
  let u=m.spike&&phase===0?J*(1-Math.exp(-t/tau)):m.spike&&phase===1?theta:m.potential(t);
  const end=m.spike?m.tint:m.horizon,pts=Array.from({length:101},(_,i)=>{const x=end*i/100;return [x,J*(1-Math.exp(-x/tau))];});
  const integrated=[...pts];if(m.spike)pts.push([m.tint,0],[m.horizon,0]);
  const trail=phase===0?Array.from({length:101},(_,i)=>{const x=t*i/100;return [x,J*(1-Math.exp(-x/tau))];}):phase===1?integrated:phase===2?[...integrated,[m.tint,0]]:phase===3?[...integrated,[m.tint,0],[t,0]]:pts;
  svg.innerHTML=a.axes+a.path([[0,theta],[m.horizon,theta]],'#e9b77c','stroke-dasharray="5 5"')+a.path(pts,'#82987a')+
   a.path(trail,'#c0f47b')+
   (m.spike?`<path d="M${a.X(m.tint)} ${a.Y(theta)}v${phase===1?-18-10*progress:-18}" stroke="#e9b77c" stroke-width="4"/>`:'')+
   a.point(t,u,phase===1?'#e9b77c':'#c0f47b',7)+a.text(0,-.12,'0')+a.text(m.horizon/2,-.12,(m.horizon/2).toFixed(1))+a.text(m.horizon,-.12,m.horizon.toFixed(1))+a.text(m.horizon*.93,theta+.06,'θ')+a.text(m.horizon*.7,Math.max(1.12,theta*1.17),'u(t)');
  q(figure,'[data-lif-values]').textContent=`t = ${t.toFixed(3)} ms · u = ${u.toFixed(3)}${m.spike?` · t* = ${m.tint.toFixed(3)} ms`:' · J ≤ θ'}`;
  figure.dataset.time=t;figure.dataset.potential=u;figure.dataset.spiking=String(m.spike);figure.dataset.reset=String(m.spike&&phase>=2);
 }
 const names=[['Integrate: u(t) = J(1 − exp(−t/τ))','积分：u(t) = J(1 − exp(−t/τ))'],['Reach the threshold and emit a spike','到达阈值，产生脉冲'],['Instantaneous reset: same event time, u → 0','瞬时复位：事件时刻不变，u → 0'],['Hold u = 0 for the 2 ms refractory period','在 2 ms 不应期内保持 u = 0']];
 root.readingDemo={target:figure,barAfter:figure,get frames(){return lifCycle(...inputs()).spike?4:1;},interval:1700,
  render(n){phase=n;label(figure,...(lifCycle(...inputs()).spike?names[n]:['Subthreshold input: integrate without firing or resetting','阈下输入：积分，不发放也不复位']));},tick(p){progress=p;draw();},complete(){phase=4;progress=1;label(figure,...(lifCycle(...inputs()).spike?['One event: integrate → spike → reset → 2 ms hold','一次事件：积分 → 脉冲 → 复位 → 保持 2 ms']:['J ≤ θ: no finite threshold crossing, no reset','J ≤ θ：有限时间内不越阈，不复位']));draw();}};
 root.addEventListener('input',()=>root.dispatchEvent(new Event('processchange')));
 all(root,'[data-lif-preset]').forEach(button=>button.addEventListener('click',()=>root.dispatchEvent(new Event('processchange'))));
 resize(figure,draw);
}
function mountWaves(root,figure){
 let time=0;
 function draw(){for(const [key,phase]of [['same',0],['opposite',Math.PI]]){
  const svg=q(figure,`[data-plot="wave-${key}"]`),a=chart(svg,200,8,-2.3,2.3),samples=Array.from({length:161},(_,i)=>i/20);
  svg.innerHTML=a.axes+a.path([[0,0],[8,0]],'#61765c')+a.path(samples.map(x=>[x,Math.sin(1.5*x-time)]),'#8dbbf8')+
   a.path(samples.map(x=>[x,Math.sin(1.5*x-time+phase)]),'#e9b77c','stroke-dasharray="5 5"')+a.path(samples.map(x=>[x,superposedWave(x,time,phase)]))+
   a.text(0,-2.8,'0')+a.text(4,-2.8,'4')+a.text(8,-2.8,'x')+a.text(.1,2,'2')+a.text(.1,-2,'−2');
 }figure.dataset.waveTime=time;}
 root.readingDemo={target:figure,barAfter:figure,frames:1,interval:4800,render(){time=0;draw();},tick(p){time=2*Math.PI*p;draw();},complete(){time=0;draw();}};
 resize(figure,draw);
}
const chartDescriptions={
 'fw-geometry':['Frank–Wolfe: a feasible triangle, target c, selected vertex s, and an iterate moving along their feasible segment.','Frank–Wolfe：三角可行域、目标 c、选中的顶点 s，以及沿可行线段移动的迭代点。'],
 'fw-search':['Objective versus step size gamma; the minimum on this segment occurs at gamma 0.8.','目标值随步长 gamma 变化；该线段上的极小值位于 gamma 为 0.8 处。'],
 'lif-cycle':['Membrane potential versus time in milliseconds. A suprathreshold input reaches threshold, emits a spike, resets at the same event time, then holds at zero. Subthreshold inputs do not fire.','膜电位随时间（毫秒）变化。超阈输入到达阈值后发放，在同一事件时刻复位，再保持为零。阈下输入不发放。'],
 'wave-same':['Equal-amplitude waves in phase: their sum has twice the amplitude. The blue and dashed orange components coincide.','等振幅的同相波：叠加振幅加倍。蓝色分量与橙色虚线分量重合。'],
 'wave-opposite':['Equal-amplitude waves with phase difference pi cancel at every position; the green sum remains zero.','等振幅的两列波相差 pi，在所有位置相消；绿色叠加曲线始终为零。']
};
function describeCharts(){const lang=document.documentElement.dataset.language;
 for(const svg of document.querySelectorAll('[data-plot]')){const labels=chartDescriptions[svg.dataset.plot];if(labels)svg.setAttribute('aria-label',lang==='en'?labels[0]:lang==='zh'?labels[1]:labels.join(' / '));}
}
export function mountProcesses(){
 describeCharts();window.addEventListener('languagechange',describeCharts);
 for(const [name,mount]of Object.entries({fw:mountFW,tensor:mountTensor,band:mountBand,symmetric:mountSymmetric,cuda:mountCUDA,lif:mountLIF,waves:mountWaves})){
  const figure=document.querySelector(`[data-process="${name}"]`);if(!figure)continue;mount(figure.closest('.lesson-lab,.tl-demo,.experiment,[data-wave-lab]'),figure);
 }
}
