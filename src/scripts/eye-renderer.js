import * as T from './three.module.min.js';
import {createEye,eyeParts} from './eye-model.js';

const pair=(en,zh)=>`<span data-lang="en" lang="en">${en}</span><span data-lang="zh" lang="zh-CN">${zh}</span>`;
export async function mountEye(root){
 const stage=root.querySelector('.eye-stage'),canvas=root.querySelector('canvas'),poster=root.querySelector('.eye-poster');
 let renderer;
 try{renderer=new T.WebGLRenderer({canvas,antialias:true,alpha:true,powerPreference:'low-power',preserveDrawingBuffer:true});}catch{return null;}
 renderer.setPixelRatio(Math.min(devicePixelRatio||1,1.5));renderer.outputColorSpace=T.SRGBColorSpace;
 renderer.toneMapping=T.ACESFilmicToneMapping;renderer.toneMappingExposure=.95;
 renderer.shadowMap.enabled=true;renderer.shadowMap.type=T.PCFSoftShadowMap;renderer.shadowMap.autoUpdate=false;renderer.shadowMap.needsUpdate=true;
 renderer.setClearColor(0xedece4,0);
 const scene=new T.Scene(),camera=new T.PerspectiveCamera(33,1,.1,40),eye=createEye();scene.add(eye.root);
 // A photographic light rig: broad key/fill shapes become real corneal reflections.
 const studio=new T.Scene();studio.background=new T.Color('#b6b9b4');
 const boxes=[];
 for(const [p,size,color,intensity] of [
  [[-3,5,3],[4,3],'#fff6e7',5],[[4,1,1],[2,5],'#e3eef7',2.6],[[0,3,-4],[3,2],'#fff7e8',4],[[0,-4,1],[3,2],'#b5b2a6',.4]
 ]){
  const mesh=new T.Mesh(new T.PlaneGeometry(...size),new T.MeshBasicMaterial({color:new T.Color(color).multiplyScalar(intensity),side:T.DoubleSide}));mesh.position.set(...p);mesh.lookAt(0,0,0);studio.add(mesh);boxes.push(mesh);
 }
 const pmrem=new T.PMREMGenerator(renderer),environment=pmrem.fromScene(studio,.035,.1,30);scene.environment=environment.texture;scene.environmentIntensity=.5;
 pmrem.dispose();for(const b of boxes){b.geometry.dispose();b.material.dispose();}
 scene.add(new T.HemisphereLight('#fff6e7','#a5ab9e',.6));
 const key=new T.DirectionalLight('#fff0dd',3.1);key.position.set(-2.8,5,4);key.castShadow=true;key.shadow.mapSize.set(1024,1024);key.shadow.camera.left=-2;key.shadow.camera.right=2;key.shadow.camera.top=2;key.shadow.camera.bottom=-2;key.shadow.camera.near=.5;key.shadow.camera.far=12;key.shadow.normalBias=.016;key.shadow.bias=-.0002;key.shadow.radius=4;scene.add(key);
 const fill=new T.DirectionalLight('#d7e6ed',.65);fill.position.set(4,1,3);scene.add(fill);
 const rim=new T.DirectionalLight('#f4e7cf',1.8);rim.position.set(0,2,-4);scene.add(rim);
 // Contact shadow is a soft studio grounding cue, not a biological structure.
 const shadowCanvas=document.createElement('canvas');shadowCanvas.width=shadowCanvas.height=128;
 const ctx=shadowCanvas.getContext('2d'),grad=ctx.createRadialGradient(64,64,5,64,64,64);grad.addColorStop(0,'rgba(56,62,56,.24)');grad.addColorStop(.35,'rgba(56,62,56,.12)');grad.addColorStop(1,'rgba(56,62,56,0)');ctx.fillStyle=grad;ctx.fillRect(0,0,128,128);
 const shadowTexture=new T.CanvasTexture(shadowCanvas),shadow=new T.Mesh(new T.PlaneGeometry(3.5,3.5),new T.MeshBasicMaterial({map:shadowTexture,transparent:true,depthWrite:false}));shadow.rotation.x=-Math.PI/2;shadow.position.y=-1.19;scene.add(shadow);
 let yaw=.77,pitch=.28,distance=4.5,mode='cut',selected='wall',visible=true,initializing=true,raf=0,lost=false,disposed=false,drag=null,moved=false;
 const target=new T.Vector3(0,-.02,-.10),pin=root.querySelector('.eye-pin'),detail=root.querySelector('.eye-detail'),title=root.querySelector('.eye-part-title');
 const clamp=T.MathUtils.clamp;
 let stageWidth=0,stageHeight=0;
 function resize(){request();}
 function draw(){
  raf=0;if(lost||disposed||!visible)return;
  const box=stage.getBoundingClientRect();if(!box.width||!box.height)return;
  // Defer buffer resizing until a visible draw, so offscreen resizing never erases the last frame.
  if(box.width!==stageWidth||box.height!==stageHeight){stageWidth=box.width;stageHeight=box.height;renderer.setSize(stageWidth,stageHeight,false);camera.aspect=stageWidth/stageHeight;camera.updateProjectionMatrix();}
  const aspect=camera.aspect,framing=aspect<1?1/Math.sqrt(aspect):1;
  camera.position.set(target.x+distance*framing*Math.sin(yaw)*Math.cos(pitch),target.y+distance*framing*Math.sin(pitch),target.z+distance*framing*Math.cos(yaw)*Math.cos(pitch));camera.lookAt(target);camera.updateMatrixWorld();
  renderer.render(scene,camera);
  root.dataset.yaw=yaw.toFixed(3);root.dataset.zoom=(4.5/distance).toFixed(2);root.dataset.mode=mode;
  const point=eye.anchor(selected),projected=point.clone().project(camera),x=(projected.x*.5+.5)*box.width,y=(-projected.y*.5+.5)*box.height;
  // Never label an occluded interior landmark through the opaque globe.
  const ray=new T.Raycaster(camera.position,point.clone().sub(camera.position).normalize());
  const hit=ray.intersectObjects(eye.root.children.filter(g=>g.visible),true).find(h=>h.object.visible&&h.object.parent.visible&&!h.object.material.transparent&&h.object.material.transmission<.5);
  const blocked=hit&&hit.distance<camera.position.distanceTo(point)-.07;
  pin.hidden=blocked||projected.z>1||x<18||x>box.width-18||y<18||y>box.height-18;
  pin.style.left=`${x}px`;pin.style.top=`${y}px`;
  root.dataset.frames=String(Number(root.dataset.frames||0)+1);
 }
 function request(){if(!initializing&&!raf&&!lost&&!disposed&&visible)raf=requestAnimationFrame(draw);}
 function select(id,focus=false){
  if(focus){if(id==='nerve'){yaw=-2.45;pitch=.18;}else if(mode==='whole'){yaw=.08;pitch=.04;}else{yaw=mode==='layers'?1.05:.77;pitch=.28;}}
  selected=id;const i=eyeParts.findIndex(p=>p.id===id),part=eyeParts[i];
  root.querySelectorAll('[data-eye-part]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.eyePart===id)));
  title.innerHTML=`<b>${String(i+1).padStart(2,'0')}</b><span>${pair(part.en,part.zh)}</span>`;
  detail.innerHTML=pair(part.detailEn,part.detailZh);pin.textContent=String(i+1);pin.style.setProperty('--eye-part',part.color);
  root.dataset.part=id;request();
 }
 function setMode(next){
  mode=next;eye.setMode(mode);renderer.shadowMap.needsUpdate=true;root.querySelectorAll('[data-eye-mode]').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.eyeMode===mode)));
  if(mode==='whole'){yaw=.08;pitch=.04;select('iris');}else{yaw=.77;pitch=.28;selected='wall';if(mode==='layers'){yaw=1.05;pitch=.24;}}
  distance=mode==='layers'?5.2:4.5;select(selected);request();
 }
 function reset(){setMode('cut');select('wall');}
 root.querySelectorAll('[data-eye-mode]').forEach(button=>button.addEventListener('click',()=>setMode(button.dataset.eyeMode)));
 root.querySelectorAll('[data-eye-part]').forEach(button=>button.addEventListener('click',()=>{if(mode==='whole'&&['lens','ciliary','wall'].includes(button.dataset.eyePart))setMode('cut');select(button.dataset.eyePart,true);}));
 root.querySelector('[data-eye-reset]').addEventListener('click',reset);
 root.querySelector('[data-eye-in]').addEventListener('click',()=>{distance=clamp(distance/1.18,2.5,7);request();});
 root.querySelector('[data-eye-out]').addEventListener('click',()=>{distance=clamp(distance*1.18,2.5,7);request();});
 canvas.addEventListener('pointerdown',event=>{
  if(event.pointerType!=='mouse'&&event.pointerType!=='pen')return;
  if(event.button!==0)return;drag={x:event.clientX,y:event.clientY,yaw,pitch};moved=false;canvas.setPointerCapture(event.pointerId);canvas.focus({preventScroll:true});
 });
 canvas.addEventListener('pointermove',event=>{if(!drag)return;const dx=event.clientX-drag.x,dy=event.clientY-drag.y;moved ||= Math.hypot(dx,dy)>4;yaw=drag.yaw-dx*.006;pitch=clamp(drag.pitch+dy*.005,-1.1,1.1);request();});
 canvas.addEventListener('pointerup',event=>{if(drag&&!moved){const box=canvas.getBoundingClientRect(),xy=new T.Vector2((event.clientX-box.left)/box.width*2-1,-(event.clientY-box.top)/box.height*2+1),ray=new T.Raycaster();ray.setFromCamera(xy,camera);const hit=ray.intersectObjects(eye.root.children.filter(g=>g.visible),true).find(h=>h.object.visible&&h.object.parent.visible);if(hit?.object.userData.part)select(hit.object.userData.part);}drag=null;});
 canvas.addEventListener('pointercancel',()=>{drag=null;});
 // One-finger vertical scrolling remains browser-native; horizontal swipes orbit.
 let touch=null;
 canvas.addEventListener('touchstart',event=>{if(event.touches.length!==1){touch=null;return;}const t=event.touches[0];touch={x:t.clientX,y:t.clientY,yaw,axis:null};},{passive:true});
 canvas.addEventListener('touchmove',event=>{if(!touch||event.touches.length!==1)return;const t=event.touches[0],dx=t.clientX-touch.x,dy=t.clientY-touch.y;if(!touch.axis&&Math.hypot(dx,dy)>9)touch.axis=Math.abs(dx)>Math.abs(dy)*1.2?'x':'y';if(touch.axis==='x'){event.preventDefault();yaw=touch.yaw-dx*.008;request();}},{passive:false});
 canvas.addEventListener('touchend',()=>{touch=null;},{passive:true});
 canvas.addEventListener('keydown',event=>{
  const actions={ArrowLeft:()=>yaw-=.12,ArrowRight:()=>yaw+=.12,ArrowUp:()=>pitch=clamp(pitch+.1,-1.1,1.1),ArrowDown:()=>pitch=clamp(pitch-.1,-1.1,1.1),'+':()=>distance=clamp(distance/1.18,2.5,7),'=':()=>distance=clamp(distance/1.18,2.5,7),'-':()=>distance=clamp(distance*1.18,2.5,7),Home:reset};
  if(actions[event.key]){event.preventDefault();actions[event.key]();request();}
 });
 const resizeObserver=new ResizeObserver(resize);resizeObserver.observe(stage);
 const visibility=new IntersectionObserver(([entry])=>{visible=entry.isIntersecting&&!document.hidden;if(visible)request();else{cancelAnimationFrame(raf);raf=0;}});visibility.observe(stage);
 const onVisibility=()=>{visible=!document.hidden&&stage.getBoundingClientRect().bottom>0&&stage.getBoundingClientRect().top<innerHeight;if(visible)request();};document.addEventListener('visibilitychange',onVisibility);
 canvas.addEventListener('webglcontextlost',event=>{event.preventDefault();lost=true;cancelAnimationFrame(raf);raf=0;root.dataset.eyeState='fallback';canvas.hidden=true;poster.hidden=false;pin.hidden=true;});
 canvas.addEventListener('webglcontextrestored',()=>{lost=false;renderer.shadowMap.needsUpdate=true;root.dataset.eyeState='ready';canvas.hidden=false;poster.hidden=true;resize();request();});
 select('wall');resize();if(renderer.extensions.has('KHR_parallel_shader_compile'))await renderer.compileAsync(scene,camera);else renderer.compile(scene,camera);initializing=false;visible=true;draw();
 canvas.hidden=false;poster.hidden=true;root.dataset.eyeState='ready';
 return {renderer,eye,camera,scene,draw:()=>{visible=true;draw();},setMode,select,dispose(){disposed=true;cancelAnimationFrame(raf);resizeObserver.disconnect();visibility.disconnect();document.removeEventListener('visibilitychange',onVisibility);eye.dispose();environment.dispose();shadowTexture.dispose();shadow.geometry.dispose();shadow.material.dispose();key.shadow.map?.dispose();renderer.dispose();}};
}
