import * as T from './three.module.min.js';

// Original educational geometry. Coordinates are normalized to the globe radius.
// Coat thickness and texture contrast are enhanced for legibility, not clinical use.
const TAU=Math.PI*2,START=.54,OPEN=Math.PI*.57,END=TAU-.12;
const clamp=T.MathUtils.clamp;
const vec=a=>new T.Vector3(...a);
export const eyeParts=[
 {id:'cornea',en:'Cornea',zh:'角膜',color:'#9bbeb9',anchor:[.04,-.12,1.062],detailEn:'The transparent anterior surface refracts incoming light. Its curvature and the air–tear interface contribute most of the eye’s optical power.',detailZh:'前部透明曲面折射入射光。角膜曲率与空气—泪膜界面贡献了眼球大部分屈光力。'},
 {id:'iris',en:'Iris / pupil',zh:'虹膜／瞳孔',color:'#6e9b8e',anchor:[-.19,-.28,.779],detailEn:'Radial stromal fibers surround the pupil. The pupil is an opening, not a black tissue layer; iris muscles regulate its diameter.',detailZh:'放射状基质纤维围绕瞳孔。瞳孔是开口，不是一层黑色组织；虹膜肌肉调节开口直径。'},
 {id:'lens',en:'Lens',zh:'晶状体',color:'#cab67d',anchor:[.13,.08,.662],detailEn:'The biconvex lens lies behind the iris. Zonular fibers attach near its equator; accommodation changes its curvature.',detailZh:'双凸晶状体位于虹膜后方。悬韧带在赤道部附近附着，调节通过改变晶状体曲率完成。'},
 {id:'ciliary',en:'Ciliary body',zh:'睫状体',color:'#b77f78',anchor:[-.30,-.41,.57],detailEn:'Ciliary processes form a folded ring behind the iris. Fine zonular fibers connect this apparatus to the lens capsule.',detailZh:'睫状突在虹膜后形成褶皱环，细的悬韧带将睫状结构与晶状体囊连接。'},
 {id:'wall',en:'Eye wall',zh:'眼球壁',color:'#c8896f',anchor:[.53,-.065,-.741],detailEn:'From the cavity outward: neural retina, a thin pigment epithelium, vascular choroid and fibrous sclera. The clear vitreous normally fills this cavity; it is omitted here to expose the tissue.',detailZh:'由腔内向外依次是神经视网膜、薄层色素上皮、富含血管的脉络膜与纤维性巩膜。透明玻璃体通常充填此腔，这里省去以显示组织。'},
 {id:'nerve',en:'Optic nerve',zh:'视神经',color:'#c4ad83',anchor:[-.62,.05,-1.51],detailEn:'Ganglion-cell axons converge at the optic disc and continue into the optic nerve. The fascicles and central vessels are enlarged in this model.',detailZh:'神经节细胞轴突在视盘汇合并延续为视神经。模型放大显示了轴突束及中央血管。'}
];

function random(seed=271828){return()=>{seed=(Math.imul(seed,1664525)+1013904223)>>>0;return seed/4294967296;};}
function canvasTexture(canvas,color=true){const t=new T.CanvasTexture(canvas);if(color)t.colorSpace=T.SRGBColorSpace;t.anisotropy=4;return t;}
function tissueTexture(rgb,seed=31){
 const rand=random(seed),c=document.createElement('canvas');c.width=512;c.height=256;
 const ctx=c.getContext('2d'),im=ctx.createImageData(c.width,c.height);
 for(let y=0;y<c.height;y++)for(let x=0;x<c.width;x++){
  const n=(rand()-.5)*9+4*Math.sin(x*.079+Math.sin(y*.092))*Math.sin(y*.065),i=(y*c.width+x)*4;
  for(let k=0;k<3;k++)im.data[i+k]=clamp(rgb[k]+n,0,255);im.data[i+3]=255;
 }
 ctx.putImageData(im,0,0);return canvasTexture(c);
}
function irisTexture(){
 const size=1024,c=document.createElement('canvas');c.width=c.height=size;
 const ctx=c.getContext('2d'),im=ctx.createImageData(size,size),rand=random(691);
 const fibers=Array.from({length:720},()=>rand());
 for(let y=0;y<size;y++)for(let x=0;x<size;x++){
  const X=(x-size/2)/(size/2),Y=(y-size/2)/(size/2),r=Math.hypot(X,Y),a=Math.atan2(Y,X),i=(y*size+x)*4;
  const v=(a+Math.PI)/TAU*720,lo=Math.floor(v)%720;
  const stria=.45*Math.sin(289*a+Math.sin(r*52+a*9)*.9)+.25*Math.sin(641*a+r*32)+.2*Math.sin(127*a-r*18);
  const warp=r+.015*Math.sin(a*37)+.009*Math.sin(a*87),crypt=Math.pow(Math.max(0,Math.sin(a*29+.8*Math.sin(a*7))),8)*Math.exp(-(((warp-.59)/.095)**2));
  const collarette=Math.exp(-(((warp-.55)/.02)**2)),radial=fibers[lo]*(.4+.6*Math.sin(r*7)**2);
  let light=1+stria*.24+radial*.32-crypt*.45+collarette*.18;
  light*=.68+.32*clamp((1-r)/.07,0,1);light*=.5+.5*clamp((r-.30)/.055,0,1);
  const amber=Math.exp(-(((r-.43)/.18)**2)),rgb=[62+amber*44,90+amber*7,78-amber*30];
  for(let k=0;k<3;k++)im.data[i+k]=r<.303?12:clamp(rgb[k]*light+(rand()-.5)*4,0,255);im.data[i+3]=255;
 }
 ctx.putImageData(im,0,0);return canvasTexture(c);
}

function surface(fn,nu=56,nv=180,flip=false){
 const pos=[],uv=[],ix=[];
 for(let i=0;i<=nu;i++)for(let j=0;j<=nv;j++){const u=i/nu,v=j/nv,q=fn(u,v);pos.push(...q.p);uv.push(...(q.uv||[v,u]));}
 for(let i=0;i<nu;i++)for(let j=0;j<nv;j++){
  const a=i*(nv+1)+j,b=a+nv+1;if(flip)ix.push(a,a+1,b,b,a+1,b+1);else ix.push(a,b,a+1,b,b+1,a+1);
 }
 const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(pos,3));g.setAttribute('uv',new T.Float32BufferAttribute(uv,2));g.setIndex(ix);g.computeVertexNormals();return g;
}
function shell(radius,start,phi,length,inside=false){
 return surface((u,v)=>{const th=start+(Math.PI-start)*u,p=phi+length*v;
  const distort=1+.012*Math.cos(th*2),r=radius*distort;
  return{p:[r*Math.sin(th)*Math.cos(p),r*Math.sin(th)*Math.sin(p),radius*1.015*Math.cos(th)],uv:[p/TAU,th/Math.PI]};
 },56,128,inside);
}
function merge(geometries){
 const positions=[],normals=[],uvs=[];
 for(const input of geometries){const g=input.index?input.toNonIndexed():input;positions.push(...g.attributes.position.array);normals.push(...g.attributes.normal.array);if(g.attributes.uv)uvs.push(...g.attributes.uv.array);if(g!==input)g.dispose();input.dispose();}
 const out=new T.BufferGeometry();out.setAttribute('position',new T.Float32BufferAttribute(positions,3));out.setAttribute('normal',new T.Float32BufferAttribute(normals,3));if(uvs.length)out.setAttribute('uv',new T.Float32BufferAttribute(uvs,2));return out;
}
function tube(points,r=.004,steps=36,segments=7){const c=new T.CatmullRomCurve3(points.map(vec));return new T.TubeGeometry(c,steps,r,segments,false);}
function taperedTube(points,r,steps=36,segments=7){
 const curve=new T.CatmullRomCurve3(points.map(vec)),g=new T.TubeGeometry(curve,steps,r,segments,false),p=g.attributes.position;
 for(let i=0;i<=steps;i++){const c=curve.getPointAt(i/steps),f=.17+.83*(1-i/steps)**.6;for(let j=0;j<=segments;j++){const k=i*(segments+1)+j;p.setXYZ(k,c.x+(p.getX(k)-c.x)*f,c.y+(p.getY(k)-c.y)*f,c.z+(p.getZ(k)-c.z)*f);}}
 g.computeVertexNormals();return g;
}
function sampled(fn,n=20){return Array.from({length:n+1},(_,i)=>fn(i/n));}

export function createEye(){
 const root=new T.Group(),groups={},resources=new Set();
 const iris=irisTexture(),sclera=tissueTexture([223,219,203]),retina=tissueTexture([203,139,111],89),choroid=tissueTexture([116,62,55],43);
 for(const t of [iris,sclera,retina,choroid])resources.add(t);
 function material(color,extra={}){const m=new T.MeshPhysicalMaterial({color,roughness:.53,metalness:0,side:T.DoubleSide,envMapIntensity:.65,...extra});resources.add(m);return m;}
 const mats={
  sclera:material('#fffaf0',{map:sclera,bumpMap:sclera,bumpScale:.009,roughness:.44,clearcoat:.3,clearcoatRoughness:.35}),
  scleraCut:material('#d8ccb0',{roughness:.74}),
  choroid:material('#f9eeeb',{map:choroid,bumpMap:choroid,bumpScale:.008,roughness:.63}),
  pigment:material('#5e3d37',{roughness:.82}),
  retina:material('#fff3e3',{map:retina,bumpMap:retina,bumpScale:.002,roughness:.48,clearcoat:.16}),
  retinaCut:material('#e7be8a',{roughness:.7}),
  cornea:material('#080e10',{transparent:true,opacity:.18,depthWrite:false,side:T.FrontSide,ior:1.376,roughness:.08,envMapIntensity:2,clearcoat:1,clearcoatRoughness:.035}),
  iris:material('#ffffff',{map:iris,bumpMap:iris,bumpScale:.0015,roughness:.68,clearcoat:.12,clearcoatRoughness:.6}),
  pupil:material('#372d25',{roughness:.84}),
  lens:material('#b8a87e',{side:T.FrontSide,transparent:true,opacity:.23,depthWrite:false,ior:1.4,roughness:.15,envMapIntensity:1,clearcoat:.9,clearcoatRoughness:.1}),
  lensCut:material('#ddd4b7',{transparent:true,opacity:.62,depthWrite:false,roughness:.25,ior:1.4}),
  lamella:material('#bcb895',{roughness:.66}),
  ciliary:material('#995d57',{roughness:.59,clearcoat:.16}),
  ciliaryRidge:material('#b77568',{roughness:.48,clearcoat:.15}),
  zonule:material('#d3c7a8',{roughness:.45}),
  vein:material('#644751',{roughness:.43,clearcoat:.2}),
  artery:material('#a4534b',{roughness:.4,clearcoat:.2}),
  superficial:material('#bd9791',{roughness:.65}),
  nerve:material('#d1b98c',{roughness:.6}),
  nerveFiber:material('#efe0b7',{roughness:.5}),
  disc:material('#e4bc8a',{roughness:.6})
 };
 function mesh(group,g,m,part){const o=new T.Mesh(g,m);o.userData.part=part;o.castShadow=!m.transmission&&!m.transparent;o.receiveShadow=true;group.add(o);resources.add(g);return o;}
 function makeGroup(name,mode,part,offset){const g=new T.Group();g.name=name;g.userData={mode,part,offset:vec(offset)};root.add(g);groups[name]=g;return g;}
 // Build an intact globe and a cutaway independently; mode changes swap visibility.
 for(const cut of [true,false]){
  const mode=cut?'cut':'whole',phi=cut?OPEN:0,len=cut?END-OPEN:TAU;
  const wall=makeGroup(mode+'-wall',mode,'wall',[0,0,-.26]);
  // In the intact globe, little light returns through the pupil from the pigmented interior.
  const retinalMaterial=cut?mats.retina:material('#080404',{roughness:1,specularIntensity:0,envMapIntensity:0});
  for(const [outer,inner,mat,cap,start] of [[1,.958,mats.sclera,mats.scleraCut,START],[.955,.930,mats.choroid,mats.choroid,.68],[.929,.923,mats.pigment,mats.pigment,.95],[.922,.907,retinalMaterial,mats.retinaCut,.95]]){
   mesh(wall,shell(outer,start,phi,len),mat,'wall');
   // Interior normals face the cavity. All coats have watertight cut-edge strips.
   mesh(wall,shell(inner,start,phi,len,true),mat,'wall');
   for(const angle of cut?[phi,phi+len]:[]){
    const edge=surface((u,v)=>{const th=start+(Math.PI-start)*v,r=inner+(outer-inner)*u,d=1+.012*Math.cos(th*2);return{p:[r*d*Math.sin(th)*Math.cos(angle),r*d*Math.sin(th)*Math.sin(angle),r*1.015*Math.cos(th)]};},3,100);
    mesh(wall,edge,cap,'wall');
   }
   mesh(wall,surface((u,v)=>{const a=phi+len*v,r=inner+(outer-inner)*u;return{p:[r*Math.sin(start)*Math.cos(a),r*Math.sin(start)*Math.sin(a),r*1.015*Math.cos(start)]};},2,160),cap,'wall');
  }
  const front=makeGroup(mode+'-cornea',mode,'cornea',[0,0,.69]);
  mesh(front,surface((u,v)=>{const th=Math.PI/2*u,a=phi+len*v,r=.519*Math.sin(th);return{p:[r*Math.cos(a),r*Math.sin(a),.867+.23*Math.cos(th)],uv:[.5+r*Math.cos(a),.5+r*Math.sin(a)]};},64,160),mats.cornea,'cornea');
  // The limbus marks the transition from sclera to cornea without a thick frame.
  mesh(front,tube(sampled(t=>{const a=phi+len*t;return[.521*Math.cos(a),.521*Math.sin(a),.864];},130),.0045,130),mats.sclera,'cornea');
  const ir=makeGroup(mode+'-iris',mode,'iris',[0,0,.46]);
  mesh(ir,surface((u,v)=>{const a=phi+len*v,r=.163+u*.365,z=.771+.014*Math.sin(u*Math.PI)+.003*Math.sin(a*89)*Math.sin(u*Math.PI);return{p:[r*Math.cos(a),r*Math.sin(a),z],uv:[.5+.5*r/.53*Math.cos(a),.5+.5*r/.53*Math.sin(a)]};},48,256),mats.iris,'iris');
  mesh(ir,tube(sampled(t=>{const a=phi+len*t;return[.163*Math.cos(a),.163*Math.sin(a),.774];},150),.005,150),mats.pupil,'iris');
  // Posterior pigment surface gives the thin iris a physical edge.
  mesh(ir,surface((u,v)=>{const a=phi+len*v,r=.163+u*.365;return{p:[r*Math.cos(a),r*Math.sin(a),.756]};},12,160,true),mats.pupil,'iris');
  const lens=makeGroup(mode+'-lens',mode,'lens',[0,0,.19]);
  const lensMaterial=cut?mats.lens:mats.lens.clone();
  if(!cut){lensMaterial.opacity=.04;lensMaterial.specularIntensity=.15;lensMaterial.envMapIntensity=.15;lensMaterial.clearcoat=.08;resources.add(lensMaterial);}
  mesh(lens,surface((u,v)=>{const th=u*Math.PI,a=phi+len*v,r=.396*Math.sin(th),z=.553+(th<Math.PI/2?.13:.185)*Math.cos(th);return{p:[r*Math.cos(a),r*Math.sin(a),z]};},64,128),lensMaterial,'lens');
  if(cut){
   for(const a of [phi,phi+len]){
    mesh(lens,surface((u,v)=>{const th=v*Math.PI,r=.395*u*Math.sin(th),z=.553+(th<Math.PI/2?.13:.185)*u*Math.cos(th);return{p:[r*Math.cos(a),r*Math.sin(a),z]};},28,100),mats.lensCut,'lens');
    const lam=[];for(let j=1;j<15;j++){const f=j/15;lam.push(tube(sampled(t=>{const th=t*Math.PI,r=.395*f*Math.sin(th);return[r*Math.cos(a+.003),r*Math.sin(a+.003),.553+(th<Math.PI/2?.13:.185)*f*Math.cos(th)];},46),.0008,46,5));}mesh(lens,merge(lam),mats.lamella,'lens');
   }
  }
  const ciliary=makeGroup(mode+'-ciliary',mode,'ciliary',[0,0,.035]);
  mesh(ciliary,surface((u,v)=>{const a=phi+len*v,r=.48+.24*u,z=.46+.13*Math.sin(u*Math.PI);return{p:[r*Math.cos(a),r*Math.sin(a),z]};},20,256),mats.ciliary,'ciliary');
  const ridges=[],fibers=[];
  for(let k=0;k<64;k++){
   const a=TAU*k/64;if(cut&&(a<phi||a>phi+len))continue;
   const g=new T.SphereGeometry(1,12,16);g.scale(.016,.020,.065);g.rotateY(.55);g.rotateZ(a);g.translate(.536*Math.cos(a),.536*Math.sin(a),.566);ridges.push(g);
   for(const d of [-.025,0,.025])fibers.push(tube([[.53*Math.cos(a),.53*Math.sin(a),.585],[.46*Math.cos(a+.012),.46*Math.sin(a+.012),.56+d*.5],[.393*Math.cos(a+.02),.393*Math.sin(a+.02),.55+d]],.00115,8,5));
  }
  mesh(ciliary,merge(ridges),mats.ciliaryRidge,'ciliary');mesh(ciliary,merge(fibers),mats.zonule,'ciliary');
  // Subtle superficial vessels follow the surface, never through the corneal window.
  const surfaceVessels=[];
  for(let k=0;k<17;k++){
   const a=.12+k*TAU/17;if(cut&&(a<phi+.08||a>phi+len-.08))continue;
   const on=(t,shift=0)=>{const th=.57+t*.75,p=a+shift+.025*Math.sin(t*12+k),r=1.002*(1+.012*Math.cos(th*2));return[r*Math.sin(th)*Math.cos(p),r*Math.sin(th)*Math.sin(p),1.017*Math.cos(th)];};
   surfaceVessels.push(taperedTube(sampled(t=>on(t),22),.0019,22,5));
   for(const tt of [.3,.6])surfaceVessels.push(taperedTube(sampled(t=>on(tt+t*.2,t*.055),9),.0012,9,5));
  }
  mesh(wall,merge(surfaceVessels),mats.superficial,'wall');
  // Vessels stay on the inner retina; four main arcades issue from the disc.
  const onRetina=(x,y)=>{const r=.905,q=Math.hypot(x,y),sc=q>r*.92?r*.92/q:1;return[x*sc,y*sc,-Math.sqrt(r*r-q*q*sc*sc)*1.015];};
  for(const [mat,shift] of [[mats.artery,0],[mats.vein,.024]]){
   const vessels=[];
   for(const sx of [-1,1])for(const sy of [-1,1]){
    const fn=t=>onRetina(-.24+sx*.53*t,.02+sy*.58*Math.sin(t*1.4)+shift);
    function append(points,r,n){
     // Exclude the removed sector; splitting prevents triangles across its gap.
     let strip=[];for(const pt of points){const a=(Math.atan2(pt[1],pt[0])+TAU)%TAU;if(!cut||a>=phi&&a<=phi+len)strip.push(pt);else{if(strip.length>2)vessels.push(taperedTube(strip,r,n,7));strip=[];}}if(strip.length>2)vessels.push(taperedTube(strip,r,n,7));
    }
    append(sampled(fn,48),.008,48);
    for(const at of [.25,.45,.67,.84]){const p=fn(at);append(sampled(t=>onRetina(p[0]+sx*.22*t,p[1]+sy*.17*t+.035*Math.sin(t*3.5)),18),.0037,18);}
   }
   mesh(wall,merge(vessels),cut?mat:retinalMaterial,'wall');
  }
  const nerve=makeGroup(mode+'-nerve',mode,'nerve',[-.035,0,-.28]);
  const nervePoints=[[-.24,.02,-.88],[-.34,.025,-1.16],[-.60,.05,-1.52]],curve=new T.CatmullRomCurve3(nervePoints.map(vec));
  mesh(nerve,new T.TubeGeometry(curve,48,.133,48,false),mats.nerve,'nerve');
  const disc=new T.SphereGeometry(1,36,24);disc.scale(.097,.116,.021);disc.translate(-.24,.02,-.869);mesh(wall,disc,cut?mats.disc:retinalMaterial,'nerve');
  const fascicles=[],end=curve.getPoint(1),axis=curve.getTangent(1),q=new T.Quaternion().setFromUnitVectors(new T.Vector3(0,0,1),axis);
  const endFace=new T.CircleGeometry(.132,64);endFace.applyQuaternion(q);endFace.translate(end.x,end.y,end.z);mesh(nerve,endFace,mats.nerve,'nerve');
  for(let i=0;i<51;i++){
   const r=.112*Math.sqrt((i+.5)/51),a=i*2.39996,offset=new T.Vector3(r*Math.cos(a),r*Math.sin(a),.001).applyQuaternion(q);
   const g=new T.SphereGeometry(1,9,7);g.scale(.011,.014,.003);g.applyQuaternion(q);g.translate(end.x+offset.x,end.y+offset.y,end.z+offset.z);fascicles.push(g);
  }
  mesh(nerve,merge(fascicles),mats.nerveFiber,'nerve');
  for(const [mat,dx] of [[mats.artery,-.013],[mats.vein,.013]])mesh(nerve,tube(nervePoints.map(p=>[p[0]+dx,p[1],p[2]-.009]),.009,40,10),mat,'nerve');
 }
 let mode='cut';
 function setMode(next){mode=next;for(const g of Object.values(groups)){g.visible=g.userData.mode===(mode==='whole'?'whole':'cut');g.position.copy(g.userData.offset).multiplyScalar(mode==='layers'?1:0);}}
 setMode(mode);
 return{root,groups,parts:eyeParts,setMode,anchor(id){const spec=eyeParts.find(x=>x.id===id),g=groups[(mode==='whole'?'whole':'cut')+'-'+id];return vec(spec.anchor).add(g?.position||new T.Vector3());},dispose(){for(const r of resources)r.dispose();}};
}
