// Teaching models: numerical results are independent of animation timing.
export function frankWolfeStep(x,target){
 const g=x.map((v,i)=>v-target[i]),vertices=[[0,0],[1,0],[0,1]];
 const dot=(a,b)=>a.reduce((sum,v,i)=>sum+v*b[i],0);
 const s=vertices.reduce((a,b)=>dot(g,b)<dot(g,a)?b:a),d=s.map((v,i)=>v-x[i]);
 const gap=-dot(g,d),gamma=dot(d,d)?Math.max(0,Math.min(1,gap/dot(d,d))):0;
 return {g,s,d,gap,gamma,next:x.map((v,i)=>v+gamma*d[i]),objective:t=>.5*x.reduce((sum,v,i)=>sum+(v+t*d[i]-target[i])**2,0)};
}
export function cudaWrites(n,width,broken=false,blocks=Math.ceil(n/width)){
 const writes=Array(n).fill(0);let masked=0;
 for(let b=0;b<blocks;b++)for(let t=0;t<width;t++){const i=broken?t:b*width+t;if(i<n)writes[i]++;else masked++;}
 return {writes,masked,missing:writes.filter(v=>!v).length,collisions:writes.filter(v=>v>1).length};
}
export function lifCycle(J,tau,theta,ref=2){
 const spike=J>theta,tint=spike?tau*Math.log(J/(J-theta)):Infinity,horizon=spike?tint+ref:4*tau;
 return {spike,tint,horizon,ref,potential:t=>spike&&t>=tint?0:J*(1-Math.exp(-t/tau)),rate:spike?1000/horizon:0};
}
export const bandAddress=(i,j,n,p,q)=>i<j-q||i>j+p?null:{row:q+i-j,column:j,offset:j*(p+q+1)+q+i-j};
export function symmetricAddress(i,j,n){const row=Math.max(i,j),column=Math.min(i,j);return {row,column,offset:column*(2*n-column+1)/2+row-column};}
export const transposeOrder=(rows,cols)=>Array.from({length:rows*cols},(_,i)=>(i%rows)*cols+Math.floor(i/rows));
export const superposedWave=(x,t,phase,k=1.5)=>Math.sin(k*x-t)+Math.sin(k*x-t+phase);
