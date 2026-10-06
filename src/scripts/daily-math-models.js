// The registry is shared by publication validation and the canvas renderer.
export const renderers=['fourier','lissajous','standing','heat','saddle','mobius','phyllotaxis','euler','taylor','beats','gibbs','convolution','gaussian','aliasing','bezier','cycloid','spiral','ellipse','involute','determinant','eigenvectors','svd','gradient','newton','bisection','riemann','divergence','curl','walk','clt','poisson','bayes'];
export function factorial(n){let v=1;for(let i=2;i<=n;i++)v*=i;return v;}
export function choose(n,k){let v=1;for(let i=1;i<=k;i++)v*= (n-i+1)/i;return v;}
export function taylorSine(x,n){let v=0;for(let k=0;k<n;k++)v+=(-1)**k*x**(2*k+1)/factorial(2*k+1);return v;}
export const convolution=t=>Math.max(0,1-Math.abs(t-1));
export const lerp=(a,b,t)=>a.map((v,i)=>v+(b[i]-v)*t);
export function bezier(points,t){let p=points;while(p.length>1)p=p.slice(1).map((v,i)=>lerp(p[i],v,t));return p[0];}
export const normal=x=>Math.exp(-x*x/2)/Math.sqrt(2*Math.PI);
export function uniformSumDensity(z,n){
 const scale=Math.sqrt(n/12),x=z*scale+n/2;
 if(x<=0||x>=n)return 0;
 let sum=0;for(let k=0;k<=Math.floor(x);k++)sum+=(-1)**k*choose(n,k)*(x-k)**(n-1);
 return Math.max(0,scale*sum/factorial(n-1));
}
export const poissonPMF=(k,lambda=3)=>Math.exp(-lambda)*lambda**k/factorial(k);
export const binomialPMF=(k,n,p)=>k>n?0:choose(n,k)*p**k*(1-p)**(n-k);
export const betaPDF=(x,a,b)=>x<0||x>1?0:x**(a-1)*(1-x)**(b-1)*factorial(a+b-1)/(factorial(a-1)*factorial(b-1));
export function randomWalk(seed,n=120){
 let state=seed>>>0;const points=[[0,0]];
 for(let i=0;i<n;i++){state=(Math.imul(1664525,state)+1013904223)>>>0;const a=state/4294967296*Math.PI*2,p=points.at(-1);points.push([p[0]+Math.cos(a),p[1]+Math.sin(a)]);}
 return points;
}
