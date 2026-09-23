self.onmessage = ({data}) => {
 try {
  const size=Number(data.size);if(![65536,262144,1048576].includes(size))throw Error('Invalid input size');
  const values=Float64Array.from({length:size},(_,i)=>i%97),warmup=12,rounds=15,batch=16;
  function plain(){let total=0;for(let i=0;i<values.length;i++)total+=values[i];return total;}
  function unrolled(){let a=0,b=0,c=0,d=0,i=0;for(;i+3<values.length;i+=4){a+=values[i];b+=values[i+1];c+=values[i+2];d+=values[i+3];}let total=a+b+c+d;for(;i<values.length;i++)total+=values[i];return total;}
  const expected=plain();if(unrolled()!==expected)throw Error('Output verification failed');
  for(let i=0;i<warmup;i++){plain();unrolled();}
  const samples={plain:[],unrolled:[]};let checksum=0;
  function measure(fn){const start=performance.now();let sum=0;for(let i=0;i<batch;i++)sum+=fn();const elapsed=performance.now()-start;if(sum!==expected*batch)throw Error('Output verification failed during timing');checksum+=sum;return elapsed/batch;}
  for(let i=0;i<rounds;i++)for(const name of i%2?['unrolled','plain']:['plain','unrolled'])samples[name].push(measure(name==='plain'?plain:unrolled));
  const summarize=a=>{const s=[...a].sort((x,y)=>x-y),quantile=q=>{const p=(s.length-1)*q,n=Math.floor(p);return s[n]+(s[Math.ceil(p)]-s[n])*(p-n);};return {median:quantile(.5),q1:quantile(.25),q3:quantile(.75),min:s[0],max:s.at(-1)};};
  self.postMessage({timestamp:new Date().toISOString(),size,warmup,rounds,batch,units:'milliseconds per summation',expected,checksum,samples,summary:{plain:summarize(samples.plain),unrolled:summarize(samples.unrolled)},notes:'Web Worker; alternating order; allocation excluded; not a statistical significance test.'});
 }catch(error){self.postMessage({error:error.message});}
};
