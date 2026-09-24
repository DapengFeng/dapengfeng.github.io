import fs from 'node:fs/promises';
import path from 'node:path';
import {spawn} from 'node:child_process';
import {fileURLToPath} from 'node:url';

// Run every suite, keeping failures visible without skipping later checks.
export async function runSuites(suites,{directory='test-results/browser',timeoutMs=300000,output=chunk=>process.stdout.write(chunk)}={}){
 if(!suites.length)throw Error('No browser suites configured');
 await fs.mkdir(directory,{recursive:true});
 const results=[];let interrupted=false;
 for(const {name,file} of suites){
  if(!/^[\w-]+$/.test(name))throw Error(`Invalid suite name: ${name}`);
  const log=await fs.open(path.join(directory,`${name}.log`),'w');
  const start=Date.now();
  output(`\n--- ${name} ---\n`);
  const result=await new Promise(resolve=>{
   let timedOut=false,error;
   const child=spawn(process.execPath,[file],{detached:process.platform!=='win32',stdio:['ignore','pipe','pipe']});
   let writes=Promise.resolve();
   const capture=chunk=>{output(chunk);writes=writes.then(()=>log.write(chunk));};
   child.stdout.on('data',capture);child.stderr.on('data',capture);
   const terminate=()=>{
    try{if(process.platform==='win32')child.kill('SIGKILL');else process.kill(-child.pid,'SIGKILL');}catch(e){if(e.code!=='ESRCH')error=e.message;}
   };
   const interrupt=()=>{interrupted=true;terminate();};
   process.once('SIGINT',interrupt);process.once('SIGTERM',interrupt);
   const timer=setTimeout(()=>{
    timedOut=true;capture(Buffer.from(`\nSuite exceeded ${timeoutMs} ms; terminating its process group.\n`));
    terminate();
   },timeoutMs);
   child.on('error',e=>{error=e.message;capture(Buffer.from(error+'\n'));});
   child.on('close',async(code,signal)=>{
    clearTimeout(timer);
    process.removeListener('SIGINT',interrupt);process.removeListener('SIGTERM',interrupt);
    await writes;await log.close();
    resolve({name,status:interrupted?'cancelled':timedOut?'timeout':code===0&&!error?'passed':'failed',code,signal,...(error?{error}:{}),durationMs:Date.now()-start});
   });
  });
  results.push(result);
  // Persist after each suite so an interrupted CI run still has useful results.
  await fs.writeFile(path.join(directory,'results.json'),JSON.stringify(results,null,2)+'\n');
  output(`${name}: ${result.status}\n`);
  if(interrupted)break;
 }
 return results;
}

if(process.argv[1]===fileURLToPath(import.meta.url)){
 const {scripts}=JSON.parse(await fs.readFile('package.json','utf8'));
 const suites=Object.entries(scripts).filter(([name])=>name.startsWith('test:')&&name!=='test:browsers').map(([name,command])=>{
  const match=/^node (tests\/[\w-]+\.mjs)$/.exec(command);
  if(!match)throw Error(`Unsupported browser suite command: ${name}`);
  return {name:name.slice(5),file:match[1]};
 });
 const directory='test-results/browser';
 await fs.rm(directory,{recursive:true,force:true});
 const results=await runSuites(suites,{directory});
 const table='| Suite | Result | Seconds |\n| --- | --- | --- |\n'+results.map(r=>`| ${r.name} | ${r.status} | ${(r.durationMs/1000).toFixed(1)} |`).join('\n')+'\n';
 console.log('\n'+table);
 if(process.env.GITHUB_STEP_SUMMARY)await fs.appendFile(process.env.GITHUB_STEP_SUMMARY,table);
 process.exitCode=results.every(r=>r.status==='passed')?0:1;
}
