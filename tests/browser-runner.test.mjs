import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {spawnSync} from 'node:child_process';
import {runSuites} from '../scripts/check-browsers.mjs';

test('browser runner retains diagnostics and continues after failures and process-group timeouts',{timeout:10000,skip:process.platform==='win32'},async()=>{
 const directory=await fs.mkdtemp(path.join(os.tmpdir(),'feng-runner-'));
 try{
  const sources={
   fail:'console.error("failure detail");process.exit(2);',
   hang:'import {spawn} from "node:child_process";spawn(process.execPath,["-e","setInterval(()=>{},1000)"],{stdio:"inherit"});setInterval(()=>{},1000);',
   pass:'console.log("later suite ran");'
  };
  const suites=[];
  for(const [name,source]of Object.entries(sources)){const file=path.join(directory,name+'.mjs');await fs.writeFile(file,source);suites.push({name,file});}
  const reports=path.join(directory,'reports');
  const results=await runSuites(suites,{directory:reports,timeoutMs:1000,output:()=>{}});
  assert.deepEqual(results.map(r=>r.status),['failed','timeout','passed']);
  assert.equal(results[0].code,2);
  assert.match(await fs.readFile(path.join(reports,'fail.log'),'utf8'),/failure detail/);
  assert.match(await fs.readFile(path.join(reports,'hang.log'),'utf8'),/terminating its process group/);
  assert.match(await fs.readFile(path.join(reports,'pass.log'),'utf8'),/later suite ran/);
  assert.deepEqual(JSON.parse(await fs.readFile(path.join(reports,'results.json'),'utf8')),results);
 }finally{await fs.rm(directory,{recursive:true,force:true});}
});

test('browser CLI discovers suites and exits nonzero even when the last suite passes',async()=>{
 const directory=await fs.mkdtemp(path.join(os.tmpdir(),'feng-runner-cli-'));
 try{
  await fs.mkdir(path.join(directory,'tests'));
  await fs.writeFile(path.join(directory,'tests/fail.mjs'),'process.exit(1);');
  await fs.writeFile(path.join(directory,'tests/pass.mjs'),'console.log("last suite passed");');
  await fs.writeFile(path.join(directory,'package.json'),JSON.stringify({scripts:{'test:browsers':'ignored','test:fail':'node tests/fail.mjs','test:pass':'node tests/pass.mjs'}}));
  const run=spawnSync(process.execPath,[path.resolve('scripts/check-browsers.mjs')],{cwd:directory,encoding:'utf8',timeout:5000,env:{...process.env,GITHUB_STEP_SUMMARY:''}});
  assert.equal(run.status,1,run.stderr);
  assert.match(run.stdout,/last suite passed/);
  const results=JSON.parse(await fs.readFile(path.join(directory,'test-results/browser/results.json'),'utf8'));
  assert.deepEqual(results.map(r=>r.status),['failed','passed']);
 }finally{await fs.rm(directory,{recursive:true,force:true});}
});

test('cancelling a run terminates the active suite and skips the remaining suites',{skip:process.platform==='win32'},async()=>{
 const directory=await fs.mkdtemp(path.join(os.tmpdir(),'feng-runner-cancel-'));
 try{
  const hang=path.join(directory,'hang.mjs'),wrapper=path.join(directory,'wrapper.mjs');
  await fs.writeFile(hang,'setInterval(()=>{},1000);');
  const moduleUrl=new URL('../scripts/check-browsers.mjs',import.meta.url).href;
  await fs.writeFile(wrapper,`import {runSuites} from ${JSON.stringify(moduleUrl)};
   setTimeout(()=>process.kill(process.pid,'SIGTERM'),500);
   const results=await runSuites([{name:'hang',file:${JSON.stringify(hang)}},{name:'must-not-start',file:'missing.mjs'}],{timeoutMs:4000});
   if(results.length!==1||results[0].status!=='cancelled')process.exitCode=1;`);
  const run=spawnSync(process.execPath,[wrapper],{cwd:directory,encoding:'utf8',timeout:6000});
  assert.equal(run.status,0,run.stderr);
  assert.match(run.stdout,/hang: cancelled/);
  assert.doesNotMatch(run.stdout,/must-not-start/);
 }finally{await fs.rm(directory,{recursive:true,force:true});}
});
