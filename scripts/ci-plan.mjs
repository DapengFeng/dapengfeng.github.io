import fs from 'node:fs/promises';
import {execFileSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {discoverSuites} from './check-browsers.mjs';

const core=['browser','language','reading'];
const article=[...core,'layout','accessibility','math','code-copy','sharing','support'];
// daily-math itself checks formulas, language modes, responsive layouts and static fallback.
const daily=[...core,'daily-math'];
const posts={
 'pytorch-01-what-is-pytorch':['pytorch','reading-demos'],
 'pytorch-02-tensor-strides-storage':['tensor','reading-demos','process-demos'],
 'pytorch-03-operator-dispatch':['dispatch'],
 'pytorch-04-autograd-engine':['autograd'],
 'pytorch-05-cuda-streams-timing':['cuda-streams'],
 'pytorch-06-compile-and-codegen':['compile'],
 'human-visual-system':['vision'],
 'cuda-rust-two-tracks-blog':['cuda','process-demos'],
 'rust-vs-cpp-blog':['compiler','reading-demos'],
 'spike_notes':['compiler','process-demos'],
 'frank-wolfe-algorithm':['lessons','reading-demos','process-demos'],
 'matrix-multiplication':['lessons','reading-demos'],
 'fast-matrix-vector-products':['lessons'],
 'band-storage-gaxpy':['lessons','process-demos'],
 'symmetric-storage-gaxpy':['lessons','process-demos'],
 'waves-and-phase':['process-demos'],
 'benchmark-with-evidence':[],
 'jiuzhaigou-water-and-mountains':[],
 'chaoshan-streets-and-sea':[]
};
const dailyFiles=new Set([
 'content/daily-math/topics.json','content/daily-math/publications.json',
 'scripts/daily-math.mjs','scripts/daily-math-cli.mjs','scripts/daily-math-views.mjs',
 'src/scripts/daily-math.js','src/scripts/daily-math-models.js',
 'src/scripts/daily-math-drawings.js','src/scripts/math-archive.js','src/styles/daily-math.css'
]);

// Rounded relative costs from Actions run 37497187554, not timing assertions.
// Long page sweeps go in different jobs; new suites still participate with a default cost.
const costs={accessibility:75,vision:55,layout:40,'daily-math':35,compile:20,lessons:15,
 autograd:15,dispatch:15,'cuda-streams':15,'reading-demos':15,'process-demos':15,
 browser:15,math:15,pytorch:15,compiler:10,language:10,tensor:10,'code-copy':10,
 discussions:10,sharing:10,support:10,cuda:5,reading:5,analytics:5};
export function partitionSuites(names){
 const groups=Array.from({length:Math.min(3,names.length)},(_,i)=>({group:`browser-${i+1}`,names:[],cost:0}));
 for(const name of [...names].sort((a,b)=>(costs[b]??20)-(costs[a]??20)||a.localeCompare(b))){
  const group=groups.reduce((best,next)=>next.cost<best.cost?next:best);
  group.names.push(name);group.cost+=costs[name]??20;
 }
 return {include:groups.map(({group,names})=>({group,suites:names.join(',')}))};
}

export function planChecks({files=null,mode='auto',suites}){
 if(!['auto','daily','full'].includes(mode))throw Error(`Unknown CI mode: ${mode}`);
 const all=suites.map(suite=>suite.name);
 if(!all.length)throw Error('No browser suites configured');
 const finish=(scope,names,reason)=>{
  const selected=[...new Set(names)];
  for(const name of selected)if(!all.includes(name))throw Error(`CI references an unconfigured suite: ${name}`);
  return {site:scope!=='docs',scope,reason,suites:selected,matrix:partitionSuites(selected)};
 };
 const full=reason=>finish('full',all,reason);
 if(mode==='full')return full('Full regression requested');
 if(mode==='daily')return finish('daily',daily,'Daily mathematics publication');
 if(!files?.length)return full('No reliable changed-file range; use full regression');
 const selected=new Set();
 for(const file of files){
  if(file==='README.md'||/^docs\/.+\.md$/.test(file))continue;
  let related;
  const post=/^content\/posts\/([^/]+)\.html$/.exec(file);
  if(post&&Object.hasOwn(posts,post[1]))related=[...article,...posts[post[1]]];
  else if(/^content\/assets\/(chaoshan|jiuzhaigou)\//.test(file))related=article;
  else if(file.startsWith('content/assets/vision/'))related=[...article,'vision'];
  else if(dailyFiles.has(file))related=daily;
  else if(['scripts/support.mjs','src/scripts/support.js','src/styles/support.css'].includes(file))related=[...core,'support','sharing'];
  else if(file==='src/scripts/analytics.js')related=[...core,'analytics'];
  else{
   const suite=suites.find(suite=>suite.file===file);
   if(suite)related=[...core,suite.name];
  }
  // Templates, common styles/scripts, dependencies, CI, new articles and unknown assets are full checks.
  if(!related)return full(`Shared or unmapped change: ${file}`);
  for(const name of related)selected.add(name);
 }
 if(!selected.size)return finish('docs',[],'Documentation-only change');
 return finish('targeted',[...selected],'Suites selected from changed files');
}

export function changedFiles({eventName,event,head='HEAD',cwd=process.cwd()}){
 const git=args=>execFileSync('git',args,{cwd,encoding:'utf8',stdio:['ignore','pipe','pipe']}).trimEnd();
 try{
  let base;
  if(eventName==='pull_request'){
   const sha=event.pull_request?.base?.sha;
   if(!/^[a-f\d]{40}$/i.test(sha||''))return null;
   base=git(['merge-base',sha,head]);
  }else if(eventName==='push'){
   base=event.before;
   if(!/^[a-f\d]{40}$/i.test(base||'')||/^0+$/.test(base))return null;
  }else return null;
  // Disable rename detection: both the old and new paths must influence test selection.
  return git(['diff','--no-renames','--name-only','-z',base,head,'--']).split('\0').filter(Boolean);
 }catch{return null;}
}

if(process.argv[1]===fileURLToPath(import.meta.url)){
 const {scripts}=JSON.parse(await fs.readFile('package.json','utf8'));
 const event=process.env.GITHUB_EVENT_PATH?JSON.parse(await fs.readFile(process.env.GITHUB_EVENT_PATH,'utf8')):{};
 const plan=planChecks({mode:process.env.CI_CHECK_MODE||'auto',suites:discoverSuites(scripts),files:changedFiles({eventName:process.env.GITHUB_EVENT_NAME,event})});
 console.log(JSON.stringify(plan,null,2));
 if(process.env.GITHUB_OUTPUT)await fs.appendFile(process.env.GITHUB_OUTPUT,`site=${plan.site}\nmatrix=${JSON.stringify(plan.matrix)}\n`);
 if(process.env.GITHUB_STEP_SUMMARY)await fs.appendFile(process.env.GITHUB_STEP_SUMMARY,`## CI scope: ${plan.scope}\n\n${plan.reason}\n\n${plan.suites.length?plan.matrix.include.map(group=>`- **${group.group}**: ${group.suites}`).join('\n'):'Documentation checks only; no build or deployment.'}\n`);
}
