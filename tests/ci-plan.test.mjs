import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {execFileSync,spawnSync} from 'node:child_process';
import {changedFiles,planChecks} from '../scripts/ci-plan.mjs';
import {discoverSuites} from '../scripts/check-browsers.mjs';
import {checkDocs} from '../scripts/check-docs.mjs';

const {scripts}=JSON.parse(await fs.readFile(new URL('../package.json',import.meta.url),'utf8'));
const suites=discoverSuites(scripts);
const plan=(files,mode='auto')=>planChecks({files,mode,suites});

test('documentation-only changes skip site jobs; mixed article edits retain relevant checks',()=>{
 const docs=plan(['README.md','docs/operations/analytics.md']);
 assert.equal(docs.site,false);assert.deepEqual(docs.matrix,{include:[]});
 const vision=plan(['docs/authoring/presentation.md','content/posts/human-visual-system.html']);
 assert.equal(vision.site,true);
 for(const name of ['browser','reading','layout','accessibility','sharing','vision'])assert.ok(vision.suites.includes(name),name);
 assert.ok(!vision.suites.includes('autograd'));
 const mixed=plan(['content/posts/pytorch-04-autograd-engine.html','content/posts/human-visual-system.html']);
 assert.ok(mixed.suites.includes('autograd')&&mixed.suites.includes('vision'));
});

test('analytics and daily publications keep basic coverage without unrelated article suites',()=>{
 assert.deepEqual(new Set(plan(['src/scripts/analytics.js']).suites),new Set(['browser','language','reading','analytics']));
 for(const selection of [plan(null,'daily'),plan(['content/daily-math/topics.json'])]){
  assert.deepEqual(new Set(selection.suites),new Set(['daily-math','browser','language','reading']));
 }
 assert.ok(plan(['content/assets/vision/eye-cutaway.webp']).suites.includes('vision'));
 assert.ok(plan(['content/assets/chaoshan/new-photo.webp']).suites.includes('sharing'));
 assert.ok(plan(['tests/compile.mjs']).suites.includes('compile'));
 for(const file of ['scripts/support.mjs','src/scripts/support.js','src/styles/support.css']){
  assert.deepEqual(new Set(plan([file]).suites),new Set(['browser','language','reading','support','sharing']));
 }
});

test('uncertain changes and shared infrastructure fall back to every configured suite',()=>{
 for(const files of [null,[],['package-lock.json'],['src/styles/site.css'],['src/templates/layout.html'],
  ['scripts/build.mjs'],['.github/workflows/deploy.yml'],['README.md','new-config.json'],
  ['content/posts/new-article.html'],['content/assets/new-series/photo.webp'],['docs/script.mjs']]){
  assert.equal(plan(files).scope,'full',JSON.stringify(files));
  assert.deepEqual(new Set(plan(files).suites),new Set(suites.map(s=>s.name)));
 }
 assert.equal(plan(['README.md'],'full').site,true);
 assert.throws(()=>plan([], 'typo'),/Unknown CI mode/);
});

test('three browser groups cover the selection exactly once, including newly registered suites',()=>{
 const extra=[...suites,{name:'new-suite',file:'tests/new-suite.mjs'}];
 const full=planChecks({mode:'full',suites:extra});
 assert.equal(full.matrix.include.length,3);
 assert.deepEqual(full.matrix.include.flatMap(group=>group.suites.split(',')).sort(),extra.map(s=>s.name).sort());
 for(const selection of [plan(['src/scripts/analytics.js']),plan(null,'daily')]){
  assert.equal(selection.matrix.include.length,3);
  assert.deepEqual(selection.matrix.include.flatMap(group=>group.suites.split(',')).sort(),[...selection.suites].sort());
 }
});

test('git comparison includes the entire push, both sides of renames, and PR changes since divergence',async()=>{
 const cwd=await fs.mkdtemp(path.join(os.tmpdir(),'feng-ci-diff-'));
 const git=(...args)=>execFileSync('git',args,{cwd,encoding:'utf8',stdio:['ignore','pipe','pipe']}).trim();
 const commit=message=>{git('add','.');git('commit','-qm',message);return git('rev-parse','HEAD');};
 try{
  git('init','-q','-b','main');git('config','user.email','ci@example.invalid');git('config','user.name','CI test');
  await fs.mkdir(path.join(cwd,'docs'));await fs.writeFile(path.join(cwd,'old-script.mjs'),'source');
  const before=commit('base');
  await fs.writeFile(path.join(cwd,'README.md'),'readme');commit('first change');
  git('mv','old-script.mjs','docs/moved.md');const head=commit('rename');
  assert.deepEqual(changedFiles({eventName:'push',event:{before},cwd}).sort(),['README.md','docs/moved.md','old-script.mjs']);
  assert.equal(plan(changedFiles({eventName:'push',event:{before},cwd})).scope,'full','moving code to docs is not docs-only');
  git('checkout','-q','-b','base-advanced',before);
  await fs.writeFile(path.join(cwd,'base-only.txt'),'unrelated');const base=commit('base advanced');
  const pr=changedFiles({eventName:'pull_request',event:{pull_request:{base:{sha:base}}},head,cwd});
  assert.deepEqual(pr.sort(),['README.md','docs/moved.md','old-script.mjs']);
  for(const bad of ['0'.repeat(40),'f'.repeat(40),'bad-ref'])assert.equal(changedFiles({eventName:'push',event:{before:bad},cwd}),null);
 }finally{await fs.rm(cwd,{recursive:true,force:true});}
});

test('documentation validation detects broken local links and unclosed fences without following external URLs',async()=>{
 const root=await fs.mkdtemp(path.join(os.tmpdir(),'feng-docs-'));
 try{
  await fs.mkdir(path.join(root,'docs'));
  await fs.writeFile(path.join(root,'README.md'),'# Docs\n\n[Guide](docs/guide.md)\n');
  await fs.writeFile(path.join(root,'docs/guide.md'),'# Guide\n\n[Home](../README.md) [Web](https://example.invalid/)\n\n```md\n[Example](missing.md)\n```\n');
  assert.equal(await checkDocs(root),2);
  await fs.appendFile(path.join(root,'docs/guide.md'),'[Missing](not-found.md)\n');
  await assert.rejects(checkDocs(root),/missing local link target not-found.md/);
  await fs.writeFile(path.join(root,'docs/guide.md'),'```js\nconst a = 1;\n');
  await assert.rejects(checkDocs(root),/unclosed code fence/);
 }finally{await fs.rm(root,{recursive:true,force:true});}
});

test('the actual workflow gate rejects failed, cancelled, or unexpectedly skipped site checks',async()=>{
 const workflow=await fs.readFile(new URL('../.github/workflows/checks.yml',import.meta.url),'utf8');
 const gate=/node <<'JS'\n([\s\S]+?)\n\s+JS/.exec(workflow)?.[1];
 assert.ok(gate,'workflow contains the gate script');
 for(const SITE of ['true','false',''])for(const PLAN_RESULT of ['success','failure','cancelled','skipped']){
  for(const [BUILD_RESULT,BROWSER_RESULT]of [['success','success'],['failure','skipped'],['success','failure'],['success','cancelled'],['success','skipped'],['skipped','skipped']]){
   const run=spawnSync(process.execPath,['-e',gate],{encoding:'utf8',env:{...process.env,SITE,PLAN_RESULT,BUILD_RESULT,BROWSER_RESULT}});
   const expected=PLAN_RESULT==='success'&&((SITE==='true'&&BUILD_RESULT==='success'&&BROWSER_RESULT==='success')||(SITE==='false'&&BUILD_RESULT==='skipped'&&BROWSER_RESULT==='skipped'));
   assert.equal(run.status===0,expected,JSON.stringify({SITE,PLAN_RESULT,BUILD_RESULT,BROWSER_RESULT}));
  }
 }
});
