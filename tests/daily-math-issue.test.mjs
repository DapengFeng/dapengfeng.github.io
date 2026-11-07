import test from 'node:test';
import assert from 'node:assert/strict';
import {checkInventory,ISSUE_MARKER,ISSUE_TITLE,syncInventoryIssue} from '../scripts/daily-math-issue.mjs';
import {inventory,loadDailyMath,MIN_FUTURE_DAYS,validateCatalogue} from '../scripts/daily-math.mjs';
import {addDays} from '../src/scripts/daily-math.js';

const repository='example/knowledge-lab',token='test-token-not-a-real-credential';
const start='2030-12-25';
const entriesFor=(count,date=start)=>Array.from({length:count},(_,i)=>({date:addDays(date,i)}));
const response=(value,status=200)=>new Response(JSON.stringify(value),{status,headers:{'Content-Type':'application/json'}});
const issue=(number,{state='open',body='Unrelated issue',...rest}={})=>({number,state,body,...rest});

function github(initial=[]){
 const issues=structuredClone(initial),calls=[];
 let next=Math.max(0,...issues.map(item=>item.number))+1;
 return {issues,calls,fetchImpl:async(url,options)=>{
  const target=new URL(url),method=options.method,body=options.body?JSON.parse(options.body):undefined;
  calls.push({method,url,body});
  assert.equal(target.origin,'https://api.github.com');
  assert.equal(options.headers.Authorization,`Bearer ${token}`);
  assert.equal(options.redirect,'error');
  if(method==='GET'){
   assert.equal(target.searchParams.get('state'),'all');
   const page=Number(target.searchParams.get('page'));
   return response([...issues].sort((a,b)=>b.number-a.number).slice((page-1)*100,page*100));
  }
  if(method==='POST'){
   const created=issue(next++,body);issues.push(created);return response(created,201);
  }
  if(method==='PATCH'){
   const found=issues.find(item=>item.number===Number(target.pathname.split('/').at(-1)));
   assert.ok(found,'only an existing issue may be patched');Object.assign(found,body);return response(found);
  }
  assert.fail(`Unexpected method ${method}`);
 }};
}
const sync=(api,stock)=>syncInventoryIssue({stock,repository,token,fetchImpl:api.fetchImpl});

test('inventory alerts below fourteen consecutive future days, including gaps and exhaustion',async()=>{
 assert.equal(MIN_FUTURE_DAYS,14);
 for(const [count,expected]of [[14,'created'],[15,'healthy'],[16,'healthy']]){
  const api=github(),result=await checkInventory({load:async()=>({entries:entriesFor(count)}),today:start,repository,token,fetchImpl:api.fetchImpl});
  assert.equal(result.stock.remaining,count-1);assert.equal(result.action,expected);
 }
 const cases=[
  {entries:entriesFor(30).filter(item=>item.date!==addDays(start,3)),remaining:2,through:addDays(start,2)},
  {entries:entriesFor(30).slice(1),remaining:0,through:null},
  {entries:entriesFor(30).filter(item=>item.date!==addDays(start,1)),remaining:0,through:start},
  {entries:entriesFor(10,addDays(start,-10)),remaining:0,through:null},
  {entries:[],remaining:0,through:null}
 ];
 for(const {entries,remaining,through}of cases){
  const result=await checkInventory({load:async()=>({entries}),today:start,dryRun:true});
  assert.equal(result.stock.remaining,remaining);assert.equal(result.stock.through,through);
  assert.ok(result.body.includes(through?addDays(through,1):start),'the reminder identifies the first uncovered date');
 }
 assert.equal(inventory(entriesFor(16,'2028-02-28'),'2028-02-29').remaining,14,'leap day does not shorten coverage');
});

test('default inventory date changes at Shanghai midnight rather than UTC midnight',async t=>{
 const entries=entriesFor(16,'2030-12-31');
 const clock=t.mock.method(Date,'now',()=>Date.parse('2030-12-31T15:59:59.999Z'));
 const before=await checkInventory({load:async()=>({entries}),dryRun:true});
 clock.mock.mockImplementation(()=>Date.parse('2030-12-31T16:00:00.000Z'));
 const after=await checkInventory({load:async()=>({entries}),dryRun:true});
 assert.equal(before.stock.today,'2030-12-31');assert.equal(before.stock.remaining,15);
 assert.equal(after.stock.today,'2031-01-01');assert.equal(after.stock.remaining,14);
});

test('one reminder is created, unchanged on retry, updated daily, closed and reopened',async()=>{
 const api=github(),entries=entriesFor(14),low=inventory(entries,start);
 let result=await sync(api,low);
 assert.equal(result.action,'created');assert.equal(api.issues.length,1);const number=result.issueNumber;
 assert.equal(api.issues[0].title,ISSUE_TITLE);assert.ok(api.issues[0].body.startsWith(ISSUE_MARKER));
 result=await sync(api,low);assert.equal(result.action,'unchanged');
 assert.equal(api.calls.filter(call=>call.method!=='GET').length,1,'retry makes no extra mutation or comment');
 result=await sync(api,inventory(entries,addDays(start,1)));
 assert.equal(result.action,'updated');assert.equal(result.issueNumber,number);
 assert.ok(api.issues[0].body.includes(addDays(start,1)));
 result=await sync(api,inventory(entriesFor(15),start));
 assert.equal(result.action,'closed');assert.equal(api.issues[0].state,'closed');assert.equal(api.issues[0].state_reason,'completed');
 result=await sync(api,inventory(entriesFor(15),start));assert.equal(result.action,'healthy');
 result=await sync(api,low);assert.equal(result.action,'reopened');assert.equal(result.issueNumber,number);
 assert.equal(api.issues[0].state,'open');assert.equal(api.issues.length,1);
 assert.equal(api.calls.filter(call=>call.method==='POST').length,1);
 assert.equal(api.calls.filter(call=>call.method==='PATCH').length,3);
});

test('the marker identifies reminders across pages without changing PRs or same-title issues',async()=>{
 const initial=Array.from({length:100},(_,i)=>issue(i+2));
 initial[0]=issue(2,{title:ISSUE_TITLE,body:null});
 initial[1]=issue(3,{body:ISSUE_MARKER,pull_request:{url:'https://api.github.com/repos/example/knowledge-lab/pulls/3'}});
 initial.push(issue(1,{state:'closed',body:ISSUE_MARKER}));
 const api=github(initial),result=await sync(api,inventory(entriesFor(14),start));
 assert.equal(result.action,'reopened');assert.equal(result.issueNumber,1);
 assert.equal(api.calls.filter(call=>call.method==='GET').length,2);
 assert.equal(api.calls.filter(call=>call.method==='POST').length,0);
 assert.equal(api.issues.find(item=>item.number===2).body,null);
 assert.equal(api.issues.find(item=>item.number===3).body,ISSUE_MARKER);
 assert.equal(api.calls.at(-1).url,`https://api.github.com/repos/${repository}/issues/1`);
 const unrelated=github([issue(1,{title:ISSUE_TITLE}),issue(2,{body:ISSUE_MARKER,pull_request:{}})]);
 assert.equal((await sync(unrelated,inventory(entriesFor(14),start))).action,'created');
});

test('an open reminder takes precedence over a more recently created closed reminder',async()=>{
 const api=github([issue(1,{body:ISSUE_MARKER}),issue(2,{body:ISSUE_MARKER,state:'closed'})]);
 const result=await sync(api,inventory(entriesFor(14),start));
 assert.equal(result.issueNumber,1);assert.equal(result.action,'updated');
 assert.equal(api.issues.find(item=>item.number===2).state,'closed');
});

test('failed or malformed listings never create or close an issue, including later pages',async()=>{
 const failures=[
  ...[401,403,429,500].map(status=>()=>response({message:token},status)),
  ()=>new Response(token,{status:200}),
  ()=>response({issues:[]}),
  ()=>response([issue(0)]),
  ()=>response([{number:1,state:'unexpected',body:ISSUE_MARKER}]),
  ()=>{throw Error(`Network failure: ${token}`);}
 ];
 for(const fail of failures)for(const count of [14,15]){
  const methods=[];
  await assert.rejects(syncInventoryIssue({stock:inventory(entriesFor(count),start),repository,token,fetchImpl:async(_url,options)=>{methods.push(options.method);return fail();}}),error=>{
   assert.ok(!error.message.includes(token),'credentials and upstream response bodies stay out of errors');return true;
  });
  assert.deepEqual(methods,['GET']);
 }
 const calls=[];
 await assert.rejects(syncInventoryIssue({stock:inventory(entriesFor(14),start),repository,token,fetchImpl:async(_url,options)=>{
  calls.push(options.method);return calls.length===1?response(Array.from({length:100},(_,i)=>issue(i+1))):response({},503);
 }}),/HTTP 503/);
 assert.deepEqual(calls,['GET','GET']);
});

test('failed mutations are reported without retrying a possibly successful creation',async()=>{
 for(const existing of [[],[issue(1,{body:ISSUE_MARKER})]]){
  const methods=[];
  await assert.rejects(syncInventoryIssue({stock:inventory(entriesFor(14),start),repository,token,fetchImpl:async(_url,options)=>{
   methods.push(options.method);if(options.method==='GET')return response(existing);throw Error('Connection lost after sending the mutation');
  }}),/request failed/);
  assert.deepEqual(methods,['GET',existing.length?'PATCH':'POST']);
 }
 const methods=[];
 await assert.rejects(syncInventoryIssue({stock:inventory(entriesFor(14),start),repository,token,fetchImpl:async(_url,options)=>{
  methods.push(options.method);return response(options.method==='GET'?[]:{number:7});
 }}),/create response/);
 assert.deepEqual(methods,['GET','POST']);
});

test('catalogue validation must succeed before either an API mutation or dry-run output',async()=>{
 const {topics,publications}=await loadDailyMath(),invalid=structuredClone(publications);
 invalid[0].concept='not-the-reviewed-concept';
 let calls=0;
 for(const dryRun of [false,true])await assert.rejects(checkInventory({
  load:async()=>({entries:validateCatalogue(topics,invalid)}),today:start,repository,token,dryRun,
  fetchImpl:async()=>{calls++;assert.fail('invalid content must never contact GitHub');}
 }),/reviewed matching topic/);
 assert.equal(calls,0);
});

test('dry-run needs no credentials or network; live synchronization validates configuration first',async()=>{
 const fetchImpl=async()=>assert.fail('no network is expected');
 const result=await checkInventory({load:async()=>({entries:entriesFor(14)}),today:start,repository:'',token:'',dryRun:true,fetchImpl});
 assert.equal(result.action,'dry-run');assert.equal(result.stock.remaining,13);
 assert.ok(result.body.includes('github.com/OWNER/REPOSITORY/blob/HEAD/docs/authoring/daily-mathematics.md'));
 for(const invalid of ['', '../owner/repo','owner/repo/extra','owner/repo?token=secret']){
  await assert.rejects(syncInventoryIssue({stock:result.stock,repository:invalid,token,fetchImpl}),/owner\/repository/);
 }
 await assert.rejects(syncInventoryIssue({stock:result.stock,repository,token:'',fetchImpl}),/GITHUB_TOKEN/);
});
