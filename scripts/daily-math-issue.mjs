import fs from 'node:fs/promises';
import {fileURLToPath} from 'node:url';
import {loadDailyMath,inventory,MIN_FUTURE_DAYS} from './daily-math.mjs';
import {addDays,dateInfo} from '../src/scripts/daily-math.js';

export const ISSUE_MARKER='<!-- daily-math-inventory:v1 -->';
export const ISSUE_TITLE='每日数学待补充 / Replenish daily mathematics';

export function issueBody(stock,repository){
 const low=stock.remaining<MIN_FUTURE_DAYS;
 return `${ISSUE_MARKER}
${low?'每日数学余量不足，请补充经过核对的新主题。 / Reviewed topics need replenishment.':'每日数学余量已恢复。 / The reviewed queue has been replenished.'}

| 检查项 / Check | 当前状态 / Status |
| --- | --- |
| 北京时间日期 / Shanghai date | ${stock.today} |
| 今日排期 / Today scheduled | ${stock.days?'有 / Yes':'缺失 / Missing'} |
| 连续未来余量（不含今天） / Consecutive future days | ${stock.remaining} |
| 连续覆盖至 / Covered through | ${stock.through||'无 / None'} |
| 首个缺档日 / First uncovered date | ${stock.through?addDays(stock.through,1):stock.today} |
| 提醒阈值 / Alert threshold | < ${MIN_FUTURE_DAYS} days |

余量从今天起连续计算；遇到缺档即停止。草稿、未排期主题和缺档后的排期不计入余量。
Stock counts consecutive scheduled days from today and stops at the first gap. Drafts, unscheduled topics and dates after a gap do not count.

补充有实质差异的主题，核对来源、公式、双语说明和图示后，再运行：
Review distinct topics, sources, formulas, both languages and diagrams, then run:

\x60\x60\x60sh
npm run math:plan
npm run math:check -- --require-today
\x60\x60\x60

按[每日数学规范 / Daily mathematics guide](https://github.com/${repository}/blob/HEAD/docs/authoring/daily-mathematics.md)完成验证，再将内容和排期合入默认分支。
Validate and merge the content and schedule into the default branch following the guide.

这条 issue 由自动检查维护；补足后自动关闭，再次不足会重新打开。 / This issue is updated automatically, closed after replenishment and reopened when stock runs low again.
`;
}

function validIssue(issue){
 return issue&&Number.isSafeInteger(issue.number)&&issue.number>0&&
  ['open','closed'].includes(issue.state)&&(issue.body===null||typeof issue.body==='string');
}

export async function syncInventoryIssue({stock,repository,token,fetchImpl=fetch}){
 if(!/^[A-Za-z0-9_.-]+\/[A-Za-z0-9_.-]+$/.test(repository||''))throw Error('GITHUB_REPOSITORY must be owner/repository.');
 if(!token)throw Error('GITHUB_TOKEN is required for issue synchronization.');
 const root=`/repos/${repository}/issues`;
 async function request(method,path,body){
  let response;
  try{
   response=await fetchImpl(`https://api.github.com${path}`,{
    method,redirect:'error',signal:AbortSignal.timeout(20000),
    headers:{Accept:'application/vnd.github+json',Authorization:`Bearer ${token}`,
     'Content-Type':'application/json','X-GitHub-Api-Version':'2022-11-28'},
    ...(body?{body:JSON.stringify(body)}:{})
   });
  }catch{throw Error(`GitHub API ${method} request failed (network or timeout).`);}
  // Never print credentials or upstream response bodies; never retry a possibly successful POST.
  if(!response.ok)throw Error(`GitHub API ${method} failed: HTTP ${response.status}. Check Issues availability and workflow permissions.`);
  try{return await response.json();}catch{throw Error(`GitHub API ${method} returned invalid JSON.`);}
 }

 // REST listing includes PRs. Search all pages and states before deciding to create anything.
 const matches=[];
 for(let page=1;;page++){
  const issues=await request('GET',`${root}?state=all&sort=created&direction=desc&per_page=100&page=${page}`);
  if(!Array.isArray(issues)||!issues.every(validIssue))throw Error('GitHub issue list has an unexpected shape; no issues were changed.');
  matches.push(...issues.filter(i=>!i.pull_request&&i.body?.includes(ISSUE_MARKER)));
  if(issues.length<100)break;
 }
 const existing=matches.find(i=>i.state==='open')||matches[0];
 const low=stock.remaining<MIN_FUTURE_DAYS,body=issueBody(stock,repository);
 const result={stock,body,issueNumber:existing?.number};
 if(!low&&(!existing||existing.state==='closed'))return {...result,action:'healthy'};
 if(low&&!existing){
  const created=await request('POST',root,{title:ISSUE_TITLE,body});
  if(!validIssue(created))throw Error('GitHub create response has an unexpected shape; check the repository before retrying.');
  return {...result,action:'created',issueNumber:created.number};
 }
 let action='updated',patch={body};
 if(!low){action='closed';patch={body,state:'closed',state_reason:'completed'};}
 else if(existing.state==='closed'){action='reopened';patch={body,state:'open',state_reason:'reopened'};}
 else if(existing.body===body)return {...result,action:'unchanged'};
 await request('PATCH',`${root}/${existing.number}`,patch);
 return {...result,action};
}

export async function checkInventory({load=loadDailyMath,today=dateInfo().date,
 repository=process.env.GITHUB_REPOSITORY,token=process.env.GITHUB_TOKEN,dryRun=false,fetchImpl=fetch}={}){
 // Catalogue validation must finish before any API call, especially an automatic closure.
 const {entries}=await load(),stock=inventory(entries,today);
 if(dryRun)return {action:'dry-run',stock,body:issueBody(stock,repository||'OWNER/REPOSITORY')};
 return syncInventoryIssue({stock,repository,token,fetchImpl});
}

if(process.argv[1]===fileURLToPath(import.meta.url)){
 try{
  const args=process.argv.slice(2);
  if(args.length&&!(args.length===1&&args[0]==='--dry-run'))throw Error('Usage: node scripts/daily-math-issue.mjs [--dry-run]');
  const result=await checkInventory({dryRun:args.includes('--dry-run')});
  const message=`Daily mathematics: ${result.stock.remaining} consecutive future days; today ${result.stock.today}; through ${result.stock.through||'none'}; ${result.action}${result.issueNumber?` #${result.issueNumber}`:''}.`;
  console.log(message);
  if(result.action==='dry-run')console.log(result.body);
  if(process.env.GITHUB_STEP_SUMMARY)await fs.appendFile(process.env.GITHUB_STEP_SUMMARY,message+'\n');
 }catch(error){console.error(error.message);process.exitCode=1;}
}
