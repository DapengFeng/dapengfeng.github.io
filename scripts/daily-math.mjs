import fs from 'node:fs/promises';
import {execFileSync} from 'node:child_process';
import {renderers} from '../src/scripts/daily-math-models.js';
import {parseDate,addDays,dateInfo} from '../src/scripts/daily-math.js';
export const topicFile='content/daily-math/topics.json',publicationFile='content/daily-math/publications.json';
const normalized=s=>s.toLowerCase().replace(/[\s\p{P}\p{S}]/gu,'');
export function validateCatalogue(topics,publications){
 const seen=new Map();
 function unique(kind,value){const key=kind+':'+value;if(seen.has(key))throw Error(`Duplicate ${kind}: ${value}`);seen.set(key,true);}
 for(const t of topics){
  if(!/^[a-z][a-z0-9-]*$/.test(t.id)||!/^[a-z][a-z0-9-]*$/.test(t.concept))throw Error('Invalid topic ID/concept');
  unique('id',t.id);unique('concept',t.concept);
  if(!['ready','draft'].includes(t.status))throw Error(`Invalid status: ${t.id}`);
  if(t.status==='draft')continue;
  if(!renderers.includes(t.renderer))throw Error(`Unknown renderer: ${t.renderer}`);
  parseDate(t.reviewedOn);
  for(const key of ['en','zh','descriptionEn','descriptionZh','readingEn','readingZh','formula'])if(typeof t[key]!=='string'||!t[key].trim())throw Error(`Missing ${key}: ${t.id}`);
  for(const key of ['en','zh','descriptionEn','descriptionZh'])unique(key,normalized(t[key]));
  if(!t.source?.en||!t.source?.zh||new URL(t.source.url).protocol!=='https:')throw Error(`Invalid source: ${t.id}`);
 }
 const byId=new Map(topics.map(t=>[t.id,t]));
 for(const [i,p] of publications.entries()){
  parseDate(p.date);unique('publication date',p.date);unique('published ID',p.id);unique('published concept',p.concept);
  const t=byId.get(p.id);
  if(!t||t.status!=='ready'||t.concept!==p.concept)throw Error(`Publication lacks a reviewed matching topic: ${p.id}`);
  if(t.reviewedOn>p.date)throw Error(`Topic reviewed after its publication: ${p.id}`);
  if(i&&p.date<=publications[i-1].date)throw Error(`Unordered publication date: ${p.date}`);
 }
 return publications.map(p=>({...byId.get(p.id),date:p.date}));
}
export function validateHistory(current,previous,today){
 parseDate(today);
 const old=new Map(previous.map(p=>[p.date,p])),updated=new Map(current.map(p=>[p.date,p]));
 for(const p of previous.filter(p=>p.date<=today)){
  const next=updated.get(p.date);
  if(!next||next.id!==p.id||next.concept!==p.concept)throw Error(`Published history is immutable: ${p.date}`);
 }
 for(const p of current)if(p.date<today&&!old.has(p.date))throw Error(`Cannot backdate a new publication: ${p.date}`);
}
export function planPublications(topics,publications,{today=dateInfo().date,count=Infinity}={}){
 validateCatalogue(topics,publications);
 if(!(count===Infinity||Number.isInteger(count)&&count>0))throw Error('Count must be a positive integer');
 const used=new Set(publications.map(p=>p.id));
 let date=publications.length?addDays(publications.at(-1).date,1):today;
 // Restart today after a lapse, preserving the honest gap in publication history.
 if(date<today)date=today;
 const next=[...publications];
 for(const t of topics.filter(t=>t.status==='ready'&&!used.has(t.id)).slice(0,count)){
  next.push({date,id:t.id,concept:t.concept});date=addDays(date,1);
 }
 validateCatalogue(topics,next);return next;
}
export async function loadDailyMath(){
 const [topics,publications]=await Promise.all([topicFile,publicationFile].map(file=>fs.readFile(file,'utf8').then(JSON.parse)));
 return {topics,publications,entries:validateCatalogue(topics,publications)};
}
export function inventory(entries,today=dateInfo().date){
 let through=today,days=0;const dates=new Set(entries.map(e=>e.date));
 while(dates.has(through)){days++;through=addDays(through,1);}
 return {today,days,through:days?addDays(through,-1):null,remaining:Math.max(0,days-1)};
}
export function previousPublications(ref){
 if(!ref||/^0+$/.test(ref))return null;
 if(!/^[\w./-]+$/.test(ref)||ref.startsWith('-'))throw Error('Invalid history ref');
 execFileSync('git',['rev-parse','--verify',ref+'^{commit}'],{stdio:['ignore','pipe','pipe']});
 const files=execFileSync('git',['ls-tree','--name-only',ref,'--',publicationFile],{encoding:'utf8'}).trim();
 return files?JSON.parse(execFileSync('git',['show',`${ref}:${publicationFile}`],{encoding:'utf8'})):null;
}

// Similar wording is a review signal, not proof that mathematical ideas are equivalent.
export function similarityWarnings(topics){
 const warnings=[],ready=topics.filter(t=>t.status==='ready');
 const wordSets=ready.map(t=>new Set(t.descriptionEn.toLowerCase().match(/[a-z]{3,}/g)||[]));
 for(let i=0;i<ready.length;i++)for(let j=0;j<i;j++){
  const a=wordSets[i],b=wordSets[j],intersection=[...a].filter(w=>b.has(w)).length;
  if(intersection/(a.size+b.size-intersection)>.6)warnings.push(`Similar descriptions require review: ${ready[i].id}, ${ready[j].id}`);
 }
 return warnings;
}
