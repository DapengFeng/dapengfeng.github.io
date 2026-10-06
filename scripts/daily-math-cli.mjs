import fs from 'node:fs/promises';
import {loadDailyMath,inventory,planPublications,publicationFile,previousPublications,validateHistory,similarityWarnings} from './daily-math.mjs';
import {dateInfo} from '../src/scripts/daily-math.js';
const [command,...args]=process.argv.slice(2),value=flag=>args[args.indexOf(flag)+1];
const {topics,publications,entries}=await loadDailyMath(),today=dateInfo().date;
if(command==='plan'){
 const next=planPublications(topics,publications,{today,count:args.includes('--count')?Number(value('--count')):Infinity});
 await fs.writeFile(publicationFile,JSON.stringify(next,null,2)+'\n');
 console.log(`Scheduled ${next.length-publications.length} topics; last date: ${next.at(-1)?.date||'none'}.`);
}else if(command==='check'){
 const ref=args.includes('--base')?value('--base'):(process.env.MATH_BASE_REF||'HEAD');
 const previous=previousPublications(ref);if(previous)validateHistory(publications,previous,today);
 for(const warning of similarityWarnings(topics))console.warn(`${process.env.GITHUB_ACTIONS?'::warning::':''}${warning}`);
 const stock=inventory(entries,today),message=`Daily mathematics: ${entries.length} unique scheduled topics; ${stock.remaining} days after today; prepared through ${stock.through||'none'}.`;
 console.log(message);
 if(process.env.GITHUB_STEP_SUMMARY)await fs.appendFile(process.env.GITHUB_STEP_SUMMARY,message+'\n');
 if(stock.remaining<14)console.warn(`${process.env.GITHUB_ACTIONS?'::warning::':''}Daily mathematics inventory is below 14 days. Add reviewed topics and run npm run math:plan.`);
 if(!stock.days&&args.includes('--require-today'))throw Error('No topic scheduled for today. Old topics will not be recycled.');
}else throw Error('Use check or plan');
