import {load} from 'cheerio';
import {escape as e} from './config.mjs';

export function collectSeries(posts) {
 const groups=new Map();
 for(const post of posts){
  const s=post.series;if(!s)continue;
  if(!/^[a-z0-9]+(?:-[a-z0-9]+)*$/.test(s.id||'')||!s.title?.trim()||!s.titleEn?.trim()||!Number.isInteger(s.part)||s.part<1)throw Error(`${post.slug}: series requires id, bilingual titles, and a positive integer part`);
  if(!groups.has(s.id))groups.set(s.id,{id:s.id,title:s.title,titleEn:s.titleEn,url:`/series/${s.id}/`,posts:[],roadmap:new Map()});
  const group=groups.get(s.id);
  if(group.title!==s.title||group.titleEn!==s.titleEn)throw Error(`Inconsistent series title: ${s.id}`);
  if(group.posts.some(p=>p.series.part===s.part))throw Error(`Duplicate series part: ${s.id}/${s.part}`);
  if(s.roadmap!==undefined&&!Array.isArray(s.roadmap))throw Error(`${post.slug}: series roadmap must be an array`);
  for(const entry of s.roadmap||[]){
   if(!Number.isInteger(entry.part)||entry.part<1||!entry.title?.trim()||!entry.titleEn?.trim())throw Error(`${post.slug}: roadmap requires a positive part and bilingual titles`);
   const previous=group.roadmap.get(entry.part);
   if(previous&&(previous.title!==entry.title||previous.titleEn!==entry.titleEn))throw Error(`Conflicting roadmap part: ${s.id}/${entry.part}`);
   group.roadmap.set(entry.part,entry);
  }
  group.posts.push(post);
 }
 return [...groups.values()].map(group=>({...group,posts:group.posts.sort((a,b)=>a.series.part-b.series.part)}));
}
export function seriesMetadata(group){
 return [group.titleEn,group.title,`Read ${group.titleEn} in order: source-guided articles, working examples, and interactive explanations.`,`按期阅读《${group.title}》：结合源码、可运行示例与交互图解循序学习。`];
}

// Resolve roadmap entries from published metadata, without article-specific URLs.
export function renderSeriesPreviews(post, posts) {
 if(!post.html.includes('data-series-preview'))return post.html;
 const $=load(post.html,null,false),groups=collectSeries(posts);
 const pair=(en,zh)=>`<span data-lang="en" lang="en">${e(en)}</span><span data-lang="zh" lang="zh-CN">${e(zh)}</span>`;
 $('[data-series-preview]').each((_,element)=>{
  const preview=$(element),id=preview.attr('data-series-preview')||post.series?.id;
  const group=groups.find(g=>g.id===id);
  if(!group)throw Error(`${post.slug}: unknown preview series ${id}`);
  const planned=new Map(group.roadmap),explicit=new Set();
  preview.children('[data-series-part]').each((_,entry)=>{
   const node=$(entry),part=Number(node.attr('data-series-part'));
   if(!Number.isInteger(part)||part<1||explicit.has(part))throw Error(`${post.slug}: invalid or duplicate preview part ${part}`);
   explicit.add(part);
   planned.set(part,{titleEn:node.find('[data-lang=en]').first().text(),title:node.find('[data-lang=zh]').first().text()});
  });
  const published=new Map(group.posts.map(p=>[p.series.part,p]));
  const parts=[...new Set([...planned.keys(),...published.keys()])].sort((a,b)=>a-b);
  preview.addClass('series-preview').html(parts.map(part=>{
   const item=published.get(part),titles=item||planned.get(part);
   const contents=`<b>${String(part).padStart(2,'0')}</b><span class="series-preview-title">${pair(titles.titleEn,titles.title)}</span><span class="series-preview-status">${item?`<time datetime="${e(item.date)}">${e(item.date)}</time>`:pair('Planned','待发布')}</span>`;
   return item?`<a data-series-part="${part}" data-state="published" href="${e(item.url)}"${item.slug===post.slug?' aria-current="page"':''}>${contents}</a>`:`<div data-series-part="${part}" data-state="planned">${contents}</div>`;
  }).join(''));
 });
 return $.html();
}
