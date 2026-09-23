export function collectSeries(posts) {
 const groups=new Map();
 for(const post of posts){
  const s=post.series;if(!s)continue;
  if(!/^[a-z0-9]+(?:-[a-z0-9]+)*$/.test(s.id||'')||!s.title?.trim()||!s.titleEn?.trim()||!Number.isInteger(s.part)||s.part<1)throw Error(`${post.slug}: series requires id, bilingual titles, and a positive integer part`);
  if(!groups.has(s.id))groups.set(s.id,{id:s.id,title:s.title,titleEn:s.titleEn,url:`/series/${s.id}/`,posts:[]});
  const group=groups.get(s.id);
  if(group.title!==s.title||group.titleEn!==s.titleEn)throw Error(`Inconsistent series title: ${s.id}`);
  if(group.posts.some(p=>p.series.part===s.part))throw Error(`Duplicate series part: ${s.id}/${s.part}`);
  group.posts.push(post);
 }
 return [...groups.values()].map(group=>({...group,posts:group.posts.sort((a,b)=>a.series.part-b.series.part)}));
}
export function seriesMetadata(group){
 return [group.titleEn,group.title,`Read ${group.titleEn} in order: source-guided articles, working examples, and interactive explanations.`,`按期阅读《${group.title}》：结合源码、可运行示例与交互图解循序学习。`];
}
