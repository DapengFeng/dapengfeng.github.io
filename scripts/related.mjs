// Build-time bilingual lexical similarity. No browser requests or visitor tracking.
const segmenter=new Intl.Segmenter('zh',{granularity:'word'});
const stop=new Set(`a an and are as at be been being by can could did do does for from had has have how if in into is it its may more most not of on one or our out should so than that the their them then there these they this those to two use used using was we were what when where which while will with would you your about also each only same such through very within without article example figure source code output input result results chapter section note notes read reading show shown following first second next data
的 了 和 是 在 与 为 将 对 把 也 中 到 有 从 可以 这个 一个 我们 你 它 其 而 或 并 但 不 会 上 下 时 就 都 等 这些 那些 其中 如果 因此 需要 使用 通过 以及 这里 下面 上面 本文 示例 代码 结果 输出 输入 数据 内容 文章 章节 图示 说明 阅读 显示`.split(/\s+/));
function tokens(text=''){
 const normalized=String(text).normalize('NFKC').toLowerCase().replace(/c\+\+/g,' cpp ').replace(/https?:\/\/\S+/g,' ');
 return [...segmenter.segment(normalized)].filter(s=>s.isWordLike).map(s=>s.segment).filter(s=>!stop.has(s)&&!/^\d+$/.test(s)&&([...s].length>1||/\p{Script=Han}/u.test(s)));
}
function counts(text){const map=new Map();for(const t of tokens(text))map.set(t,(map.get(t)||0)+1);return map;}
function unit(vector){const norm=Math.hypot(...vector.values());return new Map([...vector].map(([t,v])=>[t,norm?v/norm:0]));}
function similarity(a,b){let sum=0;for(const [t,v]of a)sum+=v*(b.get(t)||0);return sum;}
function tagSet(post){return new Set((post.tags||[]).map(t=>t.normalize('NFKC').trim().toLowerCase()));}
export function createRelatedRecommender(posts){
 const eligible=posts.filter(p=>!p.draft&&!p.archiveOnly),documents=new Map(),frequency=new Map();
 for(const post of eligible){
  const fields=[
   [4,counts(`${post.titleEn||''} ${post.title||''}`)],
   [3,counts((post.tags||[]).join(' '))],
   [2,counts(`${post.descriptionEn||''} ${post.description||''}`)],
   [1,counts((post.headings||[]).map(h=>`${h.titleEn||''} ${h.titleZh||h.title||''}`).join(' '))],
   [2,counts(post.relatedText||'')]
  ];
  documents.set(post.slug,fields);
  for(const t of new Set(fields.flatMap(([,map])=>[...map.keys()])))frequency.set(t,(frequency.get(t)||0)+1);
 }
 const vectors=new Map(),tags=new Map();
 for(const post of eligible){
  const combined=new Map();
  for(const [weight,field]of documents.get(post.slug)){
   // Normalize each field so a long article cannot win just by repeating words.
   const weighted=unit(new Map([...field].map(([t,n])=>[t,(1+Math.log(n))*Math.log(1+eligible.length/(frequency.get(t)||1))])));
   for(const [t,v]of weighted)combined.set(t,(combined.get(t)||0)+weight*v);
  }
  vectors.set(post.slug,unit(combined));tags.set(post.slug,tagSet(post));
 }
 return (post,limit=3)=>{
  const vector=vectors.get(post.slug);if(!vector)return [];
  const ownTags=tags.get(post.slug);
  return eligible.filter(p=>p.slug!==post.slug).map(candidate=>{
   const shared=[...tags.get(candidate.slug)].filter(t=>ownTags.has(t)).length;
   const union=new Set([...ownTags,...tags.get(candidate.slug)]).size;
   const content=similarity(vector,vectors.get(candidate.slug));
   const sameSeries=Boolean(post.series?.id&&post.series.id===candidate.series?.id);
   const series=sameSeries ? .24+.08/(1+Math.abs(post.series.part-candidate.series.part)):0;
   const score=.75*content+.2*(union?shared/union:0)+series+(post.category===candidate.category ? .03 : 0);
   return {candidate,score,relevant:sameSeries||shared>0||content>=.065};
  }).filter(item=>item.relevant).sort((a,b)=>b.score-a.score||a.candidate.slug.localeCompare(b.candidate.slug))
   .slice(0,Math.max(0,limit)).map(item=>item.candidate);
 };
}
