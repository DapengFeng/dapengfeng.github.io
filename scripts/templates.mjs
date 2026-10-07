import {readFileSync} from 'node:fs';
const layoutTemplate = readFileSync(new URL('../src/templates/layout.html', import.meta.url), 'utf8');
import {pair} from './i18n.mjs';
import {editorial,curatedRelated,seriesGuide} from './editorial-views.mjs';
import {mathNote,mathInitial} from './daily-math-views.mjs';
import {discussions} from './discussions.mjs';
import {shareTrigger,articleSharing} from './sharing.mjs';
import {loadSupport,supportTrigger,supportDialog} from './support.mjs';
import {collectSeries,seriesMetadata,renderSeriesPreviews} from './series.mjs';
import { categories, site, escape as e } from './config.mjs';
export function art(type, id = '') {
 const key = `g-${type}-${id}`;
 let shapes = '';
 if (type === 'spike') {
   for(let i=0;i<5;i++) {const y=45+i*27; let d=`M0 ${y}`; for(let x=0;x<600;x+=3){let phase=(x+25*i)%96; let v=phase>73&&phase<85?Math.sin((phase-73)/12*Math.PI)*25:Math.sin(x*.022+i)*3; d+=` L${x} ${y-v}`;} shapes+=`<path d="${d}" fill="none" stroke="${i===2?'#c5fb7c':'#658960'}" stroke-width="${i===2?2:1}" opacity="${i===2?1:.5}"/>`;}
   shapes+='<circle cx="313" cy="77" r="4" fill="#d8ff9e"/><circle cx="313" cy="77" r="11" fill="#bcf46e22"/>';
 } else if(type === 'vision') {
   shapes='<path d="M155 100Q300 -15 445 100Q300 215 155 100Z" fill="#202b27" stroke="#efb2bd" stroke-width="2"/><circle cx="300" cy="100" r="54" fill="#345142" stroke="#b9eaa4" stroke-width="2"/><circle cx="300" cy="100" r="24" fill="#0f1912"/><path d="M50 75H154M50 100H180M50 125H154M445 100L485 65H550M445 100L500 100H565M445 100L485 135H550" fill="none" stroke="#c0f47b" stroke-width="2"/><circle cx="320" cy="80" r="8" fill="#deefe2"/>';
 } else if(type === 'rust') {
   shapes='<circle cx="205" cy="100" r="51" fill="none" stroke="#bb96ef" stroke-width="1"/><circle cx="205" cy="100" r="59" fill="none" stroke="#bb96ef" stroke-width="1" stroke-dasharray="3 7"/><text x="205" y="113" text-anchor="middle" fill="#d1b4f5" font-size="41" font-weight="700" font-family="monospace">R</text><path d="M366 47 412 73 412 126 366 153 320 126 320 73Z" fill="none" stroke="#8ca6f1"/><text x="366" y="111" fill="#a7bffe" text-anchor="middle" font-size="29" font-family="monospace">C++</text><path d="M270 100H302M296 95l6 5-6 5" fill="none" stroke="#708072"/><text x="286" y="125" fill="#78817a" font-size="10" text-anchor="middle" font-family="monospace">vs.</text>';
 } else if(type === 'gpu'||type === 'matrix'||type==='systems') {
   for(let y=0;y<5;y++)for(let x=0;x<12;x++){let active=type==='matrix'?Math.abs(x-y-3)<2:(x+y)%4!==0;shapes+=`<rect x="${100+x*32}" y="${29+y*30}" width="24" height="22" rx="3" fill="${active?'#284845':'#161f22'}" stroke="${active?'#508580':'#263432'}"/><path d="M${107+x*32} ${37+y*30}h10" stroke="${active?'#8cbeb5':'#344b49'}"/>`;}
 } else if(type==='benchmark') {
   for(let i=0;i<9;i++)shapes+=`<rect x="${105+i*43}" y="${160-[38,70,53,103,86,121,99,135,124][i]}" width="23" height="${[38,70,53,103,86,121,99,135,124][i]}" rx="2" fill="${i===7?'#e9b77c':'#604d39'}"/><path d="M${116+i*43} ${153-[38,70,53,103,86,121,99,135,124][i]}v14" stroke="#e9b77c"/>`;
 } else if(type==='wave'||type==='physics') {
   for(let k=0;k<3;k++){let d='';for(let x=0;x<=600;x+=3)d+=`${x?'L':'M'}${x} ${100+Math.sin(x*.027+k*.8)*(45-k*8)}`;shapes+=`<path d="${d}" stroke="${['#8dbbf8','#456382','#304557'][k]}" fill="none" stroke-width="${k?1:2}"/>`;}
 } else {
   for(let i=0;i<6;i++)shapes+=`<ellipse cx="300" cy="105" rx="${35+i*27}" ry="${13+i*11}" fill="none" stroke="#677b4d" transform="rotate(-20 300 105)" opacity="${1-i*.1}"/>`;
   shapes+='<path d="M390 65 361 91 329 97 311 108 300 105" fill="none" stroke="#c4f77a" stroke-width="2"/><circle cx="300" cy="105" r="4" fill="#c4f77a"/>';
 }
 return `<svg class="card-art" viewBox="0 0 600 200" aria-hidden="true"><defs><pattern id="${key}" width="24" height="24" patternUnits="userSpaceOnUse"><path d="M24 0H0V24" fill="none" stroke="#a8bfac" stroke-opacity=".07" stroke-width=".6"/></pattern></defs><rect width="600" height="200" fill="url(#${key})"/>${shapes}</svg>`;
}
export function header(active='') {
 return `<a class="site-skip" href="#main">跳到主要内容</a><header class="lab-header"><div class="site-width header-row"><a class="lab-brand" href="/" aria-label="FENG 知识实验室首页"><span class="brand-glyph">f<span>.</span></span><span>FENG<span class="brand-divider">/</span><span class="brand-descriptor">知识实验室</span></span></a><nav class="desktop-nav" aria-label="主导航">${[['/','探索','home'],['/blog/','知识库','blog'],['/categories/','分类','categories'],['/archive/','时间线','archive'],['/about/','关于','about']].map(([href,title,id])=>`<a href="${href}" ${active===id?'aria-current="page"':''}>${title}</a>`).join('')}</nav><button class="search-trigger" aria-label="搜索知识库"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" aria-hidden="true"><circle cx="10.5" cy="10.5" r="6.5"/><path d="m16 16 5 5"/></svg><span>搜索知识</span><kbd>⌘ K</kbd></button><a class="github-link" href="https://github.com/DapengFeng" target="_blank" rel="noopener noreferrer" aria-label="GitHub（新窗口）"><svg viewBox="0 0 24 24" fill="currentColor" aria-hidden="true"><path d="M12 2a10 10 0 0 0-3.16 19.49c.5.09.68-.22.68-.48v-1.86c-2.78.61-3.37-1.18-3.37-1.18-.45-1.16-1.11-1.47-1.11-1.47-.91-.62.07-.61.07-.61 1 .07 1.53 1.03 1.53 1.03.89 1.53 2.34 1.09 2.91.83.09-.65.35-1.09.64-1.34-2.22-.25-4.56-1.11-4.56-4.94 0-1.09.39-1.99 1.03-2.69-.1-.25-.45-1.27.1-2.65 0 0 .84-.27 2.75 1.03A9.6 9.6 0 0 1 12 6.82c.85 0 1.7.11 2.5.34 1.91-1.3 2.75-1.03 2.75-1.03.55 1.38.2 2.4.1 2.65.64.7 1.03 1.6 1.03 2.69 0 3.84-2.34 4.69-4.57 4.94.36.31.68.92.68 1.85v2.75c0 .27.18.58.69.48A10 10 0 0 0 12 2Z"/></svg></a><button class="mobile-menu" aria-label="展开导航" aria-expanded="false">☰</button></div></header>`;
}
export function footer() { return `<footer class="lab-footer site-width"><div><a class="footer-name" href="/">FENG<span> / 知识实验室</span></a><p>保持好奇，把复杂的事想明白。</p></div><div class="footer-links"><a href="/archive/">文章归档</a><a href="/feed.xml">RSS</a><a href="https://github.com/DapengFeng">GitHub</a><span>© ${new Date().getFullYear()} ${pair(site.author,site.authorZh)}</span></div></footer>`; }
export function shell({title, description=site.description, body, active='', url='/', extraHead='', extraScripts=''}) {
 const token=process.env.CLOUDFLARE_WEB_ANALYTICS_TOKEN?.trim();
 if(token&&!/^[a-f0-9]{32}$/i.test(token))throw Error('CLOUDFLARE_WEB_ANALYTICS_TOKEN must be the 32-character token from the Cloudflare Web Analytics snippet, not a variable name or the full snippet.');
 const analytics=token?`<script type="module" src="/assets/analytics.js" data-analytics-token="${e(token)}" data-analytics-host="${e(new URL(site.url).hostname)}"></script>`:'';
 const values={TITLE:e(title),DESCRIPTION:e(description),CANONICAL:site.url+url,TYPE:url.endsWith('.html')?'article':'website',SITE_URL:site.url,EXTRA_HEAD:extraHead,HEADER:header(active),BODY:body,FOOTER:footer(),EXTRA_SCRIPTS:extraScripts,ANALYTICS:analytics};
 return layoutTemplate.replace(/\{\{([A-Z_]+)\}\}/g,(_,key)=>values[key]??'');
}
export function dateMeta(p) {
 const en=p.minutesEn??p.minutes,zh=p.minutesZh??p.minutes,both=p.minutesBoth??p.minutes;
 return `<span>${p.dateKind==='added'?'收录':'发布'} <time datetime="${p.date}">${p.date.replaceAll('-','.')}</time></span><span data-reading-minutes data-minutes-en="${en}" data-minutes-zh="${zh}" data-minutes-both="${both}">${pair(`${both} min read`,`${both} 分钟阅读`)}</span>`;
}
export function card(p, featured=false) {
 const c=categories.find(c=>c.id===p.category);
 return `<article class="knowledge-card ${featured?'featured-card':''}" data-category="${p.category}" data-date="${p.date}" data-title="${e(p.title)}" data-title-en="${e(p.titleEn||p.title)}" data-title-zh="${e(p.title)}" data-search="${e([p.title,p.titleEn,p.description,p.descriptionEn,p.cardDescription,p.cardDescriptionEn,...p.tags].join(' ').toLowerCase())}">
  <a class="card-link" href="${p.url}">
   <div class="card-visual visual-${p.art}">${p.cover?`<img class="card-photo" src="${e(p.cover)}" alt="" loading="lazy" decoding="async" width="640" height="480">`:art(p.art,p.slug+(featured?'-featured':''))}</div>
   <div class="card-copy">
    <div class="card-category" style="--category-color:${c.color}"><i></i>${pair(c.en,c.name)}</div>
    <h3>${pair(p.titleEn||p.title,p.title)}</h3>
    <p>${pair(p.cardDescriptionEn||p.descriptionEn||p.description,p.cardDescription||p.description)}</p>
    ${p.readingReason?`<p class="card-connection">${pair(p.readingReason.en,p.readingReason.zh)}</p>`:''}
    <div class="card-meta">${dateMeta(p)}</div>
   </div>
  </a>
 </article>`;
}
export function categoryCards(posts) {return `<div class="category-grid">${categories.map((c,i)=>`<a class="category-tile" style="--category-color:${c.color}" href="/blog/?category=${c.id}"><div class="category-top"><span class="category-symbol">${c.symbol}</span><span class="category-count">${String(posts.filter(p=>p.category===c.id&&!p.archiveOnly).length).padStart(2,'0')} 篇</span></div><h3>${c.name}</h3><p>${c.description}</p><span class="category-en">0${i+1} / ${c.en}</span></a>`).join('')}</div>`;}
export function filters(posts) {return `<div class="library-controls"><div class="filter-tabs" role="group" aria-label="按知识分类筛选"><button class="active" data-category-filter="all" aria-pressed="true">全部 <span>${posts.length}</span></button>${categories.map(c=>`<button data-category-filter="${c.id}" aria-pressed="false">${c.name}</button>`).join('')}</div><label class="sort-label"><span class="sr-only">文章排序</span><select id="article-sort"><option value="newest">最新发布</option><option value="oldest">最早发布</option><option value="title">标题 A–Z</option></select></label></div>`;}
export function home(posts,scene,date) {
 const visible=posts.filter(p=>!p.archiveOnly);
 const choices=['pytorch-01-what-is-pytorch','human-visual-system','chaoshan-streets-and-sea'];
 const selected=choices.map(slug=>visible.find(p=>p.slug===slug)).filter(Boolean);
 const latest=visible.filter(p=>!choices.includes(p.slug)).slice(0,4);
 const series=collectSeries(posts)[0];
 const routes=[
  {url:series?.url||'/blog/?category=systems',symbol:'⌘',en:'Inside PyTorch',zh:'沿着 PyTorch 深入',descriptionEn:'From tensor storage to dispatch, gradients and compiled kernels.',description:'从张量存储到算子调度、梯度与编译内核。'},
  {url:'/math/',symbol:'∿',en:'One mathematical idea',zh:'每天理解一个数学想法',descriptionEn:'A drawing, a relation, and the principle that connects them.',description:'从一幅图、一个关系式，读懂它们背后的原理。'},
  {url:'/blog/?category=travel',symbol:'⌁',en:'Places & observations',zh:'在路上的观察',descriptionEn:'Places, everyday life, and moments worth keeping in words and photographs.',description:'用照片与文字，留下一个地方的风景、日常与片刻。'}
 ];
 return shell({title:'让知识，变得可见',active:'home',body:`
 <main id="main" class="home-editorial">
  <section class="home-math-stage" data-daily-math aria-labelledby="home-title">
   <div class="home-hero site-width">
    <div class="hero-copy"><div class="overline home-overline"><span class="status-dot"></span>${pair('A PERSONAL KNOWLEDGE LAB','个人知识实验室')}<span class="hero-serial" data-current-month>VOL. ${new Intl.DateTimeFormat('en',{month:'numeric',timeZone:'Asia/Shanghai'}).format(new Date()).padStart(3,'0')}</span></div><h1 id="home-title">让知识，<br><em>变得可见。</em></h1></div>
    <div class="hero-profile">
     <div class="home-intro"><p class="home-intro-name">${pair(editorial.mission.title.en,editorial.mission.title.zh)}</p><p class="hero-description">${pair(editorial.mission.summary.en,editorial.mission.summary.zh)}</p></div>
     <div class="hero-actions"><a class="lime-button" href="#featured">${pair('Start with a note','从一篇笔记开始')}</a><a class="quiet-link" href="/blog/">${pair('Browse the notebook','浏览知识库')}</a></div>
    </div>
   </div>
   <div class="daily-math-background"><canvas id="surface-canvas" aria-hidden="true"></canvas><aside class="daily-math-note site-width" aria-labelledby="daily-math-title">${mathNote(scene)}</aside></div>${mathInitial(scene,date)}
  </section>
  <section class="site-width featured-section" id="featured" aria-labelledby="featured-title">
   <div class="section-heading"><h2 id="featured-title">${pair('Three places to begin','从这里读起')}</h2><a class="quiet-link" href="/blog/">${pair('All notes','全部文章')}</a></div>
   <div class="featured-grid">${selected.map(p=>card(p,true)).join('')}</div>
  </section>
  <section class="site-width category-section home-routes" aria-labelledby="explore-title">
   <div class="section-heading"><h2 id="explore-title">${pair('Follow a thread','沿着一个问题深入')}</h2><a class="quiet-link" href="/categories/">${pair('All topics','全部分类')}</a></div>
   <nav class="reading-routes" aria-label="Reading paths / 阅读路径">${routes.map(route=>`<a href="${route.url}"><span class="route-symbol" aria-hidden="true">${route.symbol}</span><h3>${pair(route.en,route.zh)}</h3><p>${pair(route.descriptionEn,route.description)}</p></a>`).join('')}</nav>
  </section>
  <section class="site-width recent-section" aria-labelledby="recent-title">
   <div class="section-heading"><h2 id="recent-title">${pair('Recent notes','最近更新')}</h2><a class="quiet-link" href="/archive/">${pair('Timeline','时间线')}</a></div>
   <div class="recent-notes">${latest.map(p=>`<a href="${p.url}"><time datetime="${p.date}">${p.date.replaceAll('-','.')}</time><div><h3>${pair(p.titleEn,p.title)}</h3><p>${pair(p.cardDescriptionEn||p.descriptionEn,p.cardDescription||p.description)}</p></div></a>`).join('')}</div>
  </section>
  <div class="site-width home-signoff"><p>${pair('Derivations, experiments, and moments worth keeping.','推导、实验，还有值得留下的片刻。')}</p><a href="/about/">${pair('About Dapeng','关于冯大鹏')}</a></div>
 </main>`,extraHead:'<link rel="stylesheet" href="/assets/math.css"><link rel="stylesheet" href="/assets/daily-math.css">',extraScripts:'<script type="module" src="/assets/surface.js"></script>'});
}
export function library(posts) {
 const visible=posts.filter(p=>!p.archiveOnly);
 return shell({title:'知识库',url:'/blog/',active:'blog',body:`<main id="main" class="site-width"><div class="page-heading"><span class="overline">THE NOTEBOOK / ${visible.length} NOTES</span><h1>知识库<span>.</span></h1><p>从一个问题开始，在公式、图形与实验之间找到答案。</p></div><label class="inline-search"><span>⌕</span><input type="search" id="library-search" placeholder="筛选标题、摘要与标签…" aria-label="筛选知识库"></label>${seriesLinks(posts)}${filters(visible)}<div class="article-grid" data-library>${visible.map(p=>card(p)).join('')}</div><div class="empty-state" hidden><span>∅</span><h3>没有找到匹配的笔记</h3><p>试试其他关键词，或换一个探索方向。</p><button data-reset-filters>清除筛选</button></div><div class="library-footer"><span id="result-count" aria-live="polite">共 ${visible.length} 篇笔记</span><a href="/archive/">完整归档</a></div></main>`});
}
export function categoryPage(posts) {return shell({title:'知识分类',url:'/categories/',active:'categories',body:`<main id="main" class="site-width"><div class="page-heading"><span class="overline">MAP OF KNOWLEDGE</span><h1>不同方向，无限联结<span>.</span></h1><p>建立索引，也发现知识之间的联系。</p></div>${categoryCards(posts)}${categories.map(c=>`<section class="category-group" id="${c.id}"><div class="section-heading"><h2 style="color:${c.color}"><span aria-hidden="true" class="category-symbol">${c.symbol}</span> ${pair(c.en,c.name)}</h2><a href="/blog/?category=${c.id}">浏览分类</a></div>${posts.filter(p=>p.category===c.id&&!p.archiveOnly).map(p=>`<a class="index-row" href="${p.url}"><time>${p.date}</time><h3>${e(p.title)}</h3><span>${p.tags.map(e).join(' / ')}</span></a>`).join('')}</section>`).join('')}</main>`});}
export function archive(posts) {return shell({title:'时间线',url:'/archive/',active:'archive',body:`<main id="main" class="site-width narrow-page"><div class="page-heading"><span class="overline">A GROWING COLLECTION / ${posts.length} NOTES</span><h1>想法的时间线<span>.</span></h1><p>${pair('Revisit each exploration by its first publication date.','按首次分享时间，回看每一次探索。')}</p></div>${[...new Set(posts.map(p=>p.date.slice(0,4)))].map(year=>`<section class="timeline-year"><h2>${year}<span>${posts.filter(p=>p.date.startsWith(year)).length} NOTES</span></h2><div>${posts.filter(p=>p.date.startsWith(year)).map(p=>`<a class="timeline-entry" href="${p.url}"><time datetime="${p.date}">${p.date.slice(5).replace('-','.')}<small>${p.dateKind==='added'?'收录':'发布'}</small></time><div><span class="timeline-category">${categories.find(c=>c.id===p.category).name}${p.archiveOnly?' / 早期记录':''}</span><h3>${e(p.title)}</h3><p>${e(p.description)}</p></div></a>`).join('')}</div></section>`).join('')}</main>`});}
export function about() {
 return shell({title:'关于实验室',url:'/about/',active:'about',body:`
 <main id="main" class="site-width narrow-page about-page">
  <div class="page-heading"><span class="overline">${pair('BEHIND THE NOTEBOOK','笔记背后')}</span><h1>${pair(editorial.mission.title.en,editorial.mission.title.zh)}</h1><p>${pair('A personal knowledge lab by Dapeng Feng.','冯大鹏的个人知识实验室。')}</p></div>
  <div class="about-grid"><div class="about-mark" aria-hidden="true">f<span>.</span><small>FENG</small></div><div class="prose">
   <h2>${pair('Why I write','我为什么写')}</h2><p>${pair(editorial.about.mission.en,editorial.about.mission.zh)}</p>
   <h2>${pair('How I work through an idea','怎样把一个问题想明白')}</h2><p>${pair(editorial.about.approach.en,editorial.about.approach.zh)}</p>
   <h2>${pair('A few places to begin','从这些文章开始')}</h2><ul class="about-reading"><li><a href="/series/pytorch-internals/">${pair('Trace one operation through PyTorch','跟踪一次 PyTorch 运算')}</a></li><li><a href="/blog/human-visual-system.html">${pair('Follow a red cup from light to neural signals','从红杯的光，读到神经信号')}</a></li><li><a href="/blog/chaoshan-streets-and-sea.html">${pair('Walk through Chaoshan, in photographs and words','在照片与文字中走过潮汕')}</a></li></ul>
   <h2>${pair('Sources and conversation','来源与交流')}</h2><p>${pair('Technical notes link to their sources and record assumptions and reproduction conditions. Questions and corrections belong with the article, where the explanation can be improved in context.','技术笔记链接原始来源，记录假设与复现条件。文章下的讨论区可以提问、纠错，把解释中尚未说清的地方继续补全。')}</p>
   <a class="lime-button" href="https://github.com/DapengFeng">${pair('Find me on GitHub','在 GitHub 找到我')}</a><a class="quiet-link about-rss" href="/feed.xml">${pair('Subscribe via RSS','通过 RSS 订阅')}</a>
  </div></div>
 </main>`});
}
export function article(p, posts, related = []) {
 const support=loadSupport();
 const c=categories.find(c=>c.id===p.category);
 const toc=p.headings.filter(h=>h.level===2);
 related=curatedRelated(p,posts,related);
 return shell({title:p.title,description:p.description,url:p.url,active:'blog',extraHead:`<link rel="stylesheet" href="/assets/math.css">${p.styles?`<style>${p.styles}</style><link rel="stylesheet" href="/assets/legacy.css">`:''}<link rel="stylesheet" href="/assets/reader.css"><link rel="stylesheet" href="/assets/discussions.css"><link rel="stylesheet" href="/assets/sharing.css"><link rel="stylesheet" href="/assets/support.css">`,body:`<div class="article-progress" aria-hidden="true"></div><main id="main" class="site-width article-page article-${p.category==='travel'?'travel':'technical'}"><div class="breadcrumbs"><a href="/blog/">知识库</a><span>/</span><a href="/blog/?category=${c.id}">${c.name}</a><span>/</span><span>知识笔记</span></div><header class="article-heading"><div class="overline" style="color:${c.color}">${c.en} <span> / </span> ${p.isHtml?'VISUAL ESSAY':'FIELD NOTES'}</div><h1>${e(p.title)}</h1><p>${e(p.description)}</p><div class="article-byline"><span class="author-avatar">F</span><a class="article-author" rel="author" href="/about/">${pair(site.author,site.authorZh)}</a><i></i>${dateMeta(p)}${p.updated?`<span>更新 <time datetime="${p.updated}">${p.updated}</time></span>`:''}</div><div class="article-tags">${p.tags.map(t=>`<a href="/blog/?q=${encodeURIComponent(t)}"># ${e(t)}</a>`).join('')}</div>${p.dateKind==='added'?'<p class="date-note">原文未标注发布日期；以上为本站收录日期。</p>':''}</header>${seriesNavigation(p,posts)}<div class="reading-layout"><aside class="article-toc"><details open><summary>本篇目录 <span>CONTENTS</span></summary><nav aria-label="文章目录">${toc.map((h,i)=>`<a href="#${e(h.id)}" ${h.titleEn?'':`data-lang="${h.lang||'en'}"`}><span>${String(toc.slice(0,i+1).filter(x=>x.lang===h.lang).length).padStart(2,'0')}</span>${h.titleEn?pair(h.titleEn,h.titleZh):e(h.title)}</a>`).join('')||'<span class="muted">短篇笔记</span>'}</nav></details><div class="toc-bottom"><span id="read-percentage">READING · 0%</span><a href="#main">回到顶部</a></div></aside><div class="article-reading-main"><article class="article-body " id="article-content">${renderSeriesPreviews(p,posts)}</article>${articleSharing(p)}${supportDialog(support)}</div></div>${seriesDirectory(p,posts)}${related.length?`<section class="related-section"><div class="section-heading"><h2>继续探索<span>CONNECTED IDEAS</span></h2><a href="/blog/">返回知识库</a></div><div class="featured-grid">${related.map(x=>card(x)).join('')}</div></section>`:''}${discussions(p,shareTrigger()+(support.paypal.url?supportTrigger():''))}</main>`,extraScripts:`<script src="/assets/sharing.js" defer></script><script src="/assets/support.js" defer></script><script src="/assets/discussions.js" defer></script><script src="/assets/article.js" defer></script>${p.hasCode?'<script src="/assets/syntax.js" defer></script>':''}${p.lab?'<script src="/assets/labs.js" defer></script>':''}${p.isPaired?'<script src="/assets/paired.js" defer></script>':''}${p.hasCompiler?'<script src="/assets/compiler.js" defer></script>':''}${p.hasCode?'<script src="/assets/code-copy.js" defer></script>':''}${p.hasReadingDemos?'<script type="module" src="/assets/reading-demos.js"></script>':''}`});
}

export function seriesLinks(posts) {
 const groups=collectSeries(posts);if(!groups.length)return '';
 return `<nav class="series-links" aria-label="Learning series / 学习专题">${groups.map(g=>`<a href="${g.url}"><span class="series-label">${pair('LEARNING SERIES','学习专题')}</span><strong>${pair(g.titleEn,g.title)}</strong><span>${pair(`${g.posts.length} published · Read in order`, `已发布 ${g.posts.length} 期 · 按顺序阅读`)}</span></a>`).join('')}</nav>`;
}
export function seriesPage(group){
 const [titleEn,title,descriptionEn,description]=seriesMetadata(group);
 return shell({title,url:group.url,active:'blog',body:`<main id="main" class="site-width"><div class="page-heading"><span class="overline">${pair('LEARNING SERIES','学习专题')}</span><h1>${pair(titleEn,title)}</h1><p>${pair(descriptionEn,description)}</p></div>${seriesGuide(group)}<ol class="series-episodes">${group.posts.map(p=>`<li id="part-${p.series.part}" value="${p.series.part}"><a href="${p.url}"><span class="series-label">${pair(`PART ${String(p.series.part).padStart(2,'0')}`,`第 ${p.series.part} 期`)} · <time datetime="${p.date}">${p.date}</time></span><h2>${pair(p.titleEn,p.title)}</h2><p>${pair(p.descriptionEn,p.description)}</p></a></li>`).join('')}</ol><p><a href="/blog/">${pair('Browse all notes','浏览全部笔记')}</a></p></main>`});
}
export function seriesNavigation(post,posts){
 if(!post.series)return '';
 const group=collectSeries(posts).find(g=>g.id===post.series.id),index=group.posts.findIndex(p=>p.slug===post.slug),previous=group.posts[index-1],next=group.posts[index+1];
 return `<div class="series-navigation"><nav class="article-series" aria-label="Series navigation / 专题导航"><a href="${group.url}">${pair(group.titleEn,group.title)}</a><span>${pair(`Part ${post.series.part}`,`第 ${post.series.part} 期`)}</span><div>${previous?`<a rel="prev" href="${previous.url}">${pair('Previous part','上一期')}</a>`:''}${next?`<a rel="next" href="${next.url}">${pair('Next part','下一期')}</a>`:''}</div></nav></div>`;
}

export function seriesDirectory(post,posts){
 if(!post.series)return '';
 return `<div class="series-footer"><details class="series-directory" open><summary>${pair('Series contents','专题目录')}</summary>${renderSeriesPreviews({...post,html:'<div data-series-preview></div>'},posts)}</details></div>`;
}
