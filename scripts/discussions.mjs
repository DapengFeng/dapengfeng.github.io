import {site, escape as e} from './config.mjs';
import {pair} from './i18n.mjs';

// Public GitHub node IDs, not credentials. Terms depend on URLs, never titles or UI language.
export const discussionRepo = 'DapengFeng/dapengfeng.github.io';
export const discussionRepoId = 'MDEwOlJlcG9zaXRvcnkyNDkzMTU1OTU=';
export const discussionKinds = [
 {id:'comment', en:'Comments', zh:'评论', category:'General', categoryId:'MDE4OkRpc2N1c3Npb25DYXRlZ29yeTMyMDY4NTky', icon:'M4 4h16v12H9l-5 4V4Z'},
 {id:'idea', en:'Ideas', zh:'想法', category:'Ideas', categoryId:'MDE4OkRpc2N1c3Npb25DYXRlZ29yeTMyMDY4NTk0', icon:'M9 18h6m-5 3h4M8 14a7 7 0 1 1 8 0l-1 2H9l-1-2Z'},
 {id:'discussion', en:'Discussion', zh:'讨论', category:'General', categoryId:'MDE4OkRpc2N1c3Npb25DYXRlZ29yeTMyMDY4NTky', icon:'M3 3h13v10H7l-4 3V3Zm6 14h8l4 3V8h-2'}
];
export const discussionTerm = (url,kind) => `${url} · ${kind}`;
const icon = d => `<svg viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"><path d="${d}"/></svg>`;
export function discussions(post) {
 const root=`https://github.com/${discussionRepo}/discussions`;
 return `<section class="article-discussions" id="article-discussions" aria-labelledby="discussions-title" data-repo="${discussionRepo}" data-repo-id="${discussionRepoId}" data-backlink="${e(site.url+post.url)}" data-description="${e(post.titleEn+' / '+post.title)}">
 <header class="discussion-heading"><h2 id="discussions-title">${pair('Discuss this article','讨论这篇文章')}</h2><a href="${root}" target="_blank" rel="noopener noreferrer">GitHub ↗</a></header>
 <div class="discussion-shell" id="discussion-shell">
 <div class="discussion-toolbar"><div class="discussion-kinds" role="group" aria-label="Discussion type / 讨论类型">${discussionKinds.map((k,i)=>`<button type="button" data-discussion-kind="${k.id}" aria-pressed="${!i}" aria-controls="discussion-${k.id}" aria-label="${k.en} / ${k.zh}" data-aria-label-en="${k.en}" data-aria-label-zh="${k.zh}" title="${k.en} / ${k.zh}" data-title-en="${k.en}" data-title-zh="${k.zh}">${icon(k.icon)}${pair(k.en,k.zh,true)}</button>`).join('')}</div><button type="button" class="discussion-close" hidden data-aria-label-en="Close floating discussion" data-aria-label-zh="收起浮动讨论" aria-label="Close floating discussion / 收起浮动讨论">${icon('m6 6 12 12M6 18 18 6')}</button></div>
 ${discussionKinds.map((k,i)=>`<div class="discussion-thread" id="discussion-${k.id}" ${i?'hidden':''} data-kind="${k.id}" data-term="${e(discussionTerm(post.url,k.id))}" data-category="${k.category}" data-category-id="${k.categoryId}">
 <p class="discussion-status" role="status" hidden>${pair('Loading GitHub discussion…','正在加载 GitHub 讨论…')}</p>
 <div class="discussion-embed"></div>
 <button type="button" class="discussion-retry" hidden>${pair('Reload','重新加载',true)}</button>
 </div>`).join('')}
 </div>
 <noscript><p>${pair('JavaScript is required for embedded comments. You can also join the discussion on GitHub.','内嵌评论需要 JavaScript，也可前往 GitHub 参与讨论。')}</p></noscript>
 <button type="button" class="discussion-launcher" aria-controls="discussion-shell" aria-expanded="false" hidden>${icon(discussionKinds[0].icon)}<span>${pair('Join the discussion…','参与讨论…',true)}</span>${icon('m9 5 7 7-7 7')}</button>
 </section>`;
}
