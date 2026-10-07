import {site, escape as e} from './config.mjs';
import {pair} from './i18n.mjs';

// Public GitHub node IDs, not credentials. Terms depend on URLs, never titles or UI language.
export const discussionRepo = 'DapengFeng/dapengfeng.github.io';
export const discussionRepoId = 'MDEwOlJlcG9zaXRvcnkyNDkzMTU1OTU=';
const category = 'General';
const categoryId = 'MDE4OkRpc2N1c3Npb25DYXRlZ29yeTMyMDY4NTky';
const bubble = 'M4 4h16v12H9l-5 4V4Z';
// Keep the original Comments identity so existing article threads stay connected.
export const discussionTerm = url => `${url} · comment`;
const icon = d => `<svg viewBox="0 0 24 24" aria-hidden="true" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round"><path d="${d}"/></svg>`;
export function discussions(post,extraActions='') {
 return `<section class="article-discussions" id="article-discussions" aria-labelledby="discussions-title" data-repo="${discussionRepo}" data-repo-id="${discussionRepoId}" data-backlink="${e(site.url+post.url)}" data-description="${e(post.titleEn+' / '+post.title)}">
 <header class="discussion-heading"><h2 id="discussions-title">${pair('Discussion','讨论')}</h2></header>
 <div class="discussion-shell" id="discussion-shell">
 <div class="discussion-toolbar"><span class="discussion-panel-title">${icon(bubble)}${pair('Discussion','讨论',true)}</span><button type="button" class="discussion-close" hidden data-aria-label-en="Close floating discussion" data-aria-label-zh="收起浮动讨论" aria-label="Close floating discussion / 收起浮动讨论">${icon('m6 6 12 12M6 18 18 6')}</button></div>
 <div class="discussion-thread" id="discussion-comment" data-term="${e(discussionTerm(post.url))}" data-category="${category}" data-category-id="${categoryId}">
 <p class="discussion-status" role="status" hidden>${pair('Loading GitHub discussion…','正在加载 GitHub 讨论…')}</p>
 <div class="discussion-embed"></div>
 <button type="button" class="discussion-retry" hidden aria-label="Reload discussion / 重新加载讨论" data-aria-label-en="Reload discussion" data-aria-label-zh="重新加载讨论" title="Reload discussion / 重新加载讨论" data-title-en="Reload discussion" data-title-zh="重新加载讨论">${icon('M20 7v5h-5M20 12a8 8 0 1 0-2.3 5.7')}</button>
 </div>
 </div>
 <noscript><p>${pair('Enable JavaScript to read and join the discussion below.','请启用 JavaScript，在文章下方查看和参与讨论。')}</p></noscript>
 <div class="article-action-dock" role="group" aria-label="Article actions / 文章操作" data-aria-label-en="Article actions" data-aria-label-zh="文章操作">${extraActions}<button type="button" class="discussion-launcher" aria-controls="discussion-shell" aria-expanded="false" aria-label="Open discussion / 打开讨论" data-aria-label-en="Open discussion" data-aria-label-zh="打开讨论" title="Open discussion / 打开讨论" data-title-en="Open discussion" data-title-zh="打开讨论" hidden>${icon(bubble)}</button></div>
 </section>`;
}
