import {site,escape as e} from './config.mjs';
import {pair} from './i18n.mjs';

const svg=body=>`<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${body}</svg>`;
const icons={
 share:svg('<circle cx="18" cy="5" r="3"/><circle cx="6" cy="12" r="3"/><circle cx="18" cy="19" r="3"/><path d="m9 10 6-4M9 14l6 4"/>'),
 link:svg('<path d="m10 13 4-4M8 16l-1 1a4 4 0 0 1-6-6l4-4a4 4 0 0 1 6 0m2 1 1-1a4 4 0 0 1 6 6l-4 4a4 4 0 0 1-6 0"/>'),
 copy:svg('<rect x="8" y="8" width="12" height="13" rx="2"/><path d="M16 8V3H3v13h5"/>'),
 save:svg('<path d="M12 3v12m-5-5 5 5 5-5M4 16v5h16v-5"/>'),
 close:svg('<path d="m6 6 12 12M6 18 18 6"/>'),
 wechat:svg('<path d="M13 14c-1 .5-2 .8-3 .8L5 18l.5-4C3.3 12.8 2 11 2 9c0-3.3 3.5-6 8-6s8 2.7 8 6"/><path d="M22 14c0 2-1 3.3-2.5 4l.5 3-3.5-2H16c-3.3 0-6-2.2-6-5s2.7-5 6-5 6 2.2 6 5Z"/><path d="M7 8h.01M12 8h.01M14 13h.01M18 13h.01" stroke-width="2.7"/>'),
 rednote:svg('<rect x="3" y="3" width="18" height="18" rx="5"/><path d="M8 7v10m4-7v7m4-10v10"/>'),
 x:svg('<path d="M4 3h5l11 18h-5ZM20 3l-7 8M4 21l7-8"/>'),
 linkedin:svg('<rect x="2" y="2" width="20" height="20" rx="3"/><path d="M7 10v7m0-11v.01M11 17v-7m0 3c0-4 6-4 6 0v4"/>'),
 telegram:svg('<path d="m2 10 20-7-4 18-6-6-4 3v-6l10-6-10 6Z"/>'),
 more:svg('<circle cx="5" cy="12" r="1"/><circle cx="12" cy="12" r="1"/><circle cx="19" cy="12" r="1"/>')
};
const label=(en,zh)=>`aria-label="${e(en+' / '+zh)}" data-aria-label-en="${e(en)}" data-aria-label-zh="${e(zh)}" title="${e(en+' / '+zh)}" data-title-en="${e(en)}" data-title-zh="${e(zh)}"`;
export const sharePaths=post=>({qr:`/assets/share/${post.slug}-qr.png`,poster:mode=>`/assets/share/${post.slug}-${mode}.png`});
export function shareTrigger(){return `<button type="button" class="article-share-trigger" data-share-open aria-haspopup="dialog" aria-controls="article-share-dialog" ${label('Share article','分享文章')} hidden>${icons.share}</button>`;}
export function articleSharing(post){
 const paths=sharePaths(post),url=site.url+post.url;
 const action=(name,icon,en,zh)=>`<button type="button" data-share-action="${name}" ${label(en,zh)}>${icons[icon]}<span class="sr-only">${pair(en,zh)}</span></button>`;
 const platform=(id,en,zh)=>`<button type="button" data-share-platform="${id}" aria-pressed="false">${icons[id]}${pair(en,zh)}</button>`;
 const external=(id,title,href)=>`<a data-share-platform="${id}" href="${e(href)}" target="_blank" rel="noopener noreferrer">${icons[id]}<span>${title}</span></a>`;
 return `<div class="article-share-footer">${shareTrigger()}</div>
 <dialog class="article-share-dialog" id="article-share-dialog" aria-labelledby="article-share-title" data-share-url="${e(url)}" data-title-en="${e(post.titleEn)}" data-title-zh="${e(post.title)}" data-summary-en="${e(post.descriptionEn)}" data-summary-zh="${e(post.description)}" data-poster-base="/assets/share/${e(post.slug)}">
 <div class="share-heading"><h2 id="article-share-title">${pair('Share article','分享文章',true)}</h2><button type="button" data-share-close ${label('Close sharing','关闭分享')}>${icons.close}</button></div>
 <div class="share-platforms">${platform('wechat','WeChat','微信')}${platform('rednote','RedNote','小红书')}${external('x','X','https://x.com/intent/tweet?'+new URLSearchParams({url,text:post.titleEn}))}${external('linkedin','LinkedIn','https://www.linkedin.com/sharing/share-offsite/?'+new URLSearchParams({url}))}${external('telegram','Telegram','https://t.me/share/url?'+new URLSearchParams({url,text:post.titleEn}))}<button type="button" data-share-native hidden>${icons.more}${pair('More','更多',true)}</button></div>
 <div class="share-media" data-share-media hidden>
 <figure class="share-qr" data-share-qr hidden><img data-src="${paths.qr}" width="246" height="246" alt="Article QR code / 文章二维码" data-alt-en="Article QR code" data-alt-zh="文章二维码"><figcaption>${pair('Scan in WeChat, or save the image and share it in a chat or Moments.','用微信扫码，或保存图片后分享给好友、发布到朋友圈。')}</figcaption></figure>
 <figure class="share-poster" data-share-poster hidden><img data-src="${paths.poster('both')}" width="900" height="1200" alt="Article sharing image / 文章分享图" data-alt-en="Article sharing image" data-alt-zh="文章分享图"><figcaption>${pair('Save the image, then paste the text into a new RedNote post.','保存图片，再将文案粘贴到小红书的新笔记中。')}</figcaption></figure>
 </div>
 <label class="share-text-label" for="article-share-text">${pair('Sharing text','分享文案',true)}</label>
 <textarea id="article-share-text" rows="4" spellcheck="false"></textarea>
 <div class="share-actions">${action('copy-link','link','Copy link','复制链接')}${action('copy-text','copy','Copy text','复制文案')}<a data-share-download href="${paths.poster('both')}" download="${e(post.slug)}-both.png" ${label('Save sharing image','保存分享图')}>${icons.save}<span class="sr-only">${pair('Save sharing image','保存分享图')}</span></a><button type="button" data-share-file hidden ${label('Share image to an app','将图片分享到应用')}>${icons.share}</button><output role="status" aria-live="polite" data-share-status></output></div>
 <div class="share-manual" hidden><label for="share-manual-text">${pair('Select and copy','选中并复制',true)}</label><textarea id="share-manual-text" readonly rows="3"></textarea></div>
 </dialog>`;
}
