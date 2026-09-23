import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import {mathStyles} from './math.mjs';
import {loadContent} from './content.mjs';
import {site,categories,escape as e} from './config.mjs';
import * as templates from './templates.mjs';
import {localizePage,dictionary} from './i18n.mjs';
import {optimizePage,renderSharingImage,sharingImage} from './seo.mjs';
export async function build(){
 const posts=await loadContent();
 await fs.rm('dist',{recursive:true,force:true});await fs.mkdir('dist/assets',{recursive:true});
 async function write(file,content){await fs.mkdir(path.dirname('dist/'+file),{recursive:true});await fs.writeFile('dist/'+file,content);}
 const pages=[['index.html',templates.home(posts)],['blog/index.html',templates.library(posts)],['categories/index.html',templates.categoryPage(posts)],['archive/index.html',templates.archive(posts)],['about/index.html',templates.about()],['404.html',templates.shell({title:'未找到页面',body:'<main id="main" class="site-width page-heading"><span class="overline">404 / UNCHARTED TERRITORY</span><h1>这里还没有留下笔记。</h1><a class="lime-button" href="/">返回首页 →</a></main>'})]];
 for(const p of posts)pages.push([p.url.slice(1),templates.article(p,posts)]);
 for(const [file,html]of pages){
  const url=file==='index.html'?'/':'/'+file.replace(/index\.html$/, '');
  await write(file,optimizePage(localizePage(html,posts),url,posts.find(post=>post.url===url)));
 }
 for(const post of [null,...posts])await write(sharingImage(post).slice(1),await renderSharingImage(post));
 for(const name of ['site','legacy','reader'])await fs.copyFile(`src/styles/${name}.css`,`dist/assets/${name}.css`);
 for(const name of ['site','surface','article','labs','benchmark-worker','paired','syntax','compiler'])await fs.copyFile(`src/scripts/${name}.js`,`dist/assets/${name}.js`);
 // MathJax + AMS renders self-contained SVGs at build time, with no browser runtime.
 await write('assets/math.css',mathStyles());
 await fs.copyFile('node_modules/@mathjax/src/LICENSE','dist/assets/mathjax-LICENSE');
 await write('assets/math-NOTICE','MathJax and MathJax-Newcm font, version 4.1.3.\nCopyright MathJax Consortium. Licensed under Apache-2.0; see mathjax-LICENSE.\nhttps://github.com/mathjax/MathJax-src\nhttps://github.com/mathjax/MathJax-fonts\n');
 await write('assets/favicon.svg','<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64"><rect width="64" height="64" rx="13" fill="#c0f47b"/><text x="19" y="48" font-family="Georgia" font-style="italic" font-size="53" font-weight="bold" fill="#101310">f.</text></svg>');
 const index=posts.map(({html,styles,headings,isHtml,source,...p})=>({...p,categoryEn:categories.find(c=>c.id===p.category).en,categoryZh:categories.find(c=>c.id===p.category).name,tags:[...p.tags,...p.tags.map(t=>dictionary[t]||t)]}));
 await write('search-index.json',JSON.stringify(index));
 const items=posts.map(p=>`<item><title>${e(p.titleEn||p.title)} / ${e(p.title)}</title><link>${site.url}${p.url}</link><guid isPermaLink="true">${site.url}${p.url}</guid><pubDate>${new Date(p.date+'T12:00:00+08:00').toUTCString()}</pubDate><description>${e(p.descriptionEn||p.description)}</description><category>${e(p.category)}</category></item>`).join('');
 await write('feed.xml',`<?xml version="1.0" encoding="UTF-8"?><rss version="2.0"><channel><title>FENG / Knowledge Lab</title><link>${site.url}</link><description>${e(site.description)}</description><language>en</language>${items}</channel></rss>`);
 await write('sitemap.xml',`<?xml version="1.0" encoding="UTF-8"?><urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">${['/','/blog/','/categories/','/archive/','/about/',...posts.map(p=>p.url)].map(url=>`<url><loc>${site.url}${url}</loc>${posts.find(p=>p.url===url)?`<lastmod>${posts.find(p=>p.url===url).updated||posts.find(p=>p.url===url).date}</lastmod>`:''}</url>`).join('')}</urlset>`);
 await write('robots.txt',`User-agent: *\nAllow: /\nSitemap: ${site.url}/sitemap.xml\n`);await write('.nojekyll','');
 // Preserve historical Jekyll permalinks and old section URLs.
 const redirects=[['search/index.html','/blog/'],['publications/index.html','/about/']];
 for(const p of posts.filter(p=>p.date<'2026-01-01')){const date=p.date.replaceAll('-','/');redirects.push([`${date}/${p.slug}/index.html`,p.url]);}
 for(const [file,to]of redirects)await write(file,`<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="refresh" content="0;url=${to}"><link rel="canonical" href="${site.url}${to}"><title>Moved · FENG</title><a href="${to}">Continue / 继续阅读 →</a></html>`);
 console.log(`Built ${posts.length} bilingual notes and ${pages.length} pages → dist/`);return posts;
}
if(process.argv[1]===fileURLToPath(import.meta.url))await build();
