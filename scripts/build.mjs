import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import {renderShareAssets} from './sharing-assets.mjs';
import {loadDailyMath} from './daily-math.mjs';
import {dateInfo} from '../src/scripts/daily-math.js';
import {mathNote,mathPage,mathArchive,mathMetadata,mathHead,mathScript} from './daily-math-views.mjs';
import {mathStyles} from './math.mjs';
import {loadContent} from './content.mjs';
import {site,categories,escape as e} from './config.mjs';
import * as templates from './templates.mjs';
import {collectSeries,seriesMetadata} from './series.mjs';
import {createRelatedRecommender} from './related.mjs';
import {localizePage,dictionary} from './i18n.mjs';
import {optimizePage,renderSharingImage,sharingImage} from './seo.mjs';
export async function build(){
 const {entries}=await loadDailyMath(),today=dateInfo().date;
 const posts=await loadContent(),series=collectSeries(posts),recommend=createRelatedRecommender(posts);
 await fs.rm('dist',{recursive:true,force:true});await fs.mkdir('dist/assets',{recursive:true});
 await fs.cp('content/assets','dist/assets/content',{recursive:true});
 async function write(file,content){await fs.mkdir(path.dirname('dist/'+file),{recursive:true});await fs.writeFile('dist/'+file,content);}
 const pages=[['index.html',templates.home(posts,entries.find(t=>t.date===today),today)],['blog/index.html',templates.library(posts)],['categories/index.html',templates.categoryPage(posts)],['archive/index.html',templates.archive(posts)],['about/index.html',templates.about()],['404.html',templates.shell({title:'未找到页面',body:'<main id="main" class="site-width page-heading"><span class="overline">404 / UNCHARTED TERRITORY</span><h1>这里还没有留下笔记。</h1><a class="lime-button" href="/">返回首页</a></main>'})]];
 for(const group of series)pages.push([group.url.slice(1)+'index.html',templates.seriesPage(group)]);
 for(const p of posts)pages.push([p.url.slice(1),templates.article(p,posts,recommend(p))]);
 const mathMeta=new Map(),years=[...new Set(entries.map(t=>t.date.slice(0,4)))];
 for(const year of ['',...years]){
  const url=`/math/${year?year+'/':''}`,meta=[year?`Mathematics · ${year}`:'Daily mathematics',year?`${year} 年每日数学`:'每日数学',`Mathematical principles, equations and animated diagrams${year?' from '+year:''}.`,`数学原理、公式与动态图${year?'：'+year+' 年往期':'的日期归档'}。`];
  pages.push([url.slice(1)+'index.html',templates.shell({title:meta[0],url,body:mathArchive(entries,today,year),extraHead:mathHead,extraScripts:'<script type="module" src="/assets/math-archive.js"></script>'})]);mathMeta.set(url,meta);
 }
 for(const t of entries){
  const url=`/math/${t.date}/`;
  pages.push([url.slice(1)+'index.html',templates.shell({title:t.en,url,body:mathPage(t),extraHead:mathHead,extraScripts:mathScript})]);
  mathMeta.set(url,mathMetadata(t));
  await write(`assets/daily-math/${t.date}.json`,JSON.stringify({date:t.date,id:t.id,renderer:t.renderer,html:mathNote(t)}));
 }
 for(const [file,html]of pages){
  const url=file==='index.html'?'/':'/'+file.replace(/index\.html$/, '');
  const group=series.find(group=>group.url===url);
  let output=optimizePage(localizePage(html,posts),url,posts.find(post=>post.url===url),group?seriesMetadata(group):mathMeta.get(url));
  const future=/^\/math\/(\d{4}-\d{2}-\d{2})\/$/.exec(url);
  if(future&&future[1]>today)output=output.replace(/(<meta name="robots" content=")[^"]*/, '$1noindex,follow');
  await write(file,output);
 }
 for(const post of posts)await renderShareAssets(post,write);
 for(const post of [null,...posts])await write(sharingImage(post).slice(1),await renderSharingImage(post));
 for(const name of ['editorial','daily-math','site','legacy','reader','discussions','giscus-theme','sharing','support'])await fs.copyFile(`src/styles/${name}.css`,`dist/assets/${name}.css`);
 for(const name of ['analytics','daily-math','daily-math-models','daily-math-drawings','math-archive','eye-viewer','eye-renderer','eye-model','sharing','support','discussions','site','surface','article','reading-demos','process-demos','process-models','labs','benchmark-worker','paired','syntax','compiler','code-copy'])await fs.copyFile(`src/scripts/${name}.js`,`dist/assets/${name}.js`);
 for(const name of ['three.module.min.js','three.core.min.js'])await fs.copyFile(`node_modules/three/build/${name}`,`dist/assets/${name}`);
 await fs.copyFile('node_modules/three/LICENSE','dist/assets/three-LICENSE');
 // MathJax + AMS renders self-contained SVGs at build time, with no browser runtime.
 await write('assets/math.css',mathStyles());
 await fs.copyFile('node_modules/@mathjax/src/LICENSE','dist/assets/mathjax-LICENSE');
 await write('assets/math-NOTICE','MathJax and MathJax-Newcm font, version 4.1.3.\nCopyright MathJax Consortium. Licensed under Apache-2.0; see mathjax-LICENSE.\nhttps://github.com/mathjax/MathJax-src\nhttps://github.com/mathjax/MathJax-fonts\n');
 await write('assets/favicon.svg','<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 64 64"><rect width="64" height="64" rx="13" fill="#c0f47b"/><text x="19" y="48" font-family="Georgia" font-style="italic" font-size="53" font-weight="bold" fill="#101310">f.</text></svg>');
 // Search records contain only public discovery fields, not recommendation working data.
 const index=posts.map(p=>({kind:'article',slug:p.slug,url:p.url,title:p.title,titleEn:p.titleEn,description:p.cardDescription||p.description,descriptionEn:p.cardDescriptionEn||p.descriptionEn,date:p.date,category:p.category,categoryEn:categories.find(c=>c.id===p.category).en,categoryZh:categories.find(c=>c.id===p.category).name,tags:[...p.tags,...p.tags.map(t=>dictionary[t]||t)],searchTextEn:p.searchTextEn,searchTextZh:p.searchTextZh,featured:p.featured,archiveOnly:p.archiveOnly,minutes:p.minutes,cover:p.cover}));
 index.push(...entries.filter(t=>t.date<=today).map(t=>({kind:'math',url:`/math/${t.date}/`,title:t.zh,titleEn:t.en,description:t.descriptionZh,descriptionEn:t.descriptionEn,date:t.date,category:'math',categoryEn:'DAILY MATHEMATICS',categoryZh:'每日数学',tags:[t.id],searchTextEn:[t.en,t.descriptionEn,t.readingEn].join(' '),searchTextZh:[t.zh,t.descriptionZh,t.readingZh].join(' ')})));
 index.push(...series.map(g=>{const [titleEn,title,descriptionEn,description]=seriesMetadata(g);return {kind:'series',url:g.url,title,titleEn,description,descriptionEn,category:'systems',categoryEn:'LEARNING SERIES',categoryZh:'学习专题',tags:['PyTorch','系列','series'],searchText:g.posts.map(p=>p.title+' '+p.titleEn).join(' ')};}));
 await write('search-index.json',JSON.stringify(index));
 const items=posts.map(p=>`<item><title>${e(p.titleEn||p.title)} / ${e(p.title)}</title><link>${site.url}${p.url}</link><guid isPermaLink="true">${site.url}${p.url}</guid><pubDate>${new Date(p.date+'T12:00:00+08:00').toUTCString()}</pubDate><description>${e(p.descriptionEn||p.description)}</description><category>${e(p.category)}</category></item>`).join('');
 await write('feed.xml',`<?xml version="1.0" encoding="UTF-8"?><rss version="2.0"><channel><title>FENG / Knowledge Lab</title><link>${site.url}</link><description>${e(site.description)}</description><language>en</language>${items}</channel></rss>`);
 await write('sitemap.xml',`<?xml version="1.0" encoding="UTF-8"?><urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">${['/','/blog/','/categories/','/archive/','/about/','/math/',...years.filter(y=>y<=today.slice(0,4)).map(y=>`/math/${y}/`),...entries.filter(t=>t.date<=today).map(t=>`/math/${t.date}/`),...series.map(group=>group.url),...posts.map(p=>p.url)].map(url=>`<url><loc>${site.url}${url}</loc>${posts.find(p=>p.url===url)?`<lastmod>${posts.find(p=>p.url===url).updated||posts.find(p=>p.url===url).date}</lastmod>`:''}</url>`).join('')}</urlset>`);
 await write('robots.txt',`User-agent: *\nAllow: /\nSitemap: ${site.url}/sitemap.xml\n`);await write('.nojekyll','');
 // Preserve historical Jekyll permalinks and old section URLs.
 const redirects=[['search/index.html','/blog/'],['publications/index.html','/about/']];
 for(const p of posts.filter(p=>p.date<'2026-01-01')){const date=p.date.replaceAll('-','/');redirects.push([`${date}/${p.slug}/index.html`,p.url]);}
 for(const [file,to]of redirects)await write(file,`<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta http-equiv="refresh" content="0;url=${to}"><link rel="canonical" href="${site.url}${to}"><title>Moved · FENG</title><a href="${to}">Continue / 继续阅读</a></html>`);
 console.log(`Built ${posts.length} bilingual notes and ${pages.length} pages → dist/`);return posts;
}
if(process.argv[1]===fileURLToPath(import.meta.url))await build();
