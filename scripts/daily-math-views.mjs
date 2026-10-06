import {pair} from './i18n.mjs';
import {escape as e} from './config.mjs';
import {renderDisplayMath,renderInlineMathText} from './math.mjs';
export const mathHead='<link rel="stylesheet" href="/assets/math.css"><link rel="stylesheet" href="/assets/daily-math.css">';
export const mathScript='<script type="module" src="/assets/surface.js"></script>';
const richPair=(en,zh)=>`<span class="i18n"><span data-lang="en">${renderInlineMathText(en)}</span><span data-lang="zh">${renderInlineMathText(zh)}</span></span>`;
export const emptyMathNote=()=>`<div class="daily-math-identity"><h2 id="daily-math-title">${pair('Mathematics archive','数学往期')}</h2><a class="daily-math-source" href="/math/">${pair('Browse the collection','浏览往期内容')}</a></div>`;
export function mathNote(t,{archive=false}={}){
 if(!t)return emptyMathNote();
 return `<div class="daily-math-identity"><div class="daily-math-kicker"><span>${pair('THE MATHEMATICS','背后的数学')}</span><time data-math-date datetime="${t.date}">${t.date.replaceAll('-','.')}</time></div><${archive?'h1':'h2'} id="daily-math-title" data-math-title>${pair(t.en,t.zh)}</${archive?'h1':'h2'}><div class="daily-math-formula" data-math-formula><div data-math-expression="${t.id}" data-latex="${e(t.formula)}"${t.compactFormula?' data-compact-formula':''}><div class="daily-math-full">${renderDisplayMath(t.formula)}</div>${t.compactFormula?`<div class="daily-math-compact">${renderDisplayMath(t.compactFormula)}</div>`:''}</div></div></div><div class="daily-math-explanation"><div data-math-copy="${t.id}"><p data-math-description>${richPair(t.descriptionEn,t.descriptionZh)}</p><p class="daily-math-reading" data-math-reading>${richPair(t.readingEn,t.readingZh)}</p></div><div class="daily-math-links"><a class="daily-math-source" data-math-source href="${e(t.source.url)}">${pair(t.source.en,t.source.zh)}</a><a class="daily-math-source" href="/math/">${pair('Archive','往期')}</a><a class="daily-math-source" href="/math/${t.date}/">${pair('Permalink','日期链接')}</a></div></div>`;
}
export function mathInitial(t,date){return `<script type="application/json" data-math-entry>${JSON.stringify(t?{date:t.date,id:t.id,renderer:t.renderer}:{date}).replaceAll('<','\\u003c')}</script>`;}
export function mathPage(t){return `<main id="main"><div class="site-width math-page-heading"><a href="/math/">${pair('Daily mathematics','每日数学')}</a></div><section class="home-math-stage math-entry-stage" data-daily-math data-math-fixed-date="${t.date}" aria-labelledby="daily-math-title"><div class="daily-math-background"><canvas id="surface-canvas" aria-hidden="true"></canvas><aside class="daily-math-note site-width" aria-labelledby="daily-math-title">${mathNote(t,{archive:true})}</aside></div>${mathInitial(t,t.date)}</section></main>`;}
export function mathArchive(entries,today,year){
 const years=[...new Set(entries.map(t=>t.date.slice(0,4)))].sort().reverse();
 const rows=year?entries.filter(t=>t.date.startsWith(year)):entries;
 // Year pages stay bounded to at most 366 rows. The landing page carries only recent
 // history plus this year's prepared queue; older entries remain on their year page.
 const listed=year?rows:[...rows.filter(t=>t.date<=today).slice(-30),...rows.filter(t=>t.date>today).slice(0,31)];
 return `<main id="main" class="site-width math-archive" data-math-archive${year?' data-year-archive':''}><div class="page-heading"><span class="overline">MATHEMATICS / ${year||'ARCHIVE'}</span><h1>${pair(year?`Mathematics · ${year}`:'Daily mathematics',year?`${year} 年每日数学`:'每日数学')}</h1><p>${pair('A mathematical idea, its equation, and a moving diagram.','一个数学原理，一组公式，一幅动态图。')}</p></div><nav class="math-years" aria-label="Years / 年份"><a href="/math/">${pair('Recent','最近')}</a>${years.map(y=>`<a data-math-year="${y}" href="/math/${y}/"${y>today.slice(0,4)?' hidden':''}>${y}</a>`).join('')}</nav><ol class="math-archive-list">${[...listed].reverse().map(t=>`<li data-entry-date="${t.date}"${t.date>today?' hidden':''}><a href="/math/${t.date}/"><time datetime="${t.date}">${t.date}</time><span>${pair(t.en,t.zh)}</span></a></li>`).join('')}</ol><p data-math-archive-empty${listed.some(t=>t.date<=today)?' hidden':''}>${pair('No entries published yet.','暂时没有已发布内容。')}</p></main>`;
}
export function mathMetadata(t){return [t.en,t.zh,t.descriptionEn.replace(/\\\(|\\\)/g,''),t.descriptionZh.replace(/\\\(|\\\)/g,'')];}
