import fs from 'node:fs/promises';
import path from 'node:path';
import * as cheerio from 'cheerio';
import {renderMath} from './math.mjs';
import postcss from 'postcss';
import prefixer from 'postcss-prefix-selector';
import { categories } from './config.mjs';

// Count the text a reader can see in each edition. Shared figures and code are
// included once; UI labels, hidden demo state, and metadata are not prose.
function readingText($, language) {
  const copy = $('body').clone();
  copy.find('script,style,nav,header,footer,button,input,select,textarea,svg,canvas,[hidden],[aria-hidden="true"],[data-legacy-chrome],[data-series-preview]').remove();
  if (language) copy.find(`[data-lang="${language === 'en' ? 'zh' : 'en'}"]`).remove();
  copy.find('br').replaceWith(' ');
  copy.find('p,div,section,h1,h2,h3,h4,h5,h6,li,td,th,figcaption,pre').append(' ');
  return copy.text().replace(/\s+/g, ' ').trim();
}

export function readingMinutes(text) {
  const chinese = text.match(/[\u3400-\u9fff]/g)?.length || 0;
  const words = text.match(/[a-z]+(?:['’_-][a-z]+)*/gi)?.length || 0;
  return Math.max(1, Math.ceil(chinese / 400 + words / 220));
}

export function addHeadings($) {
  const headings = [], used = new Set($('[id]').map((_, e) => $(e).attr('id')).get());
  const chapters = new Set();
  function anchor(node, fallback) {
    let id = node.attr('id');
    if (!id) { id = fallback; while (used.has(id)) id += '-x'; node.attr('id', id); }
    used.add(id);
    return id;
  }
  function languageText(node, language) {
    return node.find(`[data-lang="${language}"]`).first().text().trim();
  }
  $('h2,h3').each((i, el) => {
    const h = $(el);
    if (h.closest('header,footer,.hero,.rail,.site-header,[data-legacy-chrome]').length) return;
    const id = anchor(h, `section-${i + 1}`);
    // Imported visual articles group paired headings around shared experiments.
    const chapter = h.closest('.bilingual-chapter');
    if (chapter.length) {
      if (chapters.has(chapter[0])) return;
      chapters.add(chapter[0]);
      const titleEn = chapter.find('[data-lang="en"] h2, h2 [data-lang="en"]').first().text().trim();
      const titleZh = chapter.find('[data-lang="zh"] h2, h2 [data-lang="zh"]').first().text().trim();
      headings.push({id: anchor(chapter, `topic-${chapters.size}`), level: 2, titleEn, titleZh, title: titleEn || titleZh});
      return;
    }
    const titleEn = languageText(h, 'en') || h.attr('data-title-en');
    const titleZh = languageText(h, 'zh') || h.attr('data-title-zh');
    headings.push({id, titleEn, titleZh, title: h.text().replace(/\s+/g, ' ').trim(), level: el.tagName === 'h2' ? 2 : 3});
  });
  return headings;
}
export async function parseContent(file) {
  const raw = await fs.readFile(file, 'utf8');
  const name = path.basename(file).replace(/\.html$/, '');
  const slug = name.replace(/^\d{4}-\d{2}-\d{2}-/, '');
  const $ = cheerio.load(raw);
  let declared = {};
  const metadata = $('#article-metadata');
  if (metadata.length) { declared = JSON.parse(metadata.text()); metadata.remove(); }
  const isHtml = $('style').length > 0;
  const metaValue = key => $(`meta[name="${key}"],meta[property="${key}"]`).first().attr('content');
  const data = {
    title: $('title').text().split('|')[0].trim() || $('h1').first().text().trim(),
    description: metaValue('description'),
    date: metaValue('article:published_time')?.slice(0, 10) || name.match(/^\d{4}-\d{2}-\d{2}/)?.[0],
    category: metaValue('category'), tags: metaValue('tags')?.split(',').map(x => x.trim()),
    lang: $('html').attr('lang')?.startsWith('zh') ? 'zh' : 'en',
    ...declared,
  };
  if (data.draft) return null;
  if (data.date instanceof Date) data.date = data.date.toISOString().slice(0, 10);
  if (!data.title || !/^\d{4}-\d{2}-\d{2}$/.test(data.date || '') || new Date(data.date).toISOString().slice(0,10) !== data.date) throw new Error(`${file}: title and valid YYYY-MM-DD date are required. Use HTML meta or article-metadata JSON.`);
  if (!categories.some(c => c.id === data.category)) throw new Error(`${file}: category must be one of ${categories.map(c => c.id).join(', ')}.`);
  if (data.updated && (!/^\d{4}-\d{2}-\d{2}$/.test(data.updated) || data.updated < data.date)) throw new Error(`${file}: updated must be YYYY-MM-DD and not precede date.`);
  if (data.tags && !Array.isArray(data.tags)) throw new Error(`${file}: tags must be an array.`);
  // Both languages are authored in this file, in reading order. No translation lookup.
  const hasTranslation = $('body [data-lang="en"]').length > 0 && $('body [data-lang="zh"]').length > 0;
  if (!hasTranslation) throw new Error(`${file}: include both English and Chinese in this HTML using data-lang="en" and data-lang="zh".`);
  $('[data-lang="en"]').attr('lang', 'en');
  $('[data-lang="zh"]').attr('lang', 'zh-CN');
  const isPaired = $('#paired-terms,#article-terms').length > 0;
  const plain = readingText($);
  const plainEn = readingText($, 'en'), plainZh = readingText($, 'zh');
  const textCopy = $('body').clone(); textCopy.find('script,style,nav,header,footer,[hidden],[aria-hidden="true"]').remove();
  const relatedCopy = textCopy.clone();
  relatedCopy.find('pre,code,button,input,select,svg,[data-series-preview],[data-citation]').remove();
  const relatedText = relatedCopy.text().replace(/\s+/g, ' ').trim();
  const headings = addHeadings($);
  renderMath($, headings);
  const styles = [];
  if (isHtml) {
    for (const el of $('style').toArray()) {
      const result = await postcss([prefixer({ prefix: '.legacy-content', transform(prefix, selector, prefixed) {
        if (/^(:root|html|body)(\b|(?=[\s:[.#>]))/.test(selector)) return selector.replace(/^(:root|html|body)/, prefix).replace(/\.legacy-content\[data-theme="dark"\]/g, '.legacy-content');
        return prefixed;
      } })]).process($(el).html(), { from: undefined });
      styles.push(result.css); $(el).remove();
    }
    // Keep hidden DOM used by the original demonstrations; preserve their scripts.
    $('body > .hero, body > header, .site-header, .site-footer, .rail, main > .hero, .skip, .progress').attr('data-legacy-chrome', 'true');
    $('main').each((_,el) => {el.tagName = 'div';});
  }
  const headScripts = $('head script').toArray().map(e => $.html(e)).join('');
  let html = headScripts + $('body').html();
  if (!$('.article-edition').length) html = `<section class="article-edition ${isHtml ? 'legacy-content' : 'prose'}">${html}</section>`;
  return {
    ...data, slug, source: file, tags: data.tags || [], art: data.art || data.category,
    description: data.description || plain.slice(0, 130),
    hasTranslation, isPaired, hasCompiler: $('code[data-godbolt]').length > 0, lab: data.lab,
    hasCode: $('pre code,pre').length > 0,
    hasReadingDemos: $('#compile-route-lab,#cuda-streams-lab,#autograd-lab,#dispatch-lab,#matrix-lab,#fw-lab,#pt-route,#pt-multiply,#tensor-layout,.memory-demo,#band-lab,#symmetric-lab,#cuda-map,#lab-lif,#lab-ultra,#lab-bptt,#lab-stdp,#lab-event,#fft-lab,[data-wave-lab]').length > 0,
    minutes: readingMinutes(plain), minutesEn: readingMinutes(plainEn), minutesZh: readingMinutes(plainZh), minutesBoth: readingMinutes(plain),
    url: `/blog/${slug}.html`, html, styles: styles.join('\n'), headings, isHtml,
    relatedText, searchText: plain.slice(0, 65000), searchTextEn: plainEn.slice(0, 65000), searchTextZh: plainZh.slice(0, 65000), dateKind: data.dateKind || 'published',
  };
}
async function walk(dir) {
  const entries = await fs.readdir(dir, {withFileTypes:true});
  const nested = await Promise.all(entries.map(e => e.isDirectory() ? walk(path.join(dir,e.name)) : path.join(dir,e.name)));
  return nested.flat().filter(f => /\.html$/.test(f) && path.basename(f) !== 'index.html').sort();
}
export async function loadContent(directory='content/posts') {
  const files = await walk(directory);
  const posts = (await Promise.all(files.map(file => parseContent(file)))).filter(Boolean);
  const seen = new Set();
  for (const p of posts) { if (seen.has(p.slug)) throw new Error(`Duplicate article slug: ${p.slug}`); seen.add(p.slug); }
  return posts.sort((a,b) => b.date.localeCompare(a.date) || (a.featured || 99) - (b.featured || 99) || a.slug.localeCompare(b.slug));
}
