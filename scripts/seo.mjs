import * as cheerio from 'cheerio';
import sharp from 'sharp';
import {site, categories, escape as e} from './config.mjs';

const pageInfo = {
  '/': ['Ideas, made visible', '让知识，变得可见', 'Explore mathematics, physics, programming, and benchmarks through visual explanations and interactive experiments by Dapeng Feng.', '冯大鹏的个人知识实验室：通过图解、交互实验与基准测试，探索数学、物理和编程。'],
  '/blog/': ['Knowledge notebook', '知识库', 'Browse bilingual notes on mathematics, physics, systems programming, and reproducible benchmarks.', '浏览数学、物理、系统编程与可复现基准测试的中英双语笔记。'],
  '/categories/': ['Knowledge categories', '知识分类', 'Find visual notes by subject: mathematics and algorithms, physics and models, systems and performance, and benchmarks.', '按主题查找可视化笔记：数学与算法、物理与模型、系统与性能、基准测试。'],
  '/archive/': ['Article timeline', '文章时间线', 'Explore all notes by their first publication date, from earlier algorithm notes to interactive visual essays.', '按首次发布日期浏览全部笔记，从早期算法记录到交互式可视化长文。'],
  '/about/': ['About Dapeng Feng', '关于冯大鹏', 'Meet Dapeng Feng and learn about this personal knowledge lab for mathematics, physics, programming, and performance experiments.', '了解冯大鹏，以及这个分享数学、物理、编程和性能实验的个人知识实验室。'],
  '/404.html': ['Page not found', '页面未找到', 'This page could not be found. Browse the knowledge notebook to continue exploring.', '未找到此页面，可前往知识库继续探索。'],
};
const bilingual = (en, zh) => en === zh ? en : `${en} / ${zh}`;
const absolute = path => site.url + path;
export const sharingImage = post => `/assets/social/${post ? post.slug : 'site'}.png`;
const json = value => JSON.stringify(value).replace(/</g, '\\u003c');

// One canonical bilingual document; no fictitious hreflang alternatives for UI modes.
export function optimizePage(html, url, post, pageMetadata) {
  const $ = cheerio.load(html);
  const info = post ? [post.titleEn, post.title, post.descriptionEn, post.description] : pageMetadata || pageInfo[url];
  if (!info || info.some(value => typeof value !== 'string' || !value.trim())) throw Error(`Missing bilingual search metadata: ${url}`);
  const [titleEn, titleZh, descriptionEn, descriptionZh] = info;
  const title = bilingual(titleEn, titleZh), description = bilingual(descriptionEn, descriptionZh);
  const canonical = absolute(url), image = absolute(sharingImage(post));
  function meta(attribute, key, value) {
    $(`meta[${attribute}="${key}"]`).remove();
    $('head').append($('<meta>').attr(attribute, key).attr('content', value));
  }
  $('title').text(`${title} · FENG`).attr('data-en', `${titleEn} · FENG`).attr('data-zh', `${titleZh} · FENG`);
  $('link[rel="canonical"]').attr('href', canonical);
  meta('name', 'description', description);
  meta('name', 'author', `${site.author} / ${site.authorZh}`);
  meta('name', 'robots', url === '/404.html' ? 'noindex,follow' : 'index,follow,max-image-preview:large');
  for (const [name, value] of [['google-site-verification', process.env.GOOGLE_SITE_VERIFICATION], ['msvalidate.01', process.env.BING_SITE_VERIFICATION]]) {
    if (value?.trim()) meta('name', name, value.trim());
  }
  for (const [key, value] of Object.entries({title, description, url: canonical, type: post ? 'article' : 'website', site_name: 'FENG / Knowledge Lab', image, 'image:type': 'image/png', 'image:width': '1200', 'image:height': '630', 'image:alt': title})) meta('property', `og:${key}`, value);
  for (const [key, value] of Object.entries({card: 'summary_large_image', title, description, image, 'image:alt': title})) meta('name', `twitter:${key}`, value);
  if (post) {
    meta('property', 'article:published_time', post.date);
    meta('property', 'article:modified_time', post.updated || post.date);
    meta('property', 'article:author', absolute('/about/'));
    $('.article-byline>span').filter((_, node) => $(node).text().includes(site.author)).wrapInner('<a href="/about/" rel="author"></a>');
  }

  const personId = absolute('/about/#person'), websiteId = absolute('/#website'), pageId = `${canonical}#webpage`;
  const graph = [
    {'@type': 'Person', '@id': personId, name: site.author, alternateName: site.authorZh, url: absolute('/about/'), sameAs: ['https://github.com/DapengFeng']},
    {'@type': 'WebSite', '@id': websiteId, url: absolute('/'), name: 'FENG / Knowledge Lab', alternateName: 'FENG / 知识实验室', inLanguage: ['en', 'zh-CN'], author: {'@id': personId}},
    {'@type': url === '/about/' ? 'AboutPage' : post || url === '/' || url === '/404.html' ? 'WebPage' : 'CollectionPage', '@id': pageId, url: canonical, name: title, description, inLanguage: ['en', 'zh-CN'], isPartOf: {'@id': websiteId}, primaryImageOfPage: {'@type': 'ImageObject', url: image, width: 1200, height: 630}},
  ];
  if (url === '/about/') graph[2].mainEntity = {'@id': personId};
  if (url !== '/' && url !== '/404.html') {
    const trail = [{name: 'Home / 首页', item: absolute('/')}];
    if (post) trail.push({name: 'Knowledge notebook / 知识库', item: absolute('/blog/')});
    trail.push({name: title, item: canonical});
    graph[2].breadcrumb = {'@id': `${canonical}#breadcrumb`};
    graph.push({'@type': 'BreadcrumbList', '@id': `${canonical}#breadcrumb`, itemListElement: trail.map((item, i) => ({'@type': 'ListItem', position: i + 1, ...item}))});
  }
  if (post) {
    // Cite only actual links explicitly placed in reference sections or marked by the author.
    let sourceLinks = $('.article-body #references a[href], .article-body .sources-grid a[href], .article-body a[data-citation]');
    $('.article-body h2#references').each((_, heading) => {
      sourceLinks = sourceLinks.add($(heading).nextUntil('h2').find('a[href]'));
    });
    const citations = [...new Set(sourceLinks.map((_, node) => $(node).attr('href')).get().filter(href => /^https?:\/\//.test(href)))];
    const section = categories.find(category => category.id === post.category);
    graph[2].mainEntity = {'@id': `${canonical}#article`};
    graph.push({'@type': 'BlogPosting', '@id': `${canonical}#article`, url: canonical, headline: title, description, inLanguage: ['en', 'zh-CN'], author: {'@id': personId}, publisher: {'@id': personId}, mainEntityOfPage: {'@id': pageId}, datePublished: post.date, dateModified: post.updated || post.date, image: [image], articleSection: bilingual(section.en, section.name), keywords: post.tags, ...(citations.length ? {citation: citations} : {}), hasPart: post.headings.filter(h => h.level === 2).map(h => ({'@type': 'WebPageElement', name: h.titleEn ? bilingual(h.titleEn, h.titleZh || h.titleEn) : h.title, url: `${canonical}#${encodeURIComponent(h.id)}`}))});
  }
  $('script[type="application/ld+json"]').remove();
  $('head').append(`<script type="application/ld+json">${json({'@context': 'https://schema.org', '@graph': graph})}</script>`);
  return $.html();
}

// Build raster previews because many sharing clients do not display SVG cards.
// The English title is graphical branding; the surrounding metadata is bilingual.
export async function renderSharingImage(post) {
  const title = post?.titleEn || 'Ideas, made visible.';
  const words = title.split(/\s+/), lines = [];
  for (const word of words) {
    if (!lines.length || `${lines.at(-1)} ${word}`.length > 34) lines.push(word);
    else lines[lines.length - 1] += ` ${word}`;
  }
  const size = lines.length > 3 ? 45 : 58;
  const step = size + 16;
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="630"><rect width="1200" height="630" fill="#101310"/><path d="M72 124H1128M72 520H1128" stroke="#43513d"/><circle cx="1110" cy="79" r="10" fill="#c0f47b"/><g font-family="DejaVu Sans,Arial,sans-serif"><text x="72" y="88" fill="#c0f47b" font-size="26">FENG / KNOWLEDGE LAB</text>${lines.map((line, i) => `<text x="72" y="${215 + i * step}" fill="#f0f2ea" font-size="${size}" font-weight="bold">${e(line)}</text>`).join('')}<text x="72" y="483" fill="#b2bba9" font-size="24">${e(post ? `${post.date}   /   ${post.category.toUpperCase()}` : 'MATHEMATICS / PHYSICS / SYSTEMS / BENCHMARKS')}</text><text x="72" y="573" fill="#b2bba9" font-size="25">Dapeng Feng</text><text x="1128" y="573" text-anchor="end" fill="#c0f47b" font-size="25">dapengfeng.github.io</text></g></svg>`;
  return sharp(Buffer.from(svg)).png().toBuffer();
}
