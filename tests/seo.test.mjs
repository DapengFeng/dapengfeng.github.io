import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import * as cheerio from 'cheerio';
import sharp from 'sharp';
import {site} from '../scripts/config.mjs';
import {optimizePage} from '../scripts/seo.mjs';

const files = (await fs.readdir('dist', {recursive: true})).filter(file => file.endsWith('.html'));
const pages = [];
for (const file of files) {
  const $ = cheerio.load(await fs.readFile(`dist/${file}`, 'utf8'));
  if (!$('meta[http-equiv="refresh"]').length) pages.push({file, $});
}

test('every published page has unique bilingual metadata and a real raster sharing image', async () => {
  const titles = new Set(), descriptions = new Set();
  assert.ok(pages.length >= 16);
  for (const {file, $} of pages) {
    const title = $('title').text(), description = $('meta[name="description"]').attr('content');
    assert.match(title, /[A-Za-z].*\/.*[\u3400-\u9fff]/u, file);
    assert.match(description, /[A-Za-z].*\/.*[\u3400-\u9fff]/u, file);
    assert.ok(!titles.has(title), `Duplicate title: ${file}`); titles.add(title);
    assert.ok(!descriptions.has(description), `Duplicate description: ${file}`); descriptions.add(description);
    assert.equal($('link[rel="canonical"]').length, 1, file);
    const url = site.url + (file === 'index.html' ? '/' : '/' + file.replace(/index\.html$/, ''));
    assert.equal($('link[rel="canonical"]').attr('href'), url);
    assert.equal($('meta[property="og:url"]').attr('content'), url);
    assert.equal($('meta[property="og:description"]').attr('content'), description);
    assert.equal($('meta[name="twitter:card"]').attr('content'), 'summary_large_image');
    const image = new URL($('meta[property="og:image"]').attr('content'));
    assert.equal(image.origin, site.url);
    const data = await sharp(await fs.readFile(`dist${image.pathname}`)).metadata();
    assert.equal(data.format, 'png'); assert.equal(data.width, 1200); assert.equal(data.height, 630);
    assert.equal($('link[hreflang]').length, 0, 'Reading modes are not independent language URLs');
  }
});

test('article structured data describes visible authors, dates, references and real section anchors', () => {
  for (const {file, $} of pages.filter(page => page.$('.article-body').length)) {
    assert.equal($('script[type="application/ld+json"]').length, 1);
    const graph = JSON.parse($('script[type="application/ld+json"]').text())['@graph'];
    const post = graph.find(item => item['@type'] === 'BlogPosting');
    const person = graph.find(item => item['@type'] === 'Person');
    assert.equal(post.author['@id'], person['@id']);
    assert.equal(person.name, site.author); assert.equal(person.alternateName, site.authorZh);
    assert.equal($('.article-byline a[rel="author"]').attr('href'), '/about/');
    assert.equal(post.datePublished, $('.article-byline time').first().attr('datetime'));
    assert.equal(post.dateModified, $('.article-byline time').last().attr('datetime'));
    assert.deepEqual(post.inLanguage, ['en', 'zh-CN']);
    assert.ok($('.article-body [data-lang="en"]').length && $('.article-body [data-lang="zh"]').length, file);
    const ids = new Set($('[id]').map((_, node) => $(node).attr('id')).get());
    for (const section of post.hasPart) assert.ok(ids.has(decodeURIComponent(new URL(section.url).hash.slice(1))), section.url);
    const links = new Set($('.article-body a[href]').map((_, node) => $(node).attr('href')).get());
    for (const citation of post.citation || []) assert.ok(links.has(citation), `Invented citation: ${citation}`);
    if (file === 'blog/matrix-multiplication.html') assert.ok(post.citation.includes('https://www.netlib.org/blas/'));
  }
});

test('404 is excluded from indexing and sitemap, and metadata is safely escaped', async () => {
  const $ = pages.find(page => page.file === '404.html').$;
  assert.equal($('meta[name="robots"]').attr('content'), 'noindex,follow');
  assert.equal($('meta[property="og:type"]').attr('content'), 'website');
  assert.ok(!(await fs.readFile('dist/sitemap.xml', 'utf8')).includes('/404.html'));
  const dangerous = 'A </script><script>alert(1)</script> & "title"';
  const result = cheerio.load(optimizePage('<html><head><title></title><link rel="canonical"></head><body></body></html>', '/blog/test.html', {slug: 'test', titleEn: dangerous, title: '测试', descriptionEn: 'Test', description: '测试描述', category: 'math', date: '2026-09-23', tags: [], headings: []}));
  assert.equal(result('script').length, 1);
  assert.ok(JSON.parse(result('script').text())['@graph'].find(item => item['@type'] === 'BlogPosting').headline.includes(dangerous));
});
