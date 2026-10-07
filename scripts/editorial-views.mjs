import {readFileSync} from 'node:fs';
import {pair} from './i18n.mjs';
import {escape as e} from './config.mjs';

export const editorial = JSON.parse(readFileSync(new URL('../content/editorial.json', import.meta.url), 'utf8'));
const bilingual = value => pair(value.en, value.zh);

export function curatedRelated(post, posts, fallback = []) {
  const links = editorial.related.filter(link => link.source === post.slug);
  const curated = links.map(link => {
    const target = posts.find(candidate => candidate.slug === link.target);
    if (!target || target.slug === post.slug) throw Error(`Invalid editorial connection: ${link.source} / ${link.target}`);
    return {...target, readingReason: link.reason};
  });
  const seen = new Set([post.slug]);
  return [...curated, ...fallback].filter(item => {
    if (seen.has(item.slug)) return false;
    seen.add(item.slug);
    return true;
  }).slice(0, 3);
}

export function seriesGuide(group) {
  const guide = editorial.series.find(item => item.id === group.id);
  if (!guide) return '';
  const facts = [['For whom', '适合谁', guide.audience], ['Before you begin', '已有基础', guide.prerequisites], ['What you can explain', '读完能解释什么', guide.outcomes]];
  return `<section class="series-guide" aria-label="Reading guide / 阅读指南">
    <dl>${facts.map(([en, zh, value]) => `<div><dt>${pair(en, zh)}</dt><dd>${bilingual(value)}</dd></div>`).join('')}</dl>
    <div class="series-paths">${guide.paths.map(route => `<div class="series-path">
      <h2>${bilingual(route.title)}</h2><p>${bilingual(route.description)}</p>
      <ol>${route.parts.map(part => {
        const post = group.posts.find(item => item.series.part === part);
        if (!post) throw Error(`Unknown reading path part: ${group.id}/${part}`);
        return `<li><a href="${e(post.url)}">${pair(`Part ${part}`, `第 ${part} 期`)}</a></li>`;
      }).join('')}</ol>
    </div>`).join('')}</div>
  </section>`;
}
