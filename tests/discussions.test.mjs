import test from 'node:test';
import assert from 'node:assert/strict';
import {load} from 'cheerio';
import {discussions,discussionTerm} from '../scripts/discussions.mjs';
test('article discussions use one stable article identity and escaped metadata',()=>{
 const p={url:'/blog/example.html',titleEn:'<script>alert(1)</script>',title:'知识 " < >'};
 const $=load(discussions(p));
 assert.equal($('script').length,0);
 assert.equal($('#article-discussions').attr('data-description'),p.titleEn+' / '+p.title);
 assert.equal($('.discussion-thread').length,1);
 assert.equal($('[data-discussion-kind]').length,0);
 assert.equal(discussionTerm(p.url),p.url+' · comment','preserve the original Comments thread');
 assert.equal($('#discussion-comment').attr('data-term'),discussionTerm(p.url));
 assert.ok($('#discussion-comment').attr('data-category-id'));
 assert.equal($('.discussion-thread:not([hidden])').length,1);
 assert.equal($('.discussion-launcher').attr('aria-expanded'),'false');
 assert.equal($('iframe').length,0,'no eager third-party request');
 assert.ok($('noscript').text());
});
