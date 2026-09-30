import test from 'node:test';
import assert from 'node:assert/strict';
import {load} from 'cheerio';
import {discussions,discussionTerm,discussionKinds} from '../scripts/discussions.mjs';
test('article discussions use stable, distinct URL/type identities and escaped metadata',()=>{
 const p={url:'/blog/example.html',titleEn:'<script>alert(1)</script>',title:'知识 " < >'};
 const $=load(discussions(p));
 assert.equal($('script').length,0);
 assert.equal($('#article-discussions').attr('data-description'),p.titleEn+' / '+p.title);
 assert.equal($('.discussion-thread').length,3);
 assert.equal(new Set(discussionKinds.map(k=>discussionTerm(p.url,k.id))).size,3);
 for(const k of discussionKinds){
  assert.equal($('#discussion-'+k.id).attr('data-term'),discussionTerm(p.url,k.id));
  assert.ok($('#discussion-'+k.id).attr('data-category-id'));
 }
 assert.equal($('.discussion-thread:not([hidden])').length,1);
 assert.equal($('.discussion-launcher').attr('aria-expanded'),'false');
 assert.equal($('iframe').length,0,'no eager third-party request');
 assert.ok($('noscript').text());
});
