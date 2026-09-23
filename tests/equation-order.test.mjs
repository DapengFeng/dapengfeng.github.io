import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {load} from 'cheerio';
import {parseContent} from '../scripts/content.mjs';
import {article} from '../scripts/templates.mjs';

// Validate against the published contents, independently of the renderer's counter.
function audit(html, name) {
 const $=load(html),root=$('#article-content'),chapters=new Map();
 $('.article-toc nav a[href^="#"]').each((i,el)=>{
  const link=$(el),number=Number(link.children('span').first().text());
  assert.equal(number,i+1,`${name}: contents numbering`);
  chapters.set(decodeURIComponent(link.attr('href').slice(1)),number);
 });
 let chapter=0,total=0;
 const counts=new Map(),seen=new Set();
 root.find('*').each((_,el)=>{
  const node=$(el),id=node.attr('id');
  if(chapters.has(id))chapter=chapters.get(id);
  if(!node.hasClass('formula-block'))return;
  assert.ok(chapter>0,`${name}: display formula must belong to a contents chapter`);
  const labels=node.find('.formula-tools .eq-number');
  assert.equal(labels.length,1,`${name}: exactly one display number`);
  const sequence=(counts.get(chapter)||0)+1,expected=`${chapter}.${sequence}`,actual=labels.text();
  assert.equal(actual,expected,`${name}: formula order must follow contents and reading order`);
  assert.ok(!seen.has(actual),`${name}: duplicate ${actual}`);
  assert.equal(node.parents('[data-lang]').length,0,`${name}: share display formulas across languages`);
  counts.set(chapter,sequence);seen.add(actual);total++;
 });
 assert.equal(root.find('.formula-inline .eq-number').length,0,`${name}: inline formulas have no numbers`);
 return total;
}

test('all generated pages have unique, continuous formula numbers matching their contents',async()=>{
 const files=(await fs.readdir('dist',{recursive:true})).filter(f=>f.endsWith('.html'));
 let total=0;
 for(const file of files)total+=audit(await fs.readFile(`dist/${file}`,'utf8'),file);
 assert.ok(total>=99,'audit covers all published display equations');
});

test('the authoring example uses the same equation numbering rules as published articles',async()=>{
 const post=await parseContent('examples/article.html');
 assert.equal(audit(article(post,[post]),'examples/article.html'),3);
});

test('separate English and Chinese headings count as one chapter',async()=>{
 const post=await parseContent('content/posts/waves-and-phase.html');
 const $=load(post.html);
 assert.deepEqual($('.eq-number').map((_,el)=>$(el).text()).get(),['1.1','3.1']);
 assert.equal(audit(article(post,[post]),post.slug),2);
});
