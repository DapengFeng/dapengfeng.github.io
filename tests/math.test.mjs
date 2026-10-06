import test from 'node:test';
import assert from 'node:assert/strict';
import {load} from 'cheerio';
import {renderMath,renderDisplayMath,renderInlineMathText} from '../scripts/math.mjs';
import {loadDailyMath} from '../scripts/daily-math.mjs';
const {topics:scenes}=await loadDailyMath();

function render(latex, inline=false) {
 const $=load('<div id="formula"></div>');
 $('#formula').attr('data-math',latex).attr('data-display',inline?'inline':'display');
 renderMath($);return $;
}
test('daily display equations retain their own glyphs across topic changes',()=>{
 const allIds=new Set();
 for(const latex of scenes.flatMap(scene=>[scene.formula,scene.compactFormula].filter(Boolean))){
  const $=load(renderDisplayMath(latex));
  assert.equal($('mjx-container[display="true"][role="math"]').attr('aria-label'),latex);
  assert.equal($('[data-mml-node="merror"],script').length,0);
  $('use').each((_,el)=>{
   const id=$(el).attr('href').slice(1);
   assert.equal($(`[id="${id}"]`).length,1,'glyphs resolve within this expression, including after a daily switch');
  });
  $('[id]').each((_,el)=>{
   const id=$(el).attr('id');assert.ok(!allIds.has(id),'expressions do not collide in the same page');allIds.add(id);
  });
 }
 assert.equal(render('a=b')('.formula-block').attr('data-equation-number'),'1.1','home expressions do not alter article numbering');
 assert.throws(()=>renderDisplayMath(String.raw`\unknownCommand{x}`),/Invalid LaTeX/);
});
test('daily inline mathematics preserves prose, escapes HTML, and uses text style',()=>{
 const sources=scenes.flatMap(scene=>[scene.descriptionEn,scene.descriptionZh,scene.readingEn,scene.readingZh]);
 sources.push(String.raw`a < b & <script>literal</script>: \(\frac{1}{2}\), then \(x^2\).`);
 for(const source of sources){
  const $=load(`<p>${renderInlineMathText(source)}</p>`);
  assert.equal($('script,mjx-container[display],[data-mml-node="merror"]').length,0);
  assert.equal($('mjx-container svg:not([aria-hidden="true"])').length,0,'each expression has one accessible label, without duplicate unlabeled glyph images');
  assert.deepEqual($('mjx-container').map((_,el)=>$(el).attr('aria-label')).get(),[...source.matchAll(/\\\(([\s\S]*?)\\\)/g)].map(match=>match[1]));
  $('mjx-container').each((_,el)=>{
   const formula=$(el);
   assert.equal(formula.children('svg').length,1,'short inline expressions stay intact at line boundaries');
   formula.find('use').each((_,glyph)=>{
    const id=$(glyph).attr('href').slice(1);
    assert.equal(formula.find(`[id="${id}"]`).length,1,'inline glyphs are self-contained');
   });
   formula.replaceWith($('<span></span>').text(`\\(${formula.attr('aria-label')}\\)`));
  });
  assert.equal($('p').text(),source,'surrounding words, punctuation, and math source survive typesetting');
 }
 assert.throws(()=>renderInlineMathText(String.raw`Broken \(\unknownCommand{x}\)`),/Invalid LaTeX/);
});
test('AMS environments and commands render to standalone SVG with exact copy source',()=>{
 for(const latex of [String.raw`\begin{align} a&=b+c\\d&=e \end{align}`,String.raw`\begin{gather} a=b\\c=d \end{gather}`,String.raw`f(x)=\begin{cases}x^2 & x\ge 0\\-x & x<0\end{cases}`,String.raw`A=\begin{pmatrix}1&2\\3&4\end{pmatrix},\quad x\in\mathbb{R},\quad\operatorname{rank}(A)=2`,String.raw`\begin{aligned}a&=b\\c&=d\end{aligned}\tag{1}`]){
  const $=render(latex);
  assert.equal($('.formula-block').attr('data-latex'),latex);
  assert.equal($('.formula-copy').length,1);
  assert.equal($('mjx-container > svg').length,1);
  assert.equal($('[data-mml-node="merror"]').length,0);
  assert.equal($('script').length,0);
  assert.ok($('mjx-container use').length>0);
  $('mjx-container use').each((_,element)=>{
   const href=$(element).attr('href');
   assert.ok(href.startsWith('#'));
   assert.equal($(`[id="${href.slice(1)}"]`).length,1,'every glyph resolves inside this page');
  });
  assert.equal($('mjx-container svg:not([aria-hidden="true"])').length,0);
  assert.ok($('mjx-container[role="math"]').attr('aria-label').startsWith(latex));
 }
});
test('inline math flows without a display block or copy toolbar',()=>{
 const $=render(String.raw`x\in\mathbb{R}`,true);
 assert.equal($('.formula-inline mjx-container:not([display])').length,1);
 assert.equal($('.formula-block,.formula-copy').length,0);
});
test('imported equations retain anchors and TeX while AMS regenerates their numbers',()=>{
 const $=load('<div class="equation" id="eq-42"><span class="eq-number">(42)</span><span class="math-display" role="math" aria-label="x &lt; y"><svg></svg></span></div>');
 renderMath($);
 assert.equal($('#eq-42.formula-block').length,1);
 assert.equal($('#eq-42').attr('data-equation-number'),'1.1');
 assert.equal($('#eq-42 svg[data-labels] g[id]').attr('id'),'mjx-eqn:1.1');
 assert.equal($('#eq-42 .eq-number').length,0);
 assert.equal($('#eq-42').attr('data-latex'),'x < y');
 assert.equal($('#eq-42 mjx-container > svg').length,1);
});
test('invalid LaTeX fails the build instead of silently publishing broken formulas',()=>{
 assert.throws(()=>render(String.raw`\unknownCommand{x}`),/Invalid LaTeX/);
});

test('AMS counters restart per article, replace old labels, and exclude inline math',()=>{
 const $=load('<div data-math="a=b"></div><span data-display="inline" data-math="x"></span><div class="equation"><span class="eq-number">1.1</span><span data-math="c=d"></span></div><div data-math="e=f"></div>');
 renderMath($);
 assert.deepEqual($('.formula-block[data-equation-number]').map((_,el)=>$(el).attr('data-equation-number')).get(),['1.1','1.2','1.3']);
 assert.equal(render('x=y')('.formula-block').attr('data-equation-number'),'1.1');
 assert.equal($('.formula-inline .eq-number').length,0);
});
test('explicit AMS tags are not numbered twice',()=>{
 const $=render(String.raw`x=y\tag{A}`);
 assert.equal($('[data-mml-node="mlabeledtr"]').length,1);
 assert.equal($('.eq-number').length,0);
});

test('automatic display numbering follows h2 sections and ignores h3 subheadings',()=>{
 const $=load('<h2>First</h2><div data-math="a=b"></div><h3>Details</h3><div data-math="c=d"></div><h2>Second</h2><span data-display="inline" data-math="x"></span><div data-math="e=f"></div>');
 renderMath($);
 assert.deepEqual($('.formula-block[data-equation-number]').map((_,el)=>$(el).attr('data-equation-number')).get(),['1.1','1.2','2.1']);
});

test('an overview counts as chapter one even without display equations',()=>{
 const $=load('<h2>Overview</h2><h2>First</h2><div data-math="a=b"></div><div data-math="c=d"></div><h2>Second</h2><div data-math="e=f"></div>');
 renderMath($);
 assert.deepEqual($('.formula-block[data-equation-number]').map((_,el)=>$(el).attr('data-equation-number')).get(),['2.1','2.2','3.1']);
});

test('repeated formulas share glyphs without adding duplicate IDs or external font requests',()=>{
 const $=load('<div data-math="x+x"></div><div data-math="x+x"></div>');
 renderMath($);
 const ids=$('[id]').map((_,el)=>$(el).attr('id')).get();
 assert.equal(ids.length,new Set(ids).size);
 assert.equal($('.math-font-cache').length,1);
 assert.equal($('.math-font-cache path[id$="-1D465"]').length,1);
 assert.equal($('mjx-container use[data-c="1D465"]').length,4);
 assert.equal($('mjx-container path').length,0);
});

test('AMS numbers rows and advances the next equation without injecting tag commands',()=>{
 const $=load('<h2>Chapter</h2><div id="rows"></div><div id="next" data-math="x=y"></div>');
 $('#rows').attr('data-math',String.raw`\begin{align}a&=b\\c&=d\end{align}`);
 renderMath($);
 assert.deepEqual(JSON.parse($('#rows').attr('data-equation-numbers')),['1.1','1.2']);
 assert.equal($('#next').attr('data-equation-number'),'1.3');
 assert.equal($('#rows').attr('data-latex').includes(String.raw`\tag`),false);
 assert.equal($('[data-latex]').toArray().some(el=>($(el).attr('data-latex')||'').includes(String.raw`\tag`)),false);
});
