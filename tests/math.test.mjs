import test from 'node:test';
import assert from 'node:assert/strict';
import {load} from 'cheerio';
import {renderMath} from '../scripts/math.mjs';

function render(latex, inline=false) {
 const $=load('<div id="formula"></div>');
 $('#formula').attr('data-math',latex).attr('data-display',inline?'inline':'display');
 renderMath($);return $;
}
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
  assert.equal($('mjx-container[role="math"]').attr('aria-label'),latex);
 }
});
test('inline math flows without a display block or copy toolbar',()=>{
 const $=render(String.raw`x\in\mathbb{R}`,true);
 assert.equal($('.formula-inline mjx-container:not([display])').length,1);
 assert.equal($('.formula-block,.formula-copy').length,0);
});
test('imported equations retain their anchor, number, and TeX',()=>{
 const $=load('<div class="equation" id="eq-42"><span class="eq-number">(42)</span><span class="math-display" role="math" aria-label="x &lt; y"><svg></svg></span></div>');
 renderMath($);
 assert.equal($('#eq-42.formula-block').length,1);
 assert.equal($('#eq-42 .eq-number').text(),'(42)');
 assert.equal($('#eq-42').attr('data-latex'),'x < y');
 assert.equal($('#eq-42 mjx-container > svg').length,1);
});
test('invalid LaTeX fails the build instead of silently publishing broken formulas',()=>{
 assert.throws(()=>render(String.raw`\unknownCommand{x}`),/Invalid LaTeX/);
});

test('display numbers restart per article, skip reserved numbers, and exclude inline math',()=>{
 const $=load('<div data-math="a=b"></div><span data-display="inline" data-math="x"></span><div class="equation"><span class="eq-number">1.1</span><span data-math="c=d"></span></div><div data-math="e=f"></div>');
 renderMath($);
 assert.deepEqual($('.eq-number').map((_,el)=>$(el).text()).get(),['1.2','1.1','1.3']);
 assert.equal(render('x=y')('.eq-number').text(),'1.1');
 assert.equal($('.formula-inline .eq-number').length,0);
});
test('explicit AMS tags are not numbered twice',()=>{
 const $=render(String.raw`x=y\tag{A}`);
 assert.equal($('[data-mml-node="mlabeledtr"]').length,1);
 assert.equal($('.formula-tools .eq-number').length,0);
});

test('automatic display numbering follows h2 sections and ignores h3 subheadings',()=>{
 const $=load('<h2>First</h2><div data-math="a=b"></div><h3>Details</h3><div data-math="c=d"></div><h2>Second</h2><span data-display="inline" data-math="x"></span><div data-math="e=f"></div>');
 renderMath($);
 assert.deepEqual($('.eq-number').map((_,el)=>$(el).text()).get(),['1.1','1.2','2.1']);
});

test('an overview counts as chapter one even without display equations',()=>{
 const $=load('<h2>Overview</h2><h2>First</h2><div data-math="a=b"></div><div data-math="c=d"></div><h2>Second</h2><div data-math="e=f"></div>');
 renderMath($);
 assert.deepEqual($('.eq-number').map((_,el)=>$(el).text()).get(),['2.1','2.2','3.1']);
});

test('repeated formulas share glyphs without adding duplicate IDs or external font requests',()=>{
 const $=load('<div data-math="x+x"></div><div data-math="x+x"></div>');
 renderMath($);
 const ids=$('[id]').map((_,el)=>$(el).attr('id')).get();
 assert.equal(ids.length,new Set(ids).size);
 assert.equal($('.math-font-cache').length,1);
 assert.equal($('.math-font-cache path').length,2);
 assert.equal($('mjx-container use').length,6);
 assert.equal($('mjx-container path').length,0);
});
