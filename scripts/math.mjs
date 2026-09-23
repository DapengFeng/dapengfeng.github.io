import {mathjax} from '@mathjax/src/js/mathjax.js';
import {TeX} from '@mathjax/src/js/input/tex.js';
import {SVG} from '@mathjax/src/js/output/svg.js';
import {liteAdaptor} from '@mathjax/src/js/adaptors/liteAdaptor.js';
import {RegisterHTMLHandler} from '@mathjax/src/js/handlers/html.js';
import '@mathjax/src/js/util/asyncLoad/esm.js';
import '@mathjax/src/js/input/tex/base/BaseConfiguration.js';
import '@mathjax/src/js/input/tex/ams/AmsConfiguration.js';
import '@mathjax/src/js/input/tex/newcommand/NewcommandConfiguration.js';
import '@mathjax/src/js/input/tex/tagformat/TagFormatConfiguration.js';

const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
const svg = new SVG({fontCache: 'global', useXlink: false});
// Load font paths once so article rendering is synchronous and works offline.
await svg.font.loadDynamicFiles();
let numberingSection = 1;
const math = mathjax.document('', {
  InputJax: new TeX({packages: ['base', 'ams', 'newcommand', 'tagformat'], tags: 'ams',
    tagformat: {number: n => `${numberingSection}.${n}`},
    formatError(_jax, error) { throw new Error(`Invalid LaTeX: ${error.message}`); }}),
  OutputJax: svg
});
export function mathStyles() { return adaptor.cssText(svg.styleSheet(math)); }
function typeset(latex, inline) {
  // AMS owns the counter. Bare display expressions get a numbered equation environment.
  const environment = /^\s*\\begin\{(?:equation|align|alignat|flalign|gather|multline|eqnarray)\*?\}/.test(latex);
  const input = inline || environment ? latex : `\\begin{equation}\n${latex}\n\\end{equation}`;
  const node = math.convert(input, {display: !inline, em: 16, ex: 8, containerWidth: 1280});
  adaptor.setAttribute(node, 'role', 'math');
  adaptor.setAttribute(node, 'aria-label', latex);
  return adaptor.outerHTML(node);
}

// Normalize authored LaTeX and imported MathJax SVGs using their original TeX labels.
export function renderMath($, headings) {
  svg.clearFontCache();
  numberingSection = 1;
  math.inputJax[0].reset();
  // Use the same semantic chapters as the contents, including paired heading wrappers.
  const chapters = headings && new Map(headings.filter(h => h.level === 2).map((h, i) => [h.id, i + 1]));
  $('.math-inline[aria-label],.math-display[aria-label]').each((_, element) => {
    const node = $(element);
    if (!node.attr('data-math')) node.attr('data-math', node.attr('aria-label'));
    node.attr('data-display', node.hasClass('math-inline') ? 'inline' : 'display');
  });
  let chapter = 0;
  $('[id],h2,[data-math]').each((_, element) => {
    const current = $(element);
    if (chapters?.has(current.attr('id'))) chapter = chapters.get(current.attr('id'));
    if (!chapters && current.is('h2') && !current.closest('header,footer,nav,.hero,.rail,.site-header,[data-legacy-chrome]').length) chapter++;
    if (Math.max(1, chapter) !== numberingSection) {
      numberingSection = Math.max(1, chapter);
      math.inputJax[0].reset();
    }
    if (!current.is('[data-math]')) return;
    const node = $(element), latex = node.attr('data-math') || node.text();
    const inline = node.attr('data-display') === 'inline';
    const rendered = typeset(latex, inline);
    node.removeAttr('role aria-label');
    node.attr('data-math', latex).attr('data-display', inline ? 'inline' : 'display');
    if (inline) {
      node.addClass('formula-inline').html(rendered);
      return;
    }
    // Reuse legacy equation containers to preserve anchors and the original expression.
    let block;
    const equation = node.parent('.equation');
    if (equation.length && equation.children('.math-display,[data-math]').length === 1) block = equation;
    else if (node.is('div,figure,section')) block = node;
    else { node.wrap('<div></div>'); block = node.parent(); }
    // Imported HTML numbers are replaced by AMS numbering, never injected as \tag.
    block.children('.eq-number').remove();
    const numbers = $(rendered).find('svg[data-labels] g[id^="mjx-eqn:"]').map((_, el) => $(el).attr('id').slice(8)).get();
    if (numbers.length) block.attr('data-equation-number', numbers[0]).attr('data-equation-numbers', JSON.stringify(numbers));
    const tools = $('<div class="formula-tools"></div>');
    tools.append('<button type="button" class="formula-copy" title="Copy LaTeX / 复制 LaTeX"><svg viewBox="0 0 24 24" width="20" height="20" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><g class="copy-symbol"><rect x="8" y="8" width="12" height="12" rx="2"></rect><path d="M16 8V5a2 2 0 0 0-2-2H5a2 2 0 0 0-2 2v9a2 2 0 0 0 2 2h3"></path></g><path class="copied-symbol" d="m5 12 4 4L19 6"></path></svg><span class="formula-copy-label sr-only"><span data-lang="en" lang="en">Copy LaTeX</span><span data-lang="zh" lang="zh-CN">复制 LaTeX</span></span></button>');
    const scroll = $('<div class="formula-scroll" tabindex="0" role="region" aria-label="LaTeX"></div>');
    // Keep the source node when it is separate from the container.
    if (block[0] === node[0]) scroll.html(rendered);
    else { node.html(rendered); scroll.append(node); }
    const body = $('<div class="formula-body"></div>').append(scroll);
    block.addClass('formula-block').attr('data-latex', latex).empty().append(tools, body);
    block.append('<span class="formula-status sr-only" role="status" aria-live="polite"></span>');
  });
  const containers = $('mjx-container[jax="SVG"]');
  // The named math container is the accessible object; visual glyph SVGs are decorative.
  containers.find('svg').attr('aria-hidden', 'true').attr('focusable', 'false').removeAttr('role');
  // Collapse transform-only wrappers without changing transform order or presentation.
  containers.find('g').toArray().reverse().forEach(element => {
    const group = $(element), attrs = element.attribs;
    if (['mlabeledtr','mtable','merror','maction'].includes(attrs['data-mml-node']) || Object.keys(attrs).some(key => !['data-mml-node','data-latex','data-mjx-texclass','transform'].includes(key))) return;
    const children = group.children();
    if (children.length === 1 && children.is('g,path,use,rect,line')) {
      const transform = [attrs.transform, children.attr('transform')].filter(Boolean).join(' ');
      if (transform) children.attr('transform', transform);
      group.replaceWith(children);
    } else if (!attrs.transform) group.replaceWith(group.contents());
  });
  if (containers.length) {
    const cache = $('<svg class="math-font-cache" xmlns="http://www.w3.org/2000/svg" aria-hidden="true" focusable="false" width="0" height="0" style="position:absolute;overflow:hidden"></svg>');
    cache.html(adaptor.outerHTML(svg.fontCache.getCache()));
    $('body').prepend(cache);
  }

}
