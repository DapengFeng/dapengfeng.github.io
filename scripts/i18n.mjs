import * as cheerio from 'cheerio';
import { escape as e } from './config.mjs';
export const dictionary = {
 '让知识，变得可见':'Ideas, made visible','知识分类':'Topics','关于实验室':'About the lab','按首次分享时间，回看每一次探索。后续更新保留首次发布日期。':'Revisit each exploration by its first publication date. Later edits keep the original date.',
 '知识实验室':'Knowledge Lab','探索':'Explore','知识库':'Notebook','分类':'Topics','时间线':'Timeline','关于':'About',
 '搜索知识':'Search','跳到主要内容':'Skip to content','FENG 知识实验室首页':'FENG Knowledge Lab home','展开导航':'Open navigation','GitHub（新窗口）':'GitHub (new window)',
 '在公式与直觉之间，搭一座桥。':'A bridge between equations and intuition.', '探索数学之美、物理之理，与代码的性能边界。':'Exploring mathematics, physics, and the limits of code.',
 '开始探索':'Start exploring','浏览全部笔记':'Browse the notebook','篇知识笔记':'notes & ideas','个探索方向':'fields to explore','持续生长中':'Always growing',
 '暂停曲面动画':'Pause surface animation','交互三维阻尼波曲面，调节频率可改变波峰的密度':'Interactive damped wave surface; adjust the frequency to change peak density',
 '频率':'Frequency','拖动，改变你的视角':'Drag to change your perspective','数学 × 物理 × 代码':'MATH × PHYSICS × CODE','选择一个方向':'Follow your curiosity','沿着好奇心，继续深入':'Find a new direction to explore',
 '数学与算法':'Mathematics','物理与模型':'Physics & Models','系统与性能':'Systems & Performance','从抽象公式，到直觉与图形':'From abstract ideas to visual intuition',
 '用模型，解释世界如何运转':'Models that make sense of the world','走近底层，让每个周期有价值':'Closer to the metal. Every cycle counts.',
 '基准测试':'Benchmark',
 '少一点猜测，多一点测量':'Less guessing. More measuring.','值得展开的想法':'Ideas worth unfolding','全部文章':'All notes','交互长文':'Interactive essay',
 '知识索引':'The notebook','每一篇，都是理解世界的一小步。':'Small steps toward a deeper understanding.', '全部':'All','按知识分类筛选':'Filter by topic','文章排序':'Sort notes',
 '最新发布 ↓':'Newest first ↓','最早发布 ↑':'Oldest first ↑','标题 A–Z':'Title A–Z','这个方向还没有匹配的笔记':'No notes match this topic','查看全部文章 →':'View all notes →',
 '按时间浏览完整归档 →':'Explore the complete timeline →','理解，是一场没有终点的探索。':'Understanding is an open-ended adventure.',
 '这里存放推导、实验，以及那些值得反复琢磨的问题。':'Derivations, experiments, and questions worth coming back to.','关于这个实验室 ↗':'About the lab ↗',
 '保持好奇，把复杂的事想明白。':'Stay curious. Make sense of the complex.','文章归档':'Archive','搜索知识库':'Search the notebook','关闭搜索':'Close search',
 '搜索概念、文章、标签…':'Search concepts, notes, and tags…','搜索概念、文章、标签':'Search concepts, notes, and tags',
 '试试「矩阵」「Rust」或「脉冲」':'Try “matrix”, “Rust”, or “spikes”','全文检索 ·':'Full-text search ·','ESC 关闭':'ESC to close',
 '从一个问题开始，在公式、图形与实验之间找到答案。':'Begin with a question. Explore equations, diagrams, and experiments.',
 '筛选标题、摘要与标签…':'Filter titles, summaries, and tags…','筛选知识库':'Filter the notebook','没有找到匹配的笔记':'No matching notes found',
 '试试其他关键词，或换一个探索方向。':'Try another keyword or explore a different topic.','清除筛选 →':'Clear filters →','完整归档 →':'Complete archive →',
 '四个方向，无限联结':'Four fields. Infinite connections','建立索引，也发现知识之间的联系。':'Map the ideas. Discover the connections.','浏览分类 ↗':'Explore topic ↗',
 '想法的时间线':'A timeline of ideas','按分享时间，回看每一次探索。原始发布日期不详的文章，明确标注收录日期。':'Revisit each exploration, ordered by its first publication date.',
 '发布':'Published','收录':'Added','更新':'Updated','早期记录':'Early notes', '你好，我是冯大鹏':'Hi, I’m Dapeng Feng','这是我的个人知识实验室。':'Welcome to my personal knowledge lab.',
 '把抽象的概念，变成可以观察的东西。':'Turn abstract concepts into things you can observe.',
 '我用这个网站记录数学、物理和编程中的思考。一次推导、一张图、一个可调节参数的实验，都可以成为理解的起点。':'I use this site to explore mathematics, physics, and programming. A derivation, a diagram, or an experiment with a parameter to change can be a starting point for understanding.',
 '这里也关注代码的性能：从存储布局到 GPU 计算，从一个直觉到一次可复现的测量。既记录结论，也保留假设、条件和边界。':'Code performance is another focus: from memory layout to GPU computing, from intuition to reproducible measurement. These notes record assumptions and boundaries alongside conclusions.',
 '这个实验室里有什么？':'What’s inside the lab?',
 '数学与算法的图形解释，物理与动力系统的交互模型，系统编程的经验，以及基准测试的设计与测量方法。':'Visual explanations of mathematics and algorithms, interactive physical and dynamical models, systems programming notes, and experiments in benchmarking.',
 '在 GitHub 找到我 ↗':'Find me on GitHub ↗','通过 RSS 订阅 →':'Subscribe via RSS →','知识笔记':'Field notes','本篇目录':'In this note','文章目录':'Table of contents',
 '短篇笔记':'Short note','回到顶部 ↑':'Back to top ↑','打印 / 保存 PDF ↗':'Print / Save PDF ↗','继续探索':'Keep exploring','返回知识库 →':'Back to notebook →',
 '凸优化':'Convex optimization','矩阵计算':'Matrix computation','性能优化':'Performance','动力系统':'Dynamical systems','交互实验':'Interactive lab','内存安全':'Memory safety','线性代数':'Linear algebra','存储布局':'Memory layout','早期项目':'Early project','计算机视觉':'Computer vision','波动':'Waves','相位':'Phase',
 '原文未标注发布日期；以上为本站收录日期。':'The original did not specify a publication date; the date above records its addition to this site.',
 '未找到页面':'Page not found','这里还没有留下笔记。':'There isn’t a note here yet.','返回首页 →':'Back to the lab →',
};
export function pair(en, zh, inline=false) { return `<span class="i18n${inline?' inline':''}"><span data-lang="en">${e(en)}</span><span data-lang="zh">${e(zh)}</span></span>`; }
export function localizePage(html, posts) {
 const $=cheerio.load(html), map={...dictionary};
 for(const p of posts){map[p.title]=p.titleEn||p.title;map[p.description]=p.descriptionEn||p.description;}
 function english(text) {
  const t=text.trim();if(map[t])return map[t];
  if(/^共 \d+ 篇笔记$/.test(t))return `${t.match(/\d+/)[0]} notes in the collection`;
  if(/^\d+ 分钟阅读$/.test(t))return `${t.match(/\d+/)[0]} min read`;
  if(/^\d+ 篇$/.test(t))return `${t.match(/\d+/)[0]} notes`;
  if(t.startsWith('# '))return '# '+(map[t.slice(2)]||t.slice(2));
  if(t.startsWith('/ '))return '/ '+(map[t.slice(2)]||t.slice(2));
  if(t.includes(' / ')&&/[\u3400-\u9fff]/.test(t))return t.split(' / ').map(x=>map[x]||x).join(' / ');
  return null;
 }
 // Replace the home title as a single semantic bilingual unit.
 $('.hero-copy h1').html(`<span class="i18n"><span data-lang="en">Ideas, made<br><em>visible.</em></span><span data-lang="zh">让知识，<em>变得可见。</em></span></span>`);
 const visit = el => {
  for(const node of [...(el.children||[])]){
   if(node.type==='tag' && ['script','style','svg','pre','code'].includes(node.name))continue;
   if(node.type==='script'||node.type==='style')continue;
   if(node.type==='tag' && ($(node).is('[data-lang],.article-body,.i18n') || $(node).is('.prose')&&!$(node).closest('.about-grid').length))continue;
   if(node.type==='text'&&/[\u3400-\u9fff]/.test(node.data)){
    const raw=node.data,t=raw.trim(),en=english(t);
    if(en){
      if($(node.parent).is('title,option')) {$(node.parent).attr('data-en',en).attr('data-zh',t);node.data=`${en} / ${t}`;}
      else $(node).replaceWith((raw.startsWith(' ')?' ':'')+pair(en,t)+(raw.endsWith(' ')?' ':''));
    }
   }else if(node.children)visit(node);
  }
 };
 visit($('body')[0]);
 // English-only supporting labels are decorative, while paired content carries both meanings.
 $('.overline,.category-en,.hero-bottom>span:first-child,.section-heading h2>span:not(.i18n):not(.small-cross):not(.category-symbol),.visual-label,.axis-caption>span').not('[data-lang]').attr('data-decorative-en','true');
 $('input[placeholder],[aria-label],button[title]').each((_,el)=>{
  if($(el).closest('.article-body').length)return;
  for(const attr of ['placeholder','aria-label','title']){const value=$(el).attr(attr);if(!value)continue;const en=english(value);if(en)$(el).attr(`data-${attr}-en`,en).attr(`data-${attr}-zh`,value).attr(attr,`${en} / ${value}`);}
 });
 $('.header-row .github-link').before('<div id="site-language" class="language-switch" role="group" aria-label="Reading language / 阅读语言"><button type="button" data-language-choice="zh" aria-pressed="false" aria-label="Chinese only / 仅中文">中</button><button type="button" data-language-choice="en" aria-pressed="false" aria-label="English only / 仅英文">EN</button><button type="button" data-language-choice="both" aria-pressed="true" aria-label="Bilingual / 中英双语">中/EN</button></div>');
 $('head').append(`<script>(()=>{const valid=v=>['en','zh','both'].includes(v);let saved,automatic;try{saved=localStorage.getItem('feng-language')}catch{}try{automatic=sessionStorage.getItem('feng-auto-language')}catch{}const browser=(navigator.languages?.[0]||navigator.language||'en').toLowerCase().startsWith('zh')?'zh':'en';const value=valid(saved)?saved:valid(automatic)?automatic:browser;document.documentElement.dataset.language=value;document.documentElement.dataset.languageSource=valid(saved)?'manual':valid(automatic)?'auto':'browser';document.documentElement.lang=value==='zh'?'zh-CN':'en'})()</script>`);
 const title=$('title').text(), base=title.replace(' · FENG','');
 const en=english(base); if(en)$('title').text(`${en} / ${base} · FENG`);
 return $.html();
}
