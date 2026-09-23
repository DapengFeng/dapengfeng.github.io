# FENG / Site Maintenance Guide / 知识实验室 · 网站维护说明

A static website built with native HTML, CSS, and JavaScript, without Jekyll, React, Vue, or a client-side routing framework. Each article is one HTML file containing its English and Chinese text, metadata, and any article-specific styles or scripts. Node.js scripts generate navigation, tables of contents, categories, dates, full-text search, RSS, and the sitemap at build time. The output can be hosted directly on GitHub Pages.

这是使用原生 HTML、CSS 和 JavaScript 构建的静态网站，不使用 Jekyll、React、Vue 或客户端路由框架。每篇文章使用一个 HTML 文件，包含英文、中文、元信息以及文章自己的样式和脚本。Node.js 脚本仅在构建时生成导航、目录、分类、日期、全文搜索、RSS 与站点地图，产物可直接托管到 GitHub Pages。

## Run locally / 本地运行

```bash
npm ci
npm run dev
```

Requires Node.js 24. The default preview address is `http://localhost:4173`. Changes to HTML, CSS, JavaScript, or article files trigger a rebuild; refresh the browser to see the update.

需要 Node.js 24。默认预览地址为 `http://localhost:4173`。修改 HTML、CSS、JavaScript 或文章后会自动重新构建，刷新浏览器即可查看更新。

```bash
npm run build
npm run preview
npm test
npm run test:browsers
```

Browser checks use Playwright. Install Chromium with `npx playwright install chromium`, or set `PLAYWRIGHT_CHROMIUM_EXECUTABLE` to an existing Chromium executable. Run `npm run test:browsers` for interaction, language, sorting, layout, formula, and accessibility checks. The compiler tests use isolated API responses; add `--live` to `node tests/compiler.mjs` to verify Godbolt itself. Static output is written to `dist/` and requires no server-side application.

浏览器检查使用 Playwright。可执行 `npx playwright install chromium` 安装 Chromium，或通过 `PLAYWRIGHT_CHROMIUM_EXECUTABLE` 指定已安装的 Chromium。`npm run test:browsers` 检查交互、语言、排序、排版、公式和无障碍。编译器测试使用隔离的接口响应；运行 `node tests/compiler.mjs --live` 可验证 Godbolt 本身。静态产物输出到 `dist/`，无需服务器端程序。

## Publish an HTML article / 发布 HTML 文章

1. Copy `examples/article.html` to `content/posts/your-slug.html`. This is the single article directory, including imported visual essays. Subdirectories are supported; article filenames must be unique across the site. Published URLs remain `/blog/your-slug.html`.

   把 `examples/article.html` 复制到 `content/posts/your-slug.html`。所有文章（包括导入的交互长文）统一放在这里。支持子目录，文章文件名必须全站唯一。发布链接仍为 `/blog/your-slug.html`。

2. Set the title, summary, category, tags, and date. Use `YYYY-MM-DD` for `date`, recording the first push or sharing date. Keep it unchanged on later edits; add `updated` to record an update date.

   填写标题、摘要、分类、标签与日期。`date` 使用 `YYYY-MM-DD`，记录首次推送或分享日期；后续编辑保留该日期，如需注明更新则添加 `updated`。

3. Write both languages directly in the same HTML: English first with `data-lang="en"`, immediately followed by Chinese with `data-lang="zh"`. Use `class="parallel-text"` on paired headings or paragraphs, as in the example below. Equations, code, SVG, Canvas, and experiments can sit outside language markers to be shared.

   直接在同一个 HTML 中写两种语言：英文使用 `data-lang="en"`，紧接的中文使用 `data-lang="zh"`。配对的标题或段落使用 `class="parallel-text"`，具体写法见下方示例。公式、代码、SVG、Canvas 与实验可放在语言标记之外，由两种语言共用。

4. Build locally or push to the main branch. New articles automatically appear in the notebook, categories, timeline, search, RSS, and sitemap. Headings `h2` and `h3` receive anchors automatically; `h2` headings also populate the reading contents.

   在本地构建或推送到主分支。新文章会自动进入知识库、分类、时间线、搜索、RSS 与站点地图。`h2` 和 `h3` 自动获得锚点，`h2` 同时进入阅读目录。

Place complete metadata in `<script type="application/json" id="article-metadata">`. The builder also accepts `title`, `meta[name=description]`, `meta[name=category]`, and `meta[name=article:published_time]`. All metadata lives in the article itself; no separate catalog or translation files are needed. Missing required metadata, invalid dates, or a missing language stop the build.

完整元信息写在 `<script type="application/json" id="article-metadata">` 中。构建器也支持 `title`、`meta[name=description]`、`meta[name=category]` 和 `meta[name=article:published_time]`。所有元信息都放在文章自身，无需单独的目录配置或译文文件。必填元信息缺失、日期无效或缺少一种语言都会阻止构建。

```html
<h2 class="parallel-text">
  <span lang="en" data-lang="en">The idea</span>
  <span lang="zh-CN" data-lang="zh">核心思路</span>
</h2>
<p class="parallel-text">
  <span lang="en" data-lang="en">English explanation.</span>
  <span lang="zh-CN" data-lang="zh">对应的中文解释。</span>
</p>
<div data-math="y = Ax"></div>
```

| Category ID<br>分类 ID | Category<br>分类 |
| --- | --- |
| `math` | Mathematics & Algorithms<br>数学与算法 |
| `physics` | Physics & Models<br>物理与模型 |
| `systems` | Systems & Performance<br>系统与性能 |
| `benchmark` | Benchmarks<br>基准测试 |

`draft: true` hides a draft. `archiveOnly: true` includes an article in the archive and search, but excludes it from the homepage notebook. A numeric `featured` sets its order among featured notes. Choose `spike`, `rust`, `gpu`, `matrix`, `math`, `wave`, or `benchmark` for `art` to generate local SVG illustrations.

`draft: true` 隐藏草稿。`archiveOnly: true` 让文章仅进入归档与搜索，不进入首页知识库。数字形式的 `featured` 指定首页精选顺序。`art` 可选 `spike`、`rust`、`gpu`、`matrix`、`math`、`wave` 或 `benchmark`，用于生成本地 SVG 配图。

## Equations, diagrams, and interaction / 公式、图表与交互

To add an article to a learning series, include an optional `series` object in its article metadata. Use the same `id`, `titleEn`, and `title` across installments and a unique positive `part` number. The build generates `/series/<id>/`, links it from the homepage and notebook, adds it to the sitemap, and connects published installments with previous/next navigation. Drafts are excluded; planned articles should remain plain text until published. Duplicate part numbers and inconsistent series titles fail the build.

若要把文章加入学习专题，在元信息中填写可选的 `series` 对象。各期使用相同的 `id`、`titleEn` 和 `title`，并填写不重复的正整数 `part`。构建自动生成 `/series/<id>/`，在首页和知识库添加入口、写入站点地图，并为已发布文章生成前后期导航。草稿不计入；未发布规划应保留为普通文字。期数重复或专题名称不一致会阻止构建。

```json
"series": {
  "id": "pytorch-internals",
  "titleEn": "Inside PyTorch",
  "title": "PyTorch 源码之旅",
  "part": 2
}
```

`npm run test:pytorch` checks the first installment's C++ model, both animations, language modes, mobile layout, and no-JavaScript reading. Its Python/PyTorch examples require a separate local PyTorch 2.10.0 CPU environment; the website build does not install PyTorch. CPU outputs were checked against that version; CUDA diagrams are source-based illustrations rather than GPU measurements.

`npm run test:pytorch` 检查第一期的 C++ 模型、两组动画、语言模式、手机排版与无 JavaScript 阅读。文中的 Python／PyTorch 示例需要独立的本地 PyTorch 2.10.0 CPU 环境，网站构建不会安装 PyTorch。CPU 输出已按该版本核对；CUDA 图示依据源码，不是 GPU 测量结果。

[The article template](examples/article.html) includes bilingual paragraphs, inline and numbered AMS equations, highlighted C++ fragments, an editable complete C++ program, a shared interactive figure, a table, and explanation blocks. Copy its source into content/posts/ and build to preview the site styling and controls.

[文章模板](examples/article.html)包含双语段落、行内与编号 AMS 公式、配色 C++ 片段、可编辑完整 C++ 程序、共用交互图、表格和说明块。将源码复制到 content/posts/ 后构建，即可预览站点样式和控件。

Write SVG, Canvas, tables, and JavaScript directly. The HTML element `<div data-math="y=Ax"></div>` is rendered with MathJax and its AMS extension at build time as offline SVG. Glyphs are shared within each page to reduce repeated paths; each formula has one accessible LaTeX name, with its visual SVG hidden from screen readers. For inline mathematics, use `<span data-math="x" data-display="inline"></span>`. Display equations share one style, use automatic section-based numbers (1.1, 1.2, 2.1) generated by the AMS counter and rendered inside the formula SVG in parentheses at the right, and include one icon to copy their original LaTeX. Imported equation SVGs with LaTeX labels are normalized to the same style. Use automatic numbering throughout an article; remove old hand-written eq-number labels when migrating. Count logical chapters from 1, including the overview; paired English and Chinese headings count as one chapter. Equation prefixes use the same chapter sequence as the contents.

直接编写 SVG、Canvas、表格与 JavaScript。HTML 元素 `<div data-math="y=Ax"></div>` 会在构建时通过 MathJax 及其 AMS 扩展渲染为 SVG，支持离线显示。页内复用字形以减少重复路径；每个公式保留一个可访问的 LaTeX 名称，其视觉 SVG 对屏幕阅读器隐藏。行内公式使用 `<span data-math="x" data-display="inline"></span>`。独立公式统一样式，按章节自动编号（1.1、1.2、2.1），编号由 AMS 自动计数，在公式 SVG 内排版，带圆括号并位于右侧，并提供一个复制原始 LaTeX 的图标。带有 LaTeX 标签的旧公式 SVG 也会统一排版。整篇文章统一自动编号，迁移时移除旧的手写 eq-number 标签。逻辑章节（含概览）从 1 开始计数，成对的中英文标题只算一章；公式前缀直接使用目录的同一章节序列。

Numbering uses MathJax’s `tags: "ams"` counter, reset at each chapter. No generated `\tag` is inserted. Use `aligned` or `gathered` for one number on a multiline block; `align` or `gather` numbers individual rows and advances subsequent equations accordingly.

编号使用 MathJax 的 `tags: "ams"` 计数器，每章重置，不插入生成的 `\tag`。多行共用一个编号时使用 `aligned` 或 `gathered`；逐行编号时使用 `align` 或 `gather`，后续公式序号自动顺延。

The AMS extension supports environments such as `align`, `gather`, `cases`, and `pmatrix`, along with commands such as `\mathbb` and `\operatorname`. Put TeX directly in `data-math`, without dollar delimiters or `\usepackage`. Escape HTML attribute characters, for example `&amp;` for an alignment ampersand.

AMS 扩展支持 `align`、`gather`、`cases`、`pmatrix` 等环境，以及 `\mathbb`、`\operatorname` 等命令。在 `data-math` 中直接写 TeX，不需要美元分隔符或 `\usepackage`。HTML 属性中的特殊字符需要转义，例如对齐用的与号写为 `&amp;`。

Inline CSS in standalone articles is scoped to `.legacy-content` to avoid overriding site navigation; inline scripts are preserved. Use unique IDs or article-local selectors. Shared styles live in `src/styles/`, and interaction scripts in `src/scripts/`. Publish only HTML and scripts you trust.

独立文章的内联 CSS 会限定在 `.legacy-content` 范围，避免覆盖站点导航；内联脚本会保留。请使用唯一 ID 或文章局部选择器。共享样式位于 `src/styles/`，交互脚本位于 `src/scripts/`。仅发布自己信任的 HTML 和脚本。

All block code (including plain `pre`, highlighted fragments, editable examples, and compiler output) receives a compact toolbar with a language label on the left and one copy icon on the right, matching display equations. Copying preserves the current source and indentation without line numbers. Labels follow the language setting; if clipboard access is blocked, a selected text field allows manual copying. Inline code has no button. Do not add article-specific copy controls.

所有块级代码（包括普通 `pre`、高亮片段、可编辑示例和编译输出）会自动添加紧凑工具栏：左侧标识代码语言，右侧显示与独立公式一致的复制图标。复制保留当前源码与缩进，不包含行号。提示遵循语言选择；剪贴板不可用时提供已选中的文本框供手动复制。行内代码不加按钮，无需在文章中手写复制控件。

## Verify code without leaving the article / 在文章内验证代码

Add `data-godbolt="rust"`, `data-godbolt="c++"`, or `data-godbolt="python"` to a complete example’s `<code>` element inside `<pre>`. The article receives an inline Compiler Explorer panel automatically. The code is editable, with a reset button for the current example. Rust, C++, and Python compile and run; diagnostics, stdout, stderr, exit status, and available assembly appear separately. Pseudocode should remain unmarked.

在完整示例的 `<pre>` 内，为 `<code>` 添加 `data-godbolt="rust"`、`data-godbolt="c++"` 或 `data-godbolt="python"`，文章会自动出现站内 Compiler Explorer 面板。代码可直接编辑，并可恢复当前示例。Rust、C++ 与 Python 编译后运行，分别显示编译诊断、标准输出、标准错误、退出状态及可用的汇编。伪代码不要添加此标记。

```html
<pre><code data-godbolt="rust">fn main() {
    println!("Hello");
}</code></pre>
```

Requests go directly to the [Compiler Explorer API](https://github.com/compiler-explorer/compiler-explorer/blob/main/docs/API.md) only after a reader clicks the button. No navigation, API key, or server backend is required. The displayed compiler version and flags identify the check. Compilation alone does not prove runtime correctness; a failed network request is never shown as a passed check. Editing code or changing the example clears the old result and cancels any pending request. Edits are local to the page and are reset on reload or example switching; copy buttons copy the edited code. Optional `data-compiler`, `data-compiler-name`, and `data-compiler-options` attributes override the defaults.

只有读者点击按钮后，才会直接请求 [Compiler Explorer API](https://github.com/compiler-explorer/compiler-explorer/blob/main/docs/API.md)。无需跳转、API 密钥或后端服务；面板显示编译器版本与参数。编译成功不证明运行正确，网络失败也不会显示为验证通过。编辑代码或切换示例会清除旧结果并取消尚未完成的请求。修改仅保留在当前页面，刷新或切换示例时恢复；复制按钮复制当前编辑的代码。可通过 `data-compiler`、`data-compiler-name` 和 `data-compiler-options` 覆盖默认配置。

## Bilingual reading / 双语阅读

On a first visit, IP countries/regions CN, HK, MO, and TW default to Chinese; other locations default to English. Readers can choose English, Chinese, or both, and their saved choice always takes priority. In bilingual mode, English is followed immediately by Chinese.

首次访问时，IP 所在国家或地区为 CN、HK、MO、TW 时默认中文，其余默认英文。读者可选择英文、中文或双语，已保存的手动选择始终优先。双语模式中英文在前，对应中文紧随其后。

Automatic language detection calls [Country](https://country.is/), which receives the visitor’s network IP. The site stores only the selected language, not the IP or country response. Automatic results are cached for the tab session; a saved manual preference skips the lookup. Browser language is used immediately and remains the fallback if the request fails or exceeds 2.5 seconds. VPNs may affect the country result. This README always displays English followed by Chinese.

自动语言判断会请求 [Country](https://country.is/)，服务方会收到访客的网络 IP。本站只保存选中的语言，不保存 IP 或地区查询响应。自动结果在标签页会话中缓存；已有手动偏好时跳过查询。页面先按浏览器语言显示，查询失败或超过 2.5 秒则继续使用该语言。VPN 可能影响地区结果。本 README 始终按英文在前、中文紧随其后的顺序显示。

Interface text comes from `scripts/i18n.mjs`; article titles, summaries, and both body languages live in each article HTML. The build preserves the order you write, derives the bilingual contents from headings, and renders formulas. For shared interactive experiments, keep any language dictionary inside the same HTML, as Spike Notes does.

界面文字由 `scripts/i18n.mjs` 提供；文章标题、摘要与双语正文都在各自的 HTML 中。构建保留你编写的顺序，从标题提取双语目录并渲染公式。共享交互实验所需的语言词典也放在同一 HTML 内，脉冲长文已采用这种方式。

Maintain English and Chinese together in the same file when editing. The build does not translate or assess translation accuracy. Existing articles pair their full explanations in English and Chinese; formulas and code are shared when appropriate. New articles must include both language markers; copy the bilingual template to begin.

修改时在同一文件中同步维护英文和中文。构建不会翻译或判断译文准确性。现有文章的英文与中文完整说明成对排列，公式与代码按需共用。新文章必须包含两种语言标记，可直接复制双语模板开始编写。

## Search and AI discoverability / 搜索与 AI 可发现性

Every build generates unique English-first bilingual titles and descriptions, canonical URLs, Open Graph and Twitter cards, and 1200 × 630 PNG sharing images. Article images use the English title as graphical text. Author, publication and revision dates, chapter anchors, and explicitly cited references are connected through Schema.org JSON-LD. The same static HTML contains both languages and the full article before JavaScript runs. The 404 page is marked `noindex`; drafts and redirects stay out of the sitemap. `npm test` checks these properties and runs automatically before deployment.

每次构建自动生成各页独有、英文在前的双语标题与摘要、规范链接、Open Graph 与 Twitter 分享元信息，以及 1200 × 630 PNG 分享图；文章分享图以英文标题作为图中文字。作者、发布与更新日期、章节锚点和明确标注的参考链接通过 Schema.org JSON-LD 关联。同一个静态 HTML 在 JavaScript 运行前就包含双语完整正文。404 页面标记为 `noindex`，草稿与跳转页不进入站点地图。`npm test` 检查这些规则，并在部署前自动执行。

Keep `titleEn`, `title`, `descriptionEn`, and `description` accurate in each article's metadata. Explain the result, assumptions, methods, and limits in visible bilingual prose, including a text explanation of interactive figures. Link primary sources where claims are made. Add `data-citation` to a source link to include it in structured citations; links inside `#references` containers and `.sources-grid` are also recognized. Only existing links are included. Do not invent sources, measurements, credentials, or update dates. Language buttons are reading modes on one URL; they are not separate `hreflang` editions.

在文章元信息中准确填写 `titleEn`、`title`、`descriptionEn` 与 `description`。用可见的双语正文说明结论、假设、方法和局限，为交互图补充文字解释。在相关论述处链接原始来源；来源链接添加 `data-citation` 后会进入结构化引用，`#references` 容器与 `.sources-grid` 内的链接也会被识别。只收录实际存在的链接，不编造来源、测量结果、资历或更新日期。语言按钮是同一 URL 的阅读模式，不是独立的 `hreflang` 语言版本。

One-time setup: add `https://dapengfeng.github.io/` as a URL-prefix property in Google Search Console and as a site in Bing Webmaster Tools. Choose HTML meta-tag verification and copy only the tag's `content` value into the repository's Actions variables `GOOGLE_SITE_VERIFICATION` and `BING_SITE_VERIFICATION` (Settings → Secrets and variables → Actions → Variables). Deploy, finish verification in each platform, then submit `https://dapengfeng.github.io/sitemap.xml`. Locally, the same environment variables populate the verification tags. Empty variables add no tags. Tokens are public verification identifiers, not account passwords. Account verification and sitemap submission have not been performed by this build.

首次配置：在 Google Search Console 添加 `https://dapengfeng.github.io/` 网址前缀资源，在 Bing Webmaster Tools 添加同一站点。选择 HTML 元标签验证，只把标签的 `content` 值分别填入仓库的 Actions 变量 `GOOGLE_SITE_VERIFICATION` 与 `BING_SITE_VERIFICATION`（Settings → Secrets and variables → Actions → Variables）。部署后在平台完成验证，再提交 `https://dapengfeng.github.io/sitemap.xml`。本地构建也可使用同名环境变量；变量为空时不输出验证标签。这些值是公开验证标识，不是账号密码。构建不会代替账号验证和站点地图提交。

These improvements support crawling, understanding, and attribution; they do not guarantee indexing, rankings, or AI citations. Google says no special GEO schema or `llms.txt` file is required. Monitor actual search queries and indexed pages in the webmaster platforms. See [Google's AI search guidance](https://developers.google.com/search/docs/fundamentals/ai-optimization-guide).

这些改进帮助抓取、理解和来源归属，但不保证收录、排名或 AI 引用。Google 说明无需专用 GEO 标记或 `llms.txt` 文件。请通过站长平台观察真实搜索词与页面收录情况。参见 [Google 的 AI 搜索指南](https://developers.google.com/search/docs/fundamentals/ai-optimization-guide)。

## Deploy to GitHub Pages / 部署到 GitHub Pages

In the repository, go to Settings → Pages → Build and deployment → Source and select **GitHub Actions**. Pushes to `main` or `master` then trigger `.github/workflows/deploy.yml` to build, check, and deploy `dist/`. The workflow installs Chromium and requires all browser checks to pass before publishing. Accessibility failure details are saved as a build artifact. Manual runs are available through `workflow_dispatch`. Local edits do not push or publish themselves.

在仓库 Settings → Pages → Build and deployment → Source 中选择 **GitHub Actions**。随后向 `main` 或 `master` 推送，会触发 `.github/workflows/deploy.yml` 自动构建、检查并部署 `dist/`。工作流会安装 Chromium，浏览器检查全部通过后才发布；无障碍检查失败的明细会保存在构建附件中。也支持通过 `workflow_dispatch` 手动运行。本地修改不会自行推送或发布。

The site uses the root path `/` and targets `https://dapengfeng.github.io`. A project site under a subpath requires a consistent path prefix. Historical date-based article URLs retain redirects. The three newly added HTML articles use the first-push date `2026-09-23`; historical notes retain their recorded dates.

站点使用根路径 `/`，目标地址为 `https://dapengfeng.github.io`。若改成项目子路径站点，需要统一配置路径前缀。历史日期式文章链接保留跳转。三篇新加入的 HTML 文章采用首次推送日期 `2026-09-23`，历史笔记保留原有记录日期。

## Where the implementation lives / 实现位置

- `scripts/templates.mjs`

  Shared HTML page templates.

  共享 HTML 页面模板。

- `scripts/content.mjs`

  HTML discovery, metadata validation, tables of contents, and style scoping.

  HTML 扫描、元信息校验、目录生成与样式隔离。

- `scripts/build.mjs`

  Static pages, RSS, the search index, and the sitemap.

  静态页面、RSS、搜索索引与站点地图。

- `src/scripts/surface.js`

  The animated 3D surface drawn with native Canvas.

  使用原生 Canvas 绘制的三维动态曲面。

- `src/scripts/benchmark-worker.js`

  Measurements on the reader’s device, retaining all raw samples.

  在读者设备上实际测量，保留全部原始样本。
