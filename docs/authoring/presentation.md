# Page design and interactive examples / 页面设计与交互示例

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md) · [Writing / 写作](writing.md)

## Page purpose and hierarchy / 页面职责与层级

Start from the reader's intended gain: understanding an idea or noticing a specific detail. Give the material that delivers it the clearest position, enough room and a readable sequence. These design criteria apply to current and future subjects; they do not prescribe one layout or claim that every existing page already meets them.

从读者应获得什么出发：理解一个问题，或注意到一个具体细节。为支撑这一收获的材料提供明确的位置、足够的空间与清楚的阅读顺序。这些设计标准适用于现有与未来主题，不要求统一布局，也不代表每个现有页面都已达到。

| Page / 页面 | Reader's task and design priority / 读者任务与设计重点 |
| --- | --- |
| Homepage / 首页 | Understand the site's purpose and choose a first read. Pair a brief mission with representative work whose subjects show the range.<br>理解网站方向并选择阅读起点。简短使命配合代表作，由作品的具体问题体现内容跨度。 |
| Notebook / 知识库 | Find material relevant to a question through scannable titles, useful summaries and topic groupings.<br>通过便于浏览的标题、有信息量的简介与主题分组，找到与问题相关的内容。 |
| Series / 专题 | Understand the learning goal, prerequisites and order, including what each part contributes.<br>了解学习目标、所需基础、阅读顺序与每一期的作用。 |
| Technical article / 技术文章 | Follow a phenomenon through its mechanism and a checkable example; keep each explanation near its figure, equation or code.<br>从现象理解机制与可检验的例子，让解释靠近对应图示、公式或代码。 |
| Travel essay / 游记 | Follow an account of places, people and moments through prose and photographs.<br>通过文字与照片进入叙述，留意其中的地方、人物与片刻。 |
| About / 关于 | Understand who writes, why and where to begin. Support the mission with factual author information and representative work, leaving room for future subjects.<br>理解谁在写、为什么写、从哪里读起。真实作者信息与代表作支撑使命，为未来主题留出空间。 |

New page types may serve another reader task. Share navigation and control meanings while allowing content structure to differ. Actual recommendation and grouping rules belong in [publishing](publishing.md); this table adds no selection or filtering features.

新增页面可以承担不同的读者任务。共用导航与操作含义保持一致，内容结构允许不同。推荐与分组的实际规则见[文章与学习专题](publishing.md)，此表不引入新的取文或筛选功能。

The first screen should establish the page's subject and where reading or selection begins. Avoid oversized identity blocks or empty introductions that push all substantive content out of view; do not force a whole argument into one screen. Use headings, spacing and alignment to put the subject first, supporting evidence and local navigation next, and metadata or optional actions last. Add cards only for genuinely parallel items; empty recommendation slots need no weak substitutes.

首屏应说明页面主题，以及从哪里开始阅读或选择。避免过大的身份展示或空旷引导区把实质内容全部推到屏外，也不强求一屏讲完整个论证。用标题、间距与对齐，让主题优先、辅助证据与局部导航其次、元信息与可选操作靠后。确实并列的内容才使用卡片，推荐空位不必用牵强内容补满。

Keep text, keyboard focus and controls legible over images or mathematical backgrounds. A subject-carrying figure may be the focal point; the daily-mathematics background stays behind the reading surface. Its principle explanation appears beside the homepage animation at a quieter, slightly smaller typographic level than the main foreground copy, while formulas and prose remain readable.

文字、键盘焦点与控件在图片或数学背景上仍须清晰。承载主题的图示可以成为视觉中心，每日数学背景则处于阅读内容之后。首页原理说明与动画相邻，字号比主要前景文案稍小、层级更安静，公式与文字仍须可读。

## Typography and links / 字号与链接

Use the shared reading styles: technical articles start from 18px body text and 16px code and controls, with space between content blocks. Genre-specific typography belongs in shared styles: travel prose uses 19px on desktop and 18px on narrow screens in `src/styles/editorial.css`. Use `lesson-note` for explanations and `lesson-scroll` for wide tables. Link text should name the destination document, API, or source file and symbol. Underlines, color and hover/focus feedback identify links; do not append decorative arrows to links, navigation or buttons. Retain arrows that express direction or mappings in equations, code and explanatory diagrams.

使用共用阅读样式：技术文章以正文 18px、代码与控件文字 16px 为基础，内容块之间留出空隙。不同体裁的字体差异集中在共用样式中，例如 `src/styles/editorial.css` 将游记正文设为桌面 19px、窄屏 18px。说明使用 `lesson-note`，宽表格使用 `lesson-scroll`。链接文字应写明目标文档、API 或源码文件与符号。下划线、颜色及悬停和焦点反馈标识链接；链接、导航和按钮不附加装饰箭头，公式、代码与知识图示中表示方向或映射的箭头保留。

[The article template](../../examples/article.html) includes bilingual paragraphs, inline and numbered AMS equations, highlighted C++ fragments, an editable complete C++ program, a shared interactive figure, a table, and explanation blocks. Copy its source into `content/posts/` and build to preview the site styling and controls.

[文章模板](../../examples/article.html)包含双语段落、行内与编号 AMS 公式、配色 C++ 片段、可编辑完整 C++ 程序、共用交互图、表格和说明块。将源码复制到 `content/posts/` 后构建，即可预览站点样式和控件。

## Visual cues and compact layout / 视觉提示与紧凑布局

Use recognizable icons, consistent control dimensions, hover/focus states and input styling for familiar actions. Play, step, reset, copy, share and discussion need no repeated instruction banners; editable fields should look editable. Keep localized accessible names, keyboard access, selected/disabled states, necessary labels and units, explanations of unfamiliar actions and actionable errors.

常见操作采用易识别的图标、一致的控件尺寸、悬停与焦点状态和输入框样式。播放、单步、恢复、复制、分享与讨论无需重复操作横幅，可编辑区域应显出可输入的特征。保留双语可访问名称、键盘操作、选中与禁用状态、必要标签与单位、陌生操作的解释，以及可处理的错误说明。

Keep a figure's input, controls, result and short explanation close enough to read together, preferably within one viewport when practical. Adapt composition for both width and height: stack narrow-screen comparisons in a meaningful order with nearby labels and legends; reduce incidental spacing on short screens. Preserve these relationships in English, Chinese and bilingual modes. Compactness removes repetition, not essential text, evidence or reasoning: do not shrink labels, crop the subject, hide required steps or force fixed-height panels. Wide equations and tables can scroll within their blocks.

让图示的输入、控件、结果与简短说明彼此靠近，条件允许时组成一屏内的完整讲解单元。同时适配宽度与高度：窄屏对照按有意义的顺序纵向排列，标注和图例放在近旁；矮屏先减少无关间距。英文、中文与双语模式均应保留这些关系。紧凑用于消除重复，不用于删减关键文字、证据或推理；不缩小标注、裁掉主体、隐藏必要步骤或强制固定高度。宽公式与表格可在各自块内滚动。

Provide the initial, loading, success, empty and failure states a component needs. Slow or unavailable networks must leave the core explanation readable and distinguish an unfinished request from a result. Supply useful media fallbacks and load optional work on demand; preserve the [language](language.md), [discussions](../operations/discussions.md), [sharing](../operations/sharing.md), [support](../operations/support.md) and [analytics](../operations/analytics.md) boundaries when changing those services.

按组件需要提供初始、加载、成功、空缺与失败状态。慢网或断网时核心解释仍可读，未完成的请求与已有结果应区分。媒体提供有信息量的回退内容，可选工作按需加载；修改相关服务时遵守[语言](language.md)、[讨论](../operations/discussions.md)、[分享](../operations/sharing.md)、[赞赏](../operations/support.md)与[统计](../operations/analytics.md)的边界。

Keep reader-facing text about the subject or a necessary action. Build details, publication-date maintenance rules and claims about a model's polish belong in repository documentation or verification notes. Retain scientific assumptions, uncertainty, source attribution and failure explanations that affect understanding.

面向读者的文字聚焦内容与必要操作。构建细节、发布日期维护规则、对模型精美程度的自我描述，放在仓库文档或验证记录中。理解内容所需的科学假设、不确定性、来源和失败原因仍须保留。

## Mathematics / 数学公式

Write SVG, Canvas, tables, and JavaScript directly. The HTML element `<div data-math="y=Ax"></div>` is rendered with MathJax and its AMS extension at build time as offline SVG. Glyphs are shared within each page to reduce repeated paths; each formula has one accessible LaTeX name, with its visual SVG hidden from screen readers. For inline mathematics, use `<span data-math="x" data-display="inline"></span>`. Display equations center the expression within the block and right-align its number, use automatic section-based numbers (1.1, 1.2, 2.1) generated by the AMS counter and rendered inside the formula SVG in parentheses at the right, and include one icon to copy their original LaTeX. For long derivations, use `aligned` to place complete equations or derivation steps on separate rows and align their relation symbols. The renderer does not split an equation automatically; if a complete row is still too wide, its block scrolls horizontally without shrinking the text. Imported equation SVGs with LaTeX labels are normalized to the same style. Use automatic numbering throughout an article; remove old hand-written eq-number labels when migrating. Count logical chapters from 1, including the overview; paired English and Chinese headings count as one chapter. Equation prefixes use the same chapter sequence as the contents.

直接编写 SVG、Canvas、表格与 JavaScript。HTML 元素 `<div data-math="y=Ax"></div>` 会在构建时通过 MathJax 及其 AMS 扩展渲染为 SVG，支持离线显示。页内复用字形以减少重复路径；每个公式保留一个可访问的 LaTeX 名称，其视觉 SVG 对屏幕阅读器隐藏。行内公式使用 `<span data-math="x" data-display="inline"></span>`。独立公式的主体在块内居中、编号右对齐，按章节自动编号（1.1、1.2、2.1），编号由 AMS 自动计数，在公式 SVG 内排版，带圆括号并位于右侧，并提供一个复制原始 LaTeX 的图标。长推导使用 `aligned` 按完整等式或推导步骤分行，并对齐关系符号。渲染器不自动拆开一个等式；完整一行仍然过宽时，在块内横向滚动，不缩小字号。带有 LaTeX 标签的旧公式 SVG 也会统一排版。整篇文章统一自动编号，迁移时移除旧的手写 eq-number 标签。逻辑章节（含概览）从 1 开始计数，成对的中英文标题只算一章；公式前缀直接使用目录的同一章节序列。

Numbering uses MathJax’s `tags: "ams"` counter, reset at each chapter. No generated `\tag` is inserted. Use `aligned` or `gathered` for one number on a multiline block; `align` or `gather` numbers individual rows and advances subsequent equations accordingly.

编号使用 MathJax 的 `tags: "ams"` 计数器，每章重置，不插入生成的 `\tag`。多行共用一个编号时使用 `aligned` 或 `gathered`；逐行编号时使用 `align` 或 `gather`，后续公式序号自动顺延。

The AMS extension supports environments such as `align`, `gather`, `cases`, and `pmatrix`, along with commands such as `\mathbb` and `\operatorname`. Put TeX directly in `data-math`, without dollar delimiters or `\usepackage`. Escape HTML attribute characters, for example `&amp;` for an alignment ampersand.

AMS 扩展支持 `align`、`gather`、`cases`、`pmatrix` 等环境，以及 `\mathbb`、`\operatorname` 等命令。在 `data-math` 中直接写 TeX，不需要美元分隔符或 `\usepackage`。HTML 属性中的特殊字符需要转义，例如对齐用的与号写为 `&amp;`。

## Article-local styles / 文章局部样式

Inline CSS in standalone articles is scoped to `.legacy-content` to avoid overriding site navigation; inline scripts are preserved. Use unique IDs or article-local selectors. Shared styles live in `src/styles/`, and interaction scripts in `src/scripts/`. Publish only HTML and scripts you trust.

独立文章的内联 CSS 会限定在 `.legacy-content` 范围，避免覆盖站点导航；内联脚本会保留。请使用唯一 ID 或文章局部选择器。共享样式位于 `src/styles/`，交互脚本位于 `src/scripts/`。仅发布自己信任的 HTML 和脚本。

## Code blocks / 代码块

All block code (including plain `pre`, highlighted fragments, editable examples, and compiler output) receives a compact toolbar with one copy icon on the right, matching display equations. Named code languages appear on the left; plain text omits that label. Copying preserves the current source and indentation without line numbers. Labels follow the language setting; if clipboard access is blocked, a selected text field allows manual copying. Inline code has no button. Do not add article-specific copy controls.

所有块级代码（包括普通 `pre`、高亮片段、可编辑示例和编译输出）会自动添加紧凑工具栏，右侧显示与独立公式一致的复制图标；已标明的代码语言显示在左侧，纯文本省略该标签。复制保留当前源码与缩进，不包含行号。提示遵循语言选择；剪贴板不可用时提供已选中的文本框供手动复制。行内代码不加按钮，无需在文章中手写复制控件。

Name code languages with `data-language` or a `language-*` class. Known language names such as C++ and Python stay visible. Plain text blocks use a distinct visual style and omit the visible “Text / 文本” label.

使用 `data-language` 或 `language-*` 类名标识代码语言，C++、Python 等语言名称保持可见。纯文本块通过样式区分，不显示“Text / 文本”标签。

## Verify code without leaving the article / 在文章内验证代码

Add `data-godbolt="rust"`, `data-godbolt="c++"`, or `data-godbolt="python"` to a complete example’s `<code>` element inside `<pre>`. The article receives an inline Compiler Explorer panel automatically. Editable examples have a line-number gutter, current-line highlight, and focus border. Pencil, copy, reset, and green play icons share the code toolbar; icon tooltips and accessible names follow the selected language. Rust, C++, and Python compile and run; diagnostics, stdout, stderr, exit status, and available assembly appear separately. Pseudocode should remain unmarked.

在完整示例的 `<pre>` 内，为 `<code>` 添加 `data-godbolt="rust"`、`data-godbolt="c++"` 或 `data-godbolt="python"`，文章会自动出现站内 Compiler Explorer 面板。可编辑示例使用行号栏、当前行高亮与焦点边框；铅笔、复制、恢复和绿色播放图标集中在代码工具栏，图标的悬停提示与可访问名称跟随语言选择。Rust、C++ 与 Python 编译后运行，分别显示编译诊断、标准输出、标准错误、退出状态及可用的汇编。伪代码不要添加此标记。

```html
<pre><code data-godbolt="rust">fn main() {
    println!("Hello");
}</code></pre>
```

Requests go directly to the [Compiler Explorer API](https://github.com/compiler-explorer/compiler-explorer/blob/main/docs/API.md) only after a reader clicks the button. No navigation, API key, or server backend is required. The displayed compiler version and flags identify the check. Compilation alone does not prove runtime correctness; a failed network request is never shown as a passed check. Editing code or changing the example clears the old result and cancels any pending request. Edits are local to the page and are reset on reload or example switching; copy buttons copy the edited code. Optional `data-compiler`, `data-compiler-name`, and `data-compiler-options` attributes override the defaults.

只有读者点击按钮后，才会直接请求 [Compiler Explorer API](https://github.com/compiler-explorer/compiler-explorer/blob/main/docs/API.md)。无需跳转、API 密钥或后端服务；面板显示编译器版本与参数。编译成功不证明运行正确，网络失败也不会显示为验证通过。编辑代码或切换示例会清除旧结果并取消尚未完成的请求。修改仅保留在当前页面，刷新或切换示例时恢复；复制按钮复制当前编辑的代码。可通过 `data-compiler`、`data-compiler-name` 和 `data-compiler-options` 覆盖默认配置。

## Reader-first diagrams / 面向阅读的示意图

Give each figure one question to answer. Make its input, transformation and result identifiable, preserve units and assumptions, and compare conditions on shared scales. Keep names and colors consistent across related figures and levels of detail, including organ, tissue and cell views. Name quantities in a nearby legend and combine color with labels, line styles or shapes. Provide an accessible name and text explanation of essential results; mobile labels must remain legible.

每张图围绕一个问题，让输入、变化过程与结果各有明确位置，保留单位与假设，对照条件使用统一尺度。相邻图示及不同细节层次保持名称与配色一致，包括器官、组织与细胞视图。近旁图例注明各个量，颜色配合文字、线型或形状传递含义。提供可访问名称与关键结果的文字解释，移动端标注仍须清晰。

| Form / 形式 | Explanatory purpose / 解释任务 |
| --- | --- |
| Static figure or photograph / 静态图或照片 | A relationship, comparison, composition or observed detail that readers can inspect at their own pace.<br>可停留查看的关系、对照、构图或观察细节。 |
| Animation / 动画 | Meaningful order or change over time, with stages and results understandable without catching a fleeting frame.<br>有意义的先后顺序或时间变化，阶段与结果不依赖捕捉瞬间画面。 |
| 3D / 三维 | Depth, position or spatial structure is part of the question, and the initial view reveals it.<br>深度、位置或空间结构本身是问题的一部分，且默认视角能呈现关系。 |
| Interaction / 交互 | Changing a condition or view answers a stated question; the informative initial result and nearby explanation already convey the main point.<br>改变条件或视角能回答明确问题，默认结果与近旁说明已传递核心信息。 |

Use aligned 2D diagrams and curves for signal processing and calculations. Prefer a clear 2D figure to an awkward 3D view. Spatial models need accurate structure, an informative initial camera and a useful fallback when WebGL or media loading fails. Review scientific accuracy and visual quality separately: smooth motion or polished rendering cannot validate a model.

信号处理与计算优先采用对齐的二维图和曲线，清楚的二维图优于粗糙的三维视图。空间模型需有准确的结构、能说明问题的默认视角，以及 WebGL 或媒体加载失败时有信息量的回退。科学准确性与视觉质量分别核对，流畅动态或精致渲染不能证明模型正确。

For article demonstrations, show the complete result or comparison before controls and use short, finite animation for causal processes. The shared `reading-demos.js` controller plays once in view, pauses offscreen, yields to reader input and displays the complete diagram for reduced-motion preferences. Code execution and benchmarks remain explicitly triggered. Daily mathematics is a separate exception governed only by its [daily-operation policy](daily-mathematics.md#daily-operation--每日运行); do not apply that exception to article demonstrations or navigation transitions.

正文演示先展示完整结果或对照，再提供控件；因果过程使用短时、有限的动画。共享的 `reading-demos.js` 控制器进入视野后演示一次，离开视野时暂停，读者操作后停止自动接管，减少动态效果模式直接显示完整示意。代码执行与基准测试仍需主动触发。每日数学是单独的例外，仅遵循其[每日运行规则](daily-mathematics.md#daily-operation--每日运行)，不将例外扩展到正文演示或导航过渡。

Process diagrams interpolate positions continuously: transpose versus copying, packed-storage mappings, two-output symmetric updates, GPU address ownership, a Frank–Wolfe line search, a LIF event and wave superposition. Numerical models live in `process-models.js`; rendering lives in `process-demos.js`. Spatial animation pauses with the shared controller. Reset events retain the same physical timestamp; GPU particles represent address mappings, not execution timing.

过程图支持位置的连续插值：转置与复制、紧凑存储映射、对称矩阵的两路更新、GPU 地址归属、Frank–Wolfe 线搜索、LIF 事件与波的叠加。数值模型位于 `process-models.js`，绘图位于 `process-demos.js`。空间动画由共享控制器暂停与恢复；复位事件保持相同的物理时刻，GPU 粒子表示地址映射而非执行时序。

## Travel essays and visual character / 游记与视觉气质

Let photographs and prose set the pace. Keep a photograph beside the passage about its actual place or subject, preserve meaningful framing and use short captions. Give prose room between photographs; avoid turning each paragraph into a technical card or filling the page with tags and borders. Technical and travel articles share navigation, controls and accessibility, while body typography and image rhythm may differ. Use the [writing guide](writing.md) for narrative and factual constraints.

让照片与文字建立节奏。图片靠近讲述其真实地点或主体的段落，保留有意义的构图，图注简短。照片之间为文字留出空间，不把每段游记装进技术卡片，也不让标签与框线铺满页面。技术文章与游记共用导航、控件和可访问性，正文排版与图片节奏可不同。叙事与事实约束见[写作指南](writing.md)。

Technical precision comes from exact alignment, consistent notation, clear code and visible model relationships. Human attention comes from the author's particular question, faithful photographs, thoughtful pacing and honest uncertainty. Let these qualities arise from the subject; neon colors, dense dashboards, slogans and ornamental effects cannot establish them.

技术的精确来自准确对齐、一致符号、清楚代码与可见的模型关系；人的观察来自作者具体的问题、忠实照片、细致节奏与诚实的不确定性。让这些特征从内容中形成，霓虹配色、密集仪表盘、套话与装饰特效不能建立它们。

## Review the rendered result / 验收成品页面

Check concrete reader tasks: can a center–surround filter reveal its center value, surround value and subtraction before a slider is touched? Does the mobile photograph preserve the street detail discussed in its adjacent passage? Can readers find the homepage's purpose and reading entry on a short screen, with visible keyboard focus and a slow network? For visual changes, inspect screenshots and relevant interaction states using [visual acceptance](../development/testing.md#visual-acceptance--视觉验收); layout assertions alone do not establish clear hierarchy or an informative figure.

用具体阅读任务检查：中心—周围滤波图能否在操作滑块之前呈现中心值、周围值与相减过程？手机照片是否保留紧邻文字讨论的街景细节？矮屏、键盘操作与慢网下，读者能否找到首页方向、阅读入口及其焦点？视觉改动按[视觉验收](../development/testing.md#visual-acceptance--视觉验收)查看截图与相关交互状态，布局断言通过不足以证明层级清楚或图示传递了知识。
