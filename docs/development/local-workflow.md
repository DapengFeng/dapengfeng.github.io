# Local workflow and source map / 本地工作流与源码索引

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

A static website built with native HTML, CSS, and JavaScript, without Jekyll, React, Vue, or a client-side routing framework. Each article is one HTML file containing its English and Chinese text, metadata, and any article-specific styles or scripts. Node.js scripts generate navigation, tables of contents, categories, dates, full-text search, RSS, and the sitemap at build time. The output can be hosted directly on GitHub Pages.

这是使用原生 HTML、CSS 和 JavaScript 构建的静态网站，不使用 Jekyll、React、Vue 或客户端路由框架。每篇文章使用一个 HTML 文件，包含英文、中文、元信息以及文章自己的样式和脚本。Node.js 脚本仅在构建时生成导航、目录、分类、日期、全文搜索、RSS 与站点地图，产物可直接托管到 GitHub Pages。

## Run locally / 本地运行

```bash
npm ci
npm run dev
```

Requires Node.js 24. The default preview address is `http://localhost:4173`. Changes to HTML, CSS, JavaScript, or article files trigger a rebuild; refresh the browser to see the update.

需要 Node.js 24。默认预览地址为 `http://localhost:4173`。修改 HTML、CSS、JavaScript 或文章后会自动重新构建，刷新浏览器即可查看更新。

`npm run dev` rebuilds on source changes; `npm run preview` only serves the current `dist/`. When using preview, rebuild before checking an edit. Set `PORT` to reuse the task's existing preview address. Shared behavior belongs in shared templates, styles or scripts, while article-specific content remains in its HTML.

`npm run dev` 监听源码并重新构建，`npm run preview` 只提供当前 `dist/`。使用 preview 时，查看修改前需先构建；可通过 `PORT` 沿用任务中的预览地址。共用行为放在共用模板、样式或脚本里，文章专属内容保留在对应 HTML 中。

```bash
npm run build
npm run preview
npm test
npm run test:browsers
```

See [testing and validation](testing.md) for browser installation, test coverage, and report locations.

浏览器安装、测试范围与报告位置见[测试与验证](testing.md)。

## Where the implementation lives / 实现位置

- `content/editorial.json`, `scripts/editorial-views.mjs`, `src/styles/editorial.css`

  Mission and about copy, series reading guides, curated article connections and shared layout refinements. Recommendation selection and its configuration are documented in [publishing](../authoring/publishing.md).

  使命与关于页文案、专题阅读指引、人工文章联系，以及共用布局的样式调整。推荐取文与配置归[发布指南](../authoring/publishing.md)维护。

- `scripts/templates.mjs`

  Shared HTML page templates.

  共享 HTML 页面模板。

- `scripts/content.mjs`

  HTML discovery, metadata validation, tables of contents, and style scoping.

  HTML 扫描、元信息校验、目录生成与样式隔离。

- `scripts/build.mjs`

  Static pages, RSS, the search index, and the sitemap.

  静态页面、RSS、搜索索引与站点地图。

- `src/scripts/surface.js`, `src/scripts/daily-math*.js`, `src/styles/daily-math.css`

  Mathematical backgrounds and responsive annotations. Date selection, playback, caching and archive behavior are maintained in [Daily mathematics](../authoring/daily-mathematics.md).

  数学背景与响应式原理注释。日期选择、播放、缓存与归档行为统一在[每日数学](../authoring/daily-mathematics.md)维护。

- `src/scripts/site.js`, `src/styles/site.css`

  The sticky navigation retracts while scrolling down and returns while scrolling up or using the keyboard. Expanded menus and search keep it visible. Transforms preserve the document layout, and reduced-motion preferences disable the transition.

  顶部导航向下滚动时收起，向上滚动或使用键盘时恢复；菜单或搜索打开时保持可见。通过平移保留原有文档布局，系统开启“减少动态效果”时取消过渡动画。

- `src/scripts/benchmark-worker.js`

  Measurements on the reader’s device, retaining all raw samples.

  在读者设备上实际测量，保留全部原始样本。

- [scripts/discussions.mjs](../../scripts/discussions.mjs), [src/scripts/discussions.js](../../src/scripts/discussions.js)

  Article discussion templates and the floating giscus integration.

  文章讨论区模板与浮动 giscus 接入。

- [examples/article.html](../../examples/article.html)

  Copyable examples of the shared article markup and controls.

  可复制的文章结构与共用控件示例。

- `dist/`

  Generated output. Edit sources and rebuild; direct edits here are overwritten.

  自动生成的产物。修改源文件后重新构建；直接修改此目录的内容会被覆盖。

## Keep resources proportional to the page / 按页面需要加载资源

`scripts/content.mjs` detects capabilities such as `hasCode`, `hasCompiler` and `hasReadingDemos`; `lab` comes from article metadata. `scripts/templates.mjs` uses them to include syntax highlighting, copying, compilation and demonstration scripts only on relevant articles. When adding a new markup form or renderer, update capability detection and check both a page that needs it and a page that should not load it. Do not restore unconditional site-wide script loading to fix a missing feature in one article.

`scripts/content.mjs` 识别 `hasCode`、`hasCompiler`、`hasReadingDemos` 等功能标记，`lab` 来自文章元信息；`scripts/templates.mjs` 据此只为相关页面引入高亮、复制、编译和演示脚本。新增标记写法或绘图类型时，同步核对识别逻辑，并检查需要该功能和不应加载该功能的两种页面。不为修复一篇文章漏加载，就恢复全站无条件加载。

Keep feature-specific loading policies in [daily mathematics](../authoring/daily-mathematics.md), [discussions](../operations/discussions.md) and [sharing](../operations/sharing.md). Preserve their boundaries when optimizing; measure the affected requests and visible behavior together rather than making all readers download optional resources on arrival.

[每日数学](../authoring/daily-mathematics.md)、[讨论](../operations/discussions.md)与[分享](../operations/sharing.md)分别维护其专属加载规则。优化时保留这些边界，同时检查相关请求与实际可见行为，不让所有读者打开页面就下载可选资源。
