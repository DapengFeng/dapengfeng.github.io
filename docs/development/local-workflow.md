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

```bash
npm run build
npm run preview
npm test
npm run test:browsers
```

See [testing and validation](testing.md) for browser installation, test coverage, and report locations.

浏览器安装、测试范围与报告位置见[测试与验证](testing.md)。

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

- `src/scripts/surface.js`, `src/scripts/daily-math*.js`, `src/styles/daily-math.css`

  Date-scheduled mathematical backgrounds and responsive annotations. The canvas plays automatically at up to 30 fps, honors reduced motion, and stops off screen. Formula SVGs are typeset with MathJax/AMS during the build. The homepage loads one dated topic; permanent links and yearly archives preserve history. Read [Daily mathematics](../authoring/daily-mathematics.md) for the topic library, scheduling commands, duplicate checks, and replenishment workflow.

  按日期排期的数学背景与响应式原理注释。Canvas 自动播放，帧率不超过每秒 30 帧，尊重减少动态效果设置，在屏幕外停止。公式在构建时通过 MathJax/AMS 排版为 SVG。首页只加载一个日期主题，通过永久链接和年度归档保留历史。题库、排期命令、重复检查与补充工作流见[每日数学](../authoring/daily-mathematics.md)。

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
