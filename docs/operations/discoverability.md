# Search and AI discoverability / 搜索与 AI 可发现性

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

Every build generates unique English-first bilingual titles and descriptions, canonical URLs, Open Graph and Twitter cards, and 1200 × 630 PNG link-preview images in `/assets/social/`. Article preview images use the English title as graphical text. The separate downloadable posters are covered by the [sharing guide](sharing.md). Author, publication and revision dates, chapter anchors, and explicitly cited references are connected through Schema.org JSON-LD. The same static HTML contains both languages and the full article before JavaScript runs. The 404 page is marked `noindex`; drafts and redirects stay out of the sitemap. `npm test` checks these properties and runs automatically before deployment.

每次构建自动生成各页独有、英文在前的双语标题与摘要、规范链接、Open Graph 与 Twitter 分享元信息，以及 `/assets/social/` 下的 1200 × 630 PNG 链接预览图；文章预览图以英文标题作为图中文字。可下载海报是另一类资源，见[分享指南](sharing.md)。作者、发布与更新日期、章节锚点和明确标注的参考链接通过 Schema.org JSON-LD 关联。同一个静态 HTML 在 JavaScript 运行前就包含双语完整正文。404 页面标记为 `noindex`，草稿与跳转页不进入站点地图。`npm test` 检查这些规则，并在部署前自动执行。

Keep `titleEn`, `title`, `descriptionEn`, and `description` accurate in each article's metadata. Explain the result, assumptions, methods, and limits in visible bilingual prose, including a text explanation of interactive figures. Link primary sources where claims are made. Add `data-citation` to a source link to include it in structured citations; links inside `#references` containers and `.sources-grid` are also recognized. Only existing links are included. Do not invent sources, measurements, credentials, or update dates. Language buttons are reading modes on one URL; they are not separate `hreflang` editions.

在文章元信息中准确填写 `titleEn`、`title`、`descriptionEn` 与 `description`。用可见的双语正文说明结论、假设、方法和局限，为交互图补充文字解释。在相关论述处链接原始来源；来源链接添加 `data-citation` 后会进入结构化引用，`#references` 容器与 `.sources-grid` 内的链接也会被识别。只收录实际存在的链接，不编造来源、测量结果、资历或更新日期。语言按钮是同一 URL 的阅读模式，不是独立的 `hreflang` 语言版本。

One-time setup: add `https://dapengfeng.github.io/` as a URL-prefix property in Google Search Console and as a site in Bing Webmaster Tools. Choose HTML meta-tag verification and copy only the tag's `content` value into the repository's Actions variables `GOOGLE_SITE_VERIFICATION` and `BING_SITE_VERIFICATION` (Settings → Secrets and variables → Actions → Variables). Deploy, finish verification in each platform, then submit `https://dapengfeng.github.io/sitemap.xml`. Locally, the same environment variables populate the verification tags. Empty variables add no tags. Tokens are public verification identifiers, not account passwords. Account verification and sitemap submission have not been performed by this build.

首次配置：在 Google Search Console 添加 `https://dapengfeng.github.io/` 网址前缀资源，在 Bing Webmaster Tools 添加同一站点。选择 HTML 元标签验证，只把标签的 `content` 值分别填入仓库的 Actions 变量 `GOOGLE_SITE_VERIFICATION` 与 `BING_SITE_VERIFICATION`（Settings → Secrets and variables → Actions → Variables）。部署后在平台完成验证，再提交 `https://dapengfeng.github.io/sitemap.xml`。本地构建也可使用同名环境变量；变量为空时不输出验证标签。这些值是公开验证标识，不是账号密码。构建不会代替账号验证和站点地图提交。

These improvements support crawling, understanding, and attribution; they do not guarantee indexing, rankings, or AI citations. Google says no special GEO schema or `llms.txt` file is required. Monitor actual search queries and indexed pages in the webmaster platforms. See [Google's AI search guidance](https://developers.google.com/search/docs/fundamentals/ai-optimization-guide).

这些改进帮助抓取、理解和来源归属，但不保证收录、排名或 AI 引用。Google 说明无需专用 GEO 标记或 `llms.txt` 文件。请通过站长平台观察真实搜索词与页面收录情况。参见 [Google 的 AI 搜索指南](https://developers.google.com/search/docs/fundamentals/ai-optimization-guide)。
