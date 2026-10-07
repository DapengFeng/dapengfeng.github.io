# Website documentation / 网站文档

[Website mission / 网站使命](../README.md)

Use these guides to write articles, maintain the site, and configure publishing. Examples of article markup remain in [examples/article.html](../examples/article.html); the guides explain when and how to use them. Commands run from the repository root.

这些文档用于文章写作、网站维护与发布配置。文章结构示例保留在 [examples/article.html](../examples/article.html)，各项规范说明其使用方式与条件。所有命令均在仓库根目录执行。

Agents start with the root [AGENTS.md](../AGENTS.md), then use the owning guide below. Writing covers content decisions; presentation covers design and markup; publishing covers metadata and selection logic. Feature-specific behavior belongs in its feature guide.

Agent 从根目录 [AGENTS.md](../AGENTS.md)开始，再阅读下列归属文档。内容决策归写作，页面设计与标记归排版，元信息与取文逻辑归发布；专属功能的行为归对应功能指南。

## 1. Content authoring / 内容创作

| Guide / 文档 | Scope / 内容 |
| --- | --- |
| [Content and writing / 内容与写作](authoring/writing.md) | Mission, content structure, titles, summaries, mechanisms, literary observation, examples and revision.<br>使命、内容组织、标题摘要、机制讲解、文学观察、修改示例与审稿。 |
| [Articles, series and recommendations / 文章、专题与推荐](authoring/publishing.md) | HTML workflow, metadata, dates, categories, series navigation, and actual recommendation rules.<br>HTML 工作流、元信息、日期、分类、专题导航与实际推荐规则。 |
| [Bilingual content and reading / 双语内容与阅读](authoring/language.md) | English–Chinese pairing, shared figures, language selection, and translation maintenance.<br>中英对应、图示共用、语言选择与译文维护。 |
| [Page design and presentation / 页面设计与排版](authoring/presentation.md) | Page roles, hierarchy, responsive layout, figures, controls, formulas, code and compiler examples.<br>页面职责、层级、响应式布局、图示、控件、公式、代码与编译示例。 |
| [Daily mathematics / 每日数学](authoring/daily-mathematics.md) | Reviewed topic library, dated publishing, permanent history and inventory checks.<br>核对后的题库、按日发布、永久历史与库存检查。 |

## 2. Development and validation / 开发与验证

| Guide / 文档 | Scope / 内容 |
| --- | --- |
| [Local workflow and source map / 本地工作流与源码索引](development/local-workflow.md) | Installation, preview, build output, and implementation locations.<br>安装、预览、构建产物与实现位置。 |
| [Testing and validation / 测试与验证](development/testing.md) | Test commands, browser setup, reports, series checks, and verification limits.<br>测试命令、浏览器安装、报告、专题检查与验证范围。 |

## 3. Deployment and services / 部署与服务

| Guide / 文档 | Scope / 内容 |
| --- | --- |
| [GitHub Pages deployment / GitHub Pages 部署](operations/deployment.md) | Repository settings, Actions, branches, permissions, and historical URLs.<br>仓库设置、Actions、分支、权限与历史链接。 |
| [Search and AI discoverability / 搜索与 AI 可发现性](operations/discoverability.md) | Metadata, citations, SEO/GEO, verification variables, and sitemap submission.<br>元信息、引用、SEO/GEO、验证变量与站点地图提交。 |
| [Cloudflare Web Analytics / Cloudflare 访问统计](operations/analytics.md) | Hosted analytics, Actions variable, production-only collection, and privacy signals.<br>托管统计、Actions 变量、正式站点采集与隐私信号。 |
| [GitHub Discussions / GitHub 讨论](operations/discussions.md) | giscus setup, one thread per article, loading, login, and draft behavior.<br>giscus 配置、每文一个话题、加载、登录与草稿行为。 |
| [Reader support / 读者赞赏](operations/support.md) | PayPal support panel and Actions variable configuration.<br>PayPal 赞赏面板及 Actions 变量配置。 |
| [Article sharing / 文章分享](operations/sharing.md) | Platform sharing, QR codes, image generation, language selection, and compatibility.<br>平台分享、二维码、图片生成、语言选择与兼容性。 |

## Keep the guides current / 文档维护

Update the relevant guide when a workflow or shared component changes. Keep English and Chinese together, update the HTML example when markup changes, and add new guides to this index. These files are repository documentation; the site build publishes articles from `content/posts/`, not the `docs/` directory.

工作流或共用组件变化时，同步更新对应文档的英文与中文。标记写法变化时更新 HTML 示例，新增文档时补充本索引。这些文件是仓库文档；站点构建发布 `content/posts/` 中的文章，不发布 `docs/` 目录。

Keep one full statement of each rule in its owning guide; use links elsewhere. When merging a guide, move its unique requirements first, update inbound links and remove the superseded file. Resolve conflicts against the latest explicit user decision and the implementation; do not record an unimplemented suggestion as behavior. Commands and checks belong in development guides rather than being copied into every authoring page.

每项规则只在归属文档中完整描述，其他位置使用链接。合并时先迁移独有要求，更新所有引用，再删除被替代的文件。冲突应结合用户最新明确决定与实现核对，不把尚未实施的建议记作已有行为。命令与检查集中在开发指南，不在每份写作规范里重复。
