# Website documentation / 网站文档

[Website mission / 网站使命](../README.md)

Use these guides to write articles, maintain the site, and configure publishing. Examples of article markup remain in [examples/article.html](../examples/article.html); the guides explain when and how to use them. Commands run from the repository root.

这些文档用于文章写作、网站维护与发布配置。文章结构示例保留在 [examples/article.html](../examples/article.html)，各项规范说明其使用方式与条件。所有命令均在仓库根目录执行。

## 1. Content authoring / 内容创作

| Guide / 文档 | Scope / 内容 |
| --- | --- |
| [Articles and learning series / 文章与学习专题](authoring/publishing.md) | HTML workflow, metadata, dates, categories, series navigation, and related reading.<br>HTML 工作流、元信息、日期、分类、专题导航与相关阅读。 |
| [Bilingual content and reading / 双语内容与阅读](authoring/language.md) | English–Chinese pairing, shared figures, language selection, and translation maintenance.<br>中英对应、图示共用、语言选择与译文维护。 |
| [Presentation and interactive examples / 排版与交互示例](authoring/presentation.md) | Typography, links, AMS equations, numbering, code blocks, copying, Godbolt, and diagrams.<br>字号、链接、AMS 公式、编号、代码块、复制、Godbolt 与示意图。 |
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
| [GitHub Discussions / GitHub 讨论](operations/discussions.md) | giscus setup, article mapping, message types, login, and draft behavior.<br>giscus 配置、文章关联、发言类型、登录与草稿行为。 |
| [Article sharing / 文章分享](operations/sharing.md) | Platform sharing, QR codes, image generation, language selection, and compatibility.<br>平台分享、二维码、图片生成、语言选择与兼容性。 |

## Keep the guides current / 文档维护

Update the relevant guide when a workflow or shared component changes. Keep English and Chinese together, update the HTML example when markup changes, and add new guides to this index. These files are repository documentation; the site build publishes articles from `content/posts/`, not the `docs/` directory.

工作流或共用组件变化时，同步更新对应文档的英文与中文。标记写法变化时更新 HTML 示例，新增文档时补充本索引。这些文件是仓库文档；站点构建发布 `content/posts/` 中的文章，不发布 `docs/` 目录。
