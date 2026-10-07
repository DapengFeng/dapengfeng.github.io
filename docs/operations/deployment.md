# GitHub Pages deployment / GitHub Pages 部署

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

In the repository, go to Settings → Pages → Build and deployment → Source and select **GitHub Actions**. Pushes to `main` or `master` and pull requests into them run `deploy.yml`. `checks.yml` is the shared implementation: it selects checks from changed files, builds once, runs all unit tests, and distributes selected browser suites across up to three parallel jobs. The final **CI result** job requires every selected job to succeed. Pull requests never publish; documentation-only changes do not build or publish.

在仓库 Settings → Pages → Build and deployment → Source 中选择 **GitHub Actions**。向 `main` 或 `master` 推送，或向它们提交 PR，会触发 `deploy.yml`。共用的 `checks.yml` 根据改动文件选择检查范围，只构建一次，运行全部单元测试，并将选中的浏览器测试分给最多三个并行任务。最后的 **CI result** 要求本次选中的任务全部成功。PR 不发布，纯文档改动不构建也不发布。

| Change or trigger / 改动或触发方式 | Checks / 检查范围 |
| --- | --- |
| `README.md`, root `AGENTS.md`, `docs/**/*.md` only / 仅文档（含根目录 `AGENTS.md`） | UTF-8, local link targets, code fences, merge conflicts; no dependency installation.<br>检查编码、本地链接目标、代码围栏和冲突标记，无需安装依赖。 |
| Known articles or their photos / 已知文章或配图 | Full build and unit tests; reading, language, layout, accessibility, formulas, copying, sharing, and relevant article suites.<br>完整构建与单元测试；阅读、语言、排版、无障碍、公式、复制、分享及相关专题检查。 |
| Analytics script / 统计脚本 | Full build and unit tests; basic browser, language, reading, and analytics suites.<br>完整构建与单元测试；基础浏览、语言、阅读及统计检查。 |
| Daily mathematics content/code or daily timer / 每日数学内容、代码或每日定时任务 | Full build, unit tests, history/inventory checks; basic browsing/language/reading plus the daily-math suite, which covers formulas and responsive layouts.<br>完整构建与单元测试、历史与库存检查；基础浏览、语言、阅读及每日数学套件，后者包含公式与响应式排版验证。 |
| Shared styles/templates, dependencies, workflows, unknown paths or unavailable diff / 公共样式、模板、依赖、工作流、未知路径或无法比较改动 | All checks / 全部检查。 |

The daily publication is scheduled for **00:23 Asia/Shanghai**. `regression.yml` runs all suites every **Sunday at 02:17 Asia/Shanghai**, without deployment, and also supports manual runs from Actions → Full site regression → Run workflow. A manual run of Build and deploy knowledge lab runs all checks and can publish from `main` or `master`. GitHub may delay scheduled jobs; these are target times, not publication guarantees.

每日发布计划在北京时间 **00:23** 运行。`regression.yml` 每周日北京时间 **02:17** 执行完整回归，不部署，也可从 Actions → Full site regression → Run workflow 手动触发。手动运行 Build and deploy knowledge lab 会进行全量检查，且仅 `main` 或 `master` 允许发布。GitHub 的定时任务可能延迟，这些时间是计划时间，不保证准点发布。

Browser jobs download the same `checked-site` artifact. After their checks pass, deployment unpacks that artifact and uploads it to Pages, without rebuilding. The archive preserves hidden files such as `.nojekyll`. Site artifacts are retained for one day; per-group logs and available accessibility reports for 14 days. Each suite has a five-minute timeout, build jobs 15 minutes, browser groups 25 minutes, and deployment 10 minutes. Only deployment has Pages write permissions. Local edits do not push or publish themselves.

浏览器任务下载同一份 `checked-site` 构建产物。检查通过后，部署任务解包该产物并上传到 Pages，不重复构建；归档会保留 `.nojekyll` 等隐藏文件。站点产物保留一天，各组日志及可用的无障碍报告保留 14 天。每个测试套件限时五分钟，构建任务 15 分钟，浏览器分组 25 分钟，部署任务 10 分钟。只有部署任务具有 Pages 写权限。本地修改不会自行推送或发布。

For branch protection, require the deployment workflow's stable **CI result** check shown in the PR, rather than individual browser matrix jobs that can be skipped on documentation changes. There are no workflow-level path filters. Test selection and the reason for a full fallback are shown in the run summary. The mapping lives in `scripts/ci-plan.mjs`; unrecognized new articles get full coverage until explicitly mapped.

设置分支保护时，应选择 PR 中发布工作流的固定 **CI result** 汇总检查，不要要求文档改动时会跳过的单个浏览器分组。工作流本身不使用路径过滤。运行摘要会显示测试范围和回退全量检查的原因。映射维护在 `scripts/ci-plan.mjs` 中；未登记的新文章会先运行完整检查。

The site uses the root path `/` and targets `https://dapengfeng.github.io`. A project site under a subpath requires a consistent path prefix. Historical date-based article URLs retain redirects.

站点使用根路径 `/`，目标地址为 `https://dapengfeng.github.io`。若改成项目子路径站点，需要统一配置路径前缀。历史日期式文章链接保留跳转。

See [testing](../development/testing.md) for failed check reports, [search configuration](discoverability.md) for Actions verification variables, and [Web Analytics](analytics.md) for the optional Cloudflare beacon variable.

检查失败时的报告位置见[测试与验证](../development/testing.md)，Actions 中的站点验证变量见[搜索配置](discoverability.md)，可选的 Cloudflare 统计变量见[访问统计](analytics.md)。
