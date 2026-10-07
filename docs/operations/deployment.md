# GitHub Pages deployment / GitHub Pages 部署

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

In the repository, go to Settings → Pages → Build and deployment → Source and select **GitHub Actions**. Pushes to `main` or `master` then trigger `.github/workflows/deploy.yml` to build, check, and deploy `dist/`. Pull requests into these branches also build and check the site, without publishing. The workflow installs the pinned Chromium and requires unit and browser checks to pass before publishing. Logs, suite results, and available accessibility or layout details are retained as artifacts for 14 days. Build and deployment jobs are limited to 30 and 10 minutes respectively; Pages write permissions belong only to the deployment job. Manual runs are available through `workflow_dispatch`, and only `main` or `master` can publish. Local edits do not push or publish themselves.

在仓库 Settings → Pages → Build and deployment → Source 中选择 **GitHub Actions**。随后向 `main` 或 `master` 推送，会触发 `.github/workflows/deploy.yml` 自动构建、检查并部署 `dist/`。合入这些分支的 PR 也会构建和检查，但不会发布。工作流安装锁定版本的 Chromium，单元测试和浏览器检查全部通过后才发布；日志、逐项结果及可用的无障碍或排版明细会作为附件保留 14 天。构建和部署分别限时 30 分钟、10 分钟，只有部署任务具有 Pages 写权限。也支持通过 `workflow_dispatch` 手动运行，只有 `main` 或 `master` 分支允许发布。本地修改不会自行推送或发布。

The site uses the root path `/` and targets `https://dapengfeng.github.io`. A project site under a subpath requires a consistent path prefix. Historical date-based article URLs retain redirects.

站点使用根路径 `/`，目标地址为 `https://dapengfeng.github.io`。若改成项目子路径站点，需要统一配置路径前缀。历史日期式文章链接保留跳转。

See [testing](../development/testing.md) for failed check reports, [search configuration](discoverability.md) for Actions verification variables, and [Web Analytics](analytics.md) for the optional Cloudflare beacon variable.

检查失败时的报告位置见[测试与验证](../development/testing.md)，Actions 中的站点验证变量见[搜索配置](discoverability.md)，可选的 Cloudflare 统计变量见[访问统计](analytics.md)。
