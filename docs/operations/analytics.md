# Cloudflare Web Analytics / Cloudflare 访问统计

[Documentation index / 文档索引](../README.md) · [Deployment / 部署](deployment.md)

The blog stays on GitHub Pages. Cloudflare hosts the analytics service; no server, DNS change, or Cloudflare proxy is required. In Cloudflare **Web Analytics**, add the hostname `dapengfeng.github.io`, then open **Manage site** and copy the `token` value from the generated `data-cf-beacon` snippet.

博客继续使用 GitHub Pages，统计服务由 Cloudflare 托管，无需另建服务器、修改 DNS 或开启 Cloudflare 代理。在 Cloudflare **Web Analytics** 添加 `dapengfeng.github.io`，然后从 **Manage site** 生成的代码中复制 `data-cf-beacon` 内的 `token` 值。

In GitHub repository **Settings → Secrets and variables → Actions → Variables**, create `CLOUDFLARE_WEB_ANALYTICS_TOKEN`. Its value must be the actual 32-character hexadecimal token, without quotes or the surrounding snippet. This is a public beacon identifier, not an account API credential. The deployment workflow passes it to the build, which adds the configured loader to all shared page layouts. Re-run the deployment workflow after changing the variable; saving an Actions variable alone does not update deployed HTML.

在 GitHub 仓库 **Settings → Secrets and variables → Actions → Variables** 新建 `CLOUDFLARE_WEB_ANALYTICS_TOKEN`，值填写真实的 32 位十六进制 token，不带引号或整段脚本。这是公开的统计站点标识，不是账号 API 密钥。部署工作流会将其传入构建，并在所有共用页面模板中加入配置好的加载脚本。修改变量后需重新运行部署工作流；只保存变量不会更新已发布的 HTML。

An unset or blank variable emits no analytics loader. A malformed value fails the build with a configuration error. Remove the variable and redeploy to disable collection. Historical redirect pages have no tracker; their destination pages are measured instead.

变量未设置或为空时不输出统计加载脚本；格式错误时构建会明确报错。移除变量并重新部署即可关闭采集。历史跳转页不加载统计脚本，由跳转后的页面记录访问。

The loader requests Cloudflare's module only on the HTTPS hostname configured in `scripts/config.mjs`. Local previews and CI browser checks make no analytics requests. Browsers with Global Privacy Control or Do Not Track enabled are also excluded. Only Cloudflare's standard traffic and performance beacon is loaded; the site does not add user identifiers or custom event payloads.

加载器只在 `scripts/config.mjs` 指定的正式 HTTPS 域名请求 Cloudflare 模块，本地预览和 CI 浏览器检查不会请求统计服务。启用 Global Privacy Control 或 Do Not Track 的浏览器也不会加载。网站仅接入 Cloudflare 标准流量与性能统计，不附加用户身份或自定义事件数据。

Keep reporting aggregate: do not expose individual readers or their IP locations on the site. Do not add IP lookups for recommendations or payment routing. [Language selection](../authoring/language.md) remains local to the browser; service configuration and optional loading stay independent of visitor location.

统计保持汇总形式，不在站点上展示个体读者或其 IP 归属，不为推荐或付款分流添加 IP 查询。[语言选择](../authoring/language.md)仍在浏览器本地完成，服务配置与按需加载不依赖访客位置。

After deployment, visit the production site and check the Cloudflare dashboard after a few minutes. Network failures and blocking extensions can prevent collection; mainland China connectivity must be tested on the reader's actual network. **Visits counts visits, not deduplicated people.** See the [official setup guide](https://developers.cloudflare.com/web-analytics/get-started/), [metric definitions](https://developers.cloudflare.com/web-analytics/data-metrics/high-level-metrics/), and [FAQ](https://developers.cloudflare.com/web-analytics/faq/).

部署后访问正式网站，等待几分钟再查看 Cloudflare 后台。网络故障和拦截扩展可能使统计缺失，大陆连通性需在读者实际网络测试。**Visits 表示访问次数，不是去重人数。** 配置与统计口径见上述官方文档。

Run `npm run build` and `npm run test:analytics` to verify the integration. The browser suite intercepts every request and uses a mock beacon, so it checks loading and configuration without reporting test visits or verifying a real Cloudflare account.

运行 `npm run build` 和 `npm run test:analytics` 验证接入。浏览器测试拦截全部请求并模拟统计模块，用于检查加载和配置，不上报测试访问，也不验证真实 Cloudflare 账号中的数据。
