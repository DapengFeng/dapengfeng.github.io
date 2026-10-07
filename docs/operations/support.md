# Reader support / 读者赞赏

## Entry points / 入口

Every blog article has a fixed bottom-right toolbar for sharing and opening the discussion. When `PAYPAL_ME_URL` is configured, a coffee icon opens a compact PayPal panel, without payment-method tabs or unavailable placeholders. The toolbar hides while a panel is open and accommodates narrow screens and device safe areas.

每篇博客右下角提供固定的浮动工具栏，可分享文章或打开讨论。配置 `PAYPAL_ME_URL` 后显示咖啡图标，打开紧凑的 PayPal 面板，不显示支付方式切换或未开通提示。面板打开时工具栏隐藏，并适配窄屏和设备安全区域。

## Configuration / 配置

In repository **Settings → Secrets and variables → Actions → Variables**, create the repository variable `PAYPAL_ME_URL` with your complete HTTPS `paypal.me` profile URL. The shared [checks workflow](../../.github/workflows/checks.yml) passes it to the build for deployment and regression runs. Leave it unset or empty to hide the support entry point. Do not append an amount, query string or redirect. This is a public payment link, visible in the generated website; it is not an API credential and does not need to be a Secret.

在仓库 **Settings → Secrets and variables → Actions → Variables** 中创建 Repository variable `PAYPAL_ME_URL`，值为完整的 HTTPS `paypal.me` 个人主页链接。共用的[检查工作流](../../.github/workflows/checks.yml)会将变量传入部署和回归构建。未配置或留空则隐藏赞赏入口，不附加金额、查询参数或跳转地址。这是公开收款链接，会出现在生成的网站中，不是 API 凭据，无需使用 Secret。

The build rejects invalid links. After changing the variable, run **Build and deploy knowledge lab** from Actions to rebuild and publish; saving the variable alone does not change the deployed site.

构建时检查链接格式。修改变量后，在 Actions 手动运行 **Build and deploy knowledge lab** 以重新构建并发布；仅保存变量不会改变线上网站。

For local previews, put `PAYPAL_ME_URL=https://paypal.me/yourprofile` in the ignored root `.env.local`. `npm run build` and `npm run dev` load it; an existing environment variable takes precedence. Keep this file out of commits.

本地预览可在仓库根目录的 `.env.local` 中填写 `PAYPAL_ME_URL=https://paypal.me/yourprofile`。`npm run build` 和 `npm run dev` 自动读取，已有环境变量优先。该文件已被 Git 忽略，不应提交。

## Behavior / 行为

The panel offers US$1, US$3 (initially selected), US$5 and an always-visible custom input starting at US$10. Clicking or editing that input selects the custom amount directly; tabbing through it leaves the selected amount unchanged. It accepts positive amounts with up to two decimal places. The selected amount and explicit `USD` currency are appended to the profile link; invalid custom amounts cannot open checkout.

面板提供 US$1、US$3（初始选中）、US$5 和始终可见的自定义输入框，初始显示 US$10。点击或编辑输入框即可直接使用自定义金额；仅按 Tab 经过输入框不会改变已选金额。可填写最多两位小数的正数。选定金额与明确的 `USD` 币种会附加到个人主页链接；自定义金额无效时无法前往付款。

PayPal opens its official page in a new tab. Payment confirmation remains on PayPal’s side. This feature adds no third-party scripts, geolocation requests, transaction database or payment-success claims. With JavaScript disabled, a direct PayPal link remains available after the article.

PayPal 在新标签页打开，付款确认由 PayPal 完成。此功能不添加第三方脚本、地理位置查询、交易数据库，也不显示无法核验的“付款成功”。关闭 JavaScript 后，正文末尾仍有 PayPal 直达链接。

The custom input keeps the same outer dimensions as the presets, with a recessed number area and underline. Its initial `10` is selected on focus/click for replacement; subsequent edits retain normal caret behavior. Mobile devices use a decimal keyboard.

自定义输入框与预设金额保持相同外部尺寸，数字区略微内凹并带下划线。初始的 `10` 在聚焦或点击时全选，方便直接替换；编辑后恢复普通光标操作。移动设备使用小数键盘。

## Verification / 验证

```sh
npm run build
node --test tests/support.test.mjs
npm run test:support
npm run test:sharing
npm run test:discussions
```

The browser suites check the payment target, floating toolbar, keyboard focus, three language modes, narrow screens, accessibility and the no-JavaScript fallback. They use a test-only profile injected into the served page, so forks and unconfigured builds exercise the same UI without publishing a sample payment link or contacting PayPal.

浏览器检查覆盖收款链接、浮动工具栏、键盘焦点、三种语言模式、窄屏、无障碍和无 JavaScript 的回退。测试只在测试页面注入示例收款链接，因此 fork 和未配置变量的构建也能验证完整界面，不会发布示例链接，也不联系 PayPal。
