# Article sharing / 文章分享

[Documentation index / 文档索引](../README.md)

Every blog article has a share icon in the fixed bottom-right toolbar, alongside discussion and, when configured, [reader support](support.md). The dialog offers WeChat, RedNote, X, LinkedIn, Telegram, link/text copying, and a downloadable image. The system share option appears when the browser supports it; the available destination apps depend on the device. Native file sharing appears after a local poster has loaded and the browser reports support. Opening a share sheet does not mean a post has been published.

每篇 blog 右下角的浮动工具栏提供分享图标，与讨论及配置后的[赞赏入口](support.md)并列。面板提供微信、小红书、X、LinkedIn、Telegram、链接／文案复制和图片下载。浏览器支持时显示系统分享，可选目标应用由设备决定。分享图加载完成且浏览器支持文件分享后，才显示图片系统分享按钮。打开系统分享面板不代表已经发表内容。

WeChat uses a scannable article QR code and a saved sharing image. RedNote uses a saved image plus editable, copyable text. These are manual publishing workflows; this implementation does not call platform SDKs, authenticate to social accounts, or submit posts automatically. X, LinkedIn, and Telegram open their own sharing pages for the reader to review and submit. There are no third-party sharing scripts or QR services.

微信通过文章二维码和保存的分享图完成分享；小红书通过保存图片及编辑、复制文案来准备笔记。这些流程由读者手动发布；本实现不调用平台 SDK、不登录社交账号，也不自动发帖。X、LinkedIn 和 Telegram 打开各自的分享页面，由读者检查并提交。站点不加载第三方分享脚本，也不使用外部二维码服务。

Titles, summaries, images, and labels follow the language selection: English, Chinese, or English followed by Chinese. Edited text is retained separately for each language mode while the page remains open. Link copying, QR codes, and platform links always use the article’s production URL from `site.url`, excluding local preview addresses, hashes, and query parameters. A newly written article becomes reachable at that URL only after deployment.

标题、摘要、图片与标签遵循语言选择：英文、中文，或英文紧接中文。页面保持打开时，各语言模式分别保留修改后的文案。复制链接、二维码与平台链接始终使用 `site.url` 下的文章正式地址，不带本地预览地址、锚点或查询参数。新文章需要部署之后，正式地址才可访问。

## Build and maintenance / 构建与维护

The build generates `/assets/share/<slug>-qr.png` and three 900 × 1200 PNG posters ending in `-en.png`, `-zh.png`, and `-both.png`. QR encoding runs locally with `qrcode`; `sharp` renders article covers or existing illustrations with titles, summaries, author, date, and QR. Chinese text requires Noto Sans CJK SC. GitHub Actions installs `fonts-noto-cjk` before building; on Debian/Ubuntu, install it with `sudo apt-get install fonts-noto-cjk` for local builds. No browser-side framework is needed. Images load only after the share dialog is opened, and a selected QR/poster view loads its preview on demand.

构建自动生成 `/assets/share/<slug>-qr.png` 和以 `-en.png`、`-zh.png`、`-both.png` 结尾的三张 900 × 1200 PNG 分享图。二维码由本地 `qrcode` 编码；`sharp` 将文章封面或已有插图与标题、摘要、作者、日期、二维码合成。中文需要 Noto Sans CJK SC 字体。GitHub Actions 在构建前安装 `fonts-noto-cjk`；本地 Debian／Ubuntu 可执行 `sudo apt-get install fonts-noto-cjk` 安装。不依赖浏览器端框架。打开分享面板后才会加载图片，二维码／海报预览在选中相应平台时加载。

`npm test` independently decodes the QR images and representative posters, checks all article assets, and tests escaping. `npm run test:sharing` checks language selection, canonical URLs, draft preservation, clipboard denial, native share cancellation/failure, downloads, keyboard focus, and mobile layout. Native sharing is mocked: these checks do not publish to platforms or certify the behavior of a particular WeChat/RedNote app version.

`npm test` 独立解码二维码与代表性海报，检查所有文章的资源及转义。`npm run test:sharing` 检查语言选择、正式链接、文案保留、剪贴板拒绝、系统分享取消／失败、下载、键盘焦点及手机排版。系统分享使用模拟接口：这些检查不会向平台发帖，也不代表已经验证特定微信／小红书客户端版本的行为。
