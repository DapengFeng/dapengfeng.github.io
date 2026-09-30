# Bilingual content and reading / 双语内容与阅读

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

On a first visit, IP countries/regions CN, HK, MO, and TW default to Chinese; other locations default to English. Readers can choose English, Chinese, or both, and their saved choice always takes priority. In bilingual mode, English is followed immediately by Chinese.

首次访问时，IP 所在国家或地区为 CN、HK、MO、TW 时默认中文，其余默认英文。读者可选择英文、中文或双语，已保存的手动选择始终优先。双语模式中英文在前，对应中文紧随其后。

Automatic language detection calls [Country](https://country.is/), which receives the visitor’s network IP. The site stores only the selected language, not the IP or country response. Automatic results are cached for the tab session; a saved manual preference skips the lookup. Browser language is used immediately and remains the fallback if the request fails or exceeds 2.5 seconds. VPNs may affect the country result. Repository documentation always displays English followed by Chinese.

自动语言判断会请求 [Country](https://country.is/)，服务方会收到访客的网络 IP。本站只保存选中的语言，不保存 IP 或地区查询响应。自动结果在标签页会话中缓存；已有手动偏好时跳过查询。页面先按浏览器语言显示，查询失败或超过 2.5 秒则继续使用该语言。VPN 可能影响地区结果。仓库文档始终按英文在前、中文紧随其后的顺序显示。

Interface text comes from `scripts/i18n.mjs`; article titles, summaries, and both body languages live in each article HTML. The build preserves the order you write, derives the bilingual contents from headings, and renders formulas. For shared interactive experiments, keep any language dictionary inside the same HTML, as Spike Notes does.

界面文字由 `scripts/i18n.mjs` 提供；文章标题、摘要与双语正文都在各自的 HTML 中。构建保留你编写的顺序，从标题提取双语目录并渲染公式。共享交互实验所需的语言词典也放在同一 HTML 内，脉冲长文已采用这种方式。

Maintain English and Chinese together in the same file when editing. The build does not translate or assess translation accuracy. Existing articles pair their full explanations in English and Chinese; formulas and code are shared when appropriate. New articles must include both language markers; copy the bilingual template to begin.

修改时在同一文件中同步维护英文和中文。构建不会翻译或判断译文准确性。现有文章的英文与中文完整说明成对排列，公式与代码按需共用。新文章必须包含两种语言标记，可直接复制双语模板开始编写。
