# Bilingual content and reading / 双语内容与阅读

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

On a first visit, the browser’s preferred language determines the edition: Chinese language tags use Chinese; other languages use English. Readers can choose English, Chinese, or both, and their saved manual choice always takes priority. In bilingual mode, English is followed immediately by Chinese.

首次访问时按浏览器的首选语言显示：中文语言标记使用中文，其他语言使用英文。读者可选择英文、中文或双语，已保存的手动选择始终优先。双语模式中英文在前，对应中文紧随其后。

Language is resolved before the page is painted, using the saved `feng-language` preference or `navigator.languages`. No location service or IP lookup is involved, and no delayed network response changes the chosen edition. Old automatic-language session values are ignored. If storage is blocked, the browser preference still works and a manual switch applies to the current page. Repository documentation always displays English followed by Chinese.

语言在页面绘制前确定，使用已保存的 `feng-language` 偏好或 `navigator.languages`。判断过程不访问定位服务，不查询 IP，也不会因稍后返回的网络结果改变语言；旧版本会话中保存的自动语言结果会被忽略。存储被禁用时仍按浏览器语言显示，手动切换对当前页面有效。仓库文档始终按英文在前、中文紧随其后的顺序显示。

Reading time follows the selected edition. The build estimates English at 220 words per minute and Chinese at 400 characters per minute; shared prose and code are counted once in bilingual mode. Page metadata, navigation, controls, scripts, and hidden demo state do not contribute to the estimate. Language fragments carry native `lang` attributes for assistive technology.

阅读时长随所选语言更新。构建按每分钟 220 个英文单词、400 个汉字估算，双语模式中的共用正文与代码只计一次。页面元数据、导航、操作控件、脚本和隐藏的演示状态不计入阅读时长。各语言片段带有原生 `lang` 属性，便于辅助技术正确朗读。

Interface text comes from `scripts/i18n.mjs`; article titles, summaries, and both body languages live in each article HTML. The build preserves the order you write, derives the bilingual contents from headings, and renders formulas. For shared interactive experiments, keep any language dictionary inside the same HTML, as Spike Notes does.

界面文字由 `scripts/i18n.mjs` 提供；文章标题、摘要与双语正文都在各自的 HTML 中。构建保留你编写的顺序，从标题提取双语目录并渲染公式。共享交互实验所需的语言词典也放在同一 HTML 内，脉冲长文已采用这种方式。

Maintain English and Chinese together in the same file when editing. The build does not translate or assess translation accuracy. Existing articles pair their full explanations in English and Chinese; formulas and code are shared when appropriate. New articles must include both language markers; copy the bilingual template to begin.

修改时在同一文件中同步维护英文和中文。构建不会翻译或判断译文准确性。现有文章的英文与中文完整说明成对排列，公式与代码按需共用。新文章必须包含两种语言标记，可直接复制双语模板开始编写。
