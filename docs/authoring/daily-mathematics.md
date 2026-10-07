# Daily mathematics / 每日数学

[Documentation index / 文档索引](../README.md)

The homepage publishes one reviewed mathematical topic for each scheduled date in Asia/Shanghai. Dates are absolute; there is no annual reset, modulo selection, or cycling of old topics. A dated archive preserves earlier entries. Use `npm run math:check` for the current inventory and final scheduled date rather than a count copied into documentation.

首页按照北京时间的排期，每天展示一个经过核对的数学主题。排期使用完整日期，不按年度重置、不取模、不轮播旧内容，往期保留在日期归档中。当前库存与排期末日以 `npm run math:check` 输出为准，不在文档中维护一份易过时的数字。

## Files / 文件

| File / 文件 | Role / 用途 |
| --- | --- |
| `content/daily-math/topics.json` | Bilingual topic library, sources, review status and renderer IDs.<br>双语题库、来源、核对状态和绘图标识。 |
| `content/daily-math/publications.json` | Permanent date–topic–concept assignments, including the upcoming queue.<br>永久的日期、主题、知识点对应关系，包含待发布队列。 |
| `scripts/daily-math.mjs` | Validation, history protection, inventory and scheduling.<br>校验、历史保护、库存计算和排期。 |
| `src/scripts/daily-math-models.js` | Renderer registry and numerical models.<br>绘图注册表与数值模型。 |
| `src/scripts/daily-math-drawings.js`, `src/scripts/surface.js` | Formula-driven canvas drawings, date loading and animation lifecycle.<br>按公式绘图、按日期加载和动画生命周期。 |

## Add and schedule a topic / 新增与排期

1. Add a new object to `topics.json` with a permanent, unique `id` and `concept`. The concept identifies the actual mathematical question, not a color palette or a parameter choice. Start with `status: "draft"` while preparing it.

   在 `topics.json` 新增对象，填写永久且唯一的 `id` 和 `concept`。知识点标识对应实际数学问题，不以换配色或改参数作为新知识点。准备期间使用 `status: "draft"`。

2. Supply `en`, `zh`, `descriptionEn`, `descriptionZh`, `readingEn`, `readingZh`, `formula`, and `source: {url, en, zh}`. The description explains the principle; the reading text maps the diagram to that principle. Use `\(...\)` for inline mathematics. `compactFormula` may provide an equivalent narrow-screen layout for a long display equation.

   填写中英标题、原理说明、读图说明、LaTeX 公式和来源。原理说明解释知识，读图说明将画面与知识对应。行内公式使用 `\(...\)`；较长的 display 公式可用 `compactFormula` 提供等价的窄屏排版。

   An optional `related: [{url, en, zh}]` links to an existing, relevant article under `/blog/`. Use a short bilingual label that states what the article covers. Prefer one useful next step; do not imply that a related article proves or explains material it does not contain. These links supplement the primary source and appear on both the homepage and dated entry.

   可选的 `related: [{url, en, zh}]` 指向 `/blog/` 下已有且相关的文章。双语短标题应说明链接内容，优先提供一个有用的下一步，不要暗示文章解释或证明了其中并未涉及的知识。这些链接补充原始来源，在首页与日期页面均显示。

3. Implement a meaningful renderer and register its ID. Check the source, mathematical identities, both translations, continuous animation, and mobile formula layout. Set `status: "ready"` and the actual `reviewedOn` date only after review. The renderer may be reused for a different substantive question; reusing a topic's identity or concept is forbidden.

   实现能表达原理的绘图并注册标识。核对来源、数学关系、双语对应、持续播放的动画，以及移动端公式排版。检查完成后才设为 `status: "ready"`，并填写实际 `reviewedOn` 日期。不同实质问题可以复用绘图代码，但不能重复使用主题或知识点标识。

4. Append ready topics to the schedule, then build and verify:

   将已核对的主题追加到排期，再构建与验证：

   ```sh
   npm run math:plan
   npm run math:check -- --require-today
   npm run build
   npm test
   npm run test:daily-math
   npm run preview
   ```

   `npm run math:plan -- --count 5` schedules only five unscheduled ready topics. Existing assignments stay unchanged. If the queue has expired, new entries start today, leaving a truthful gap rather than inventing past publications. Commit the topic library, drawings and publication ledger together.

   `npm run math:plan -- --count 5` 只追加五个未排期的就绪主题。已有记录保持不变。若队列已经用完，新内容从今天开始，保留真实的空缺日期，不补造过去的发布记录。题库、绘图和发布记录应一起提交。

## What is checked / 检查范围

Duplicate topic IDs, concept IDs, titles, descriptions, publication dates and reused published topics fail validation. Similar English descriptions raise a review warning. Changing or deleting an already published date–topic–concept assignment fails CI when compared with the previous commit or PR base. Corrections to the explanation or drawing keep the original date and URL. The initial migration has no previous ledger to compare.

重复的主题标识、知识点标识、标题、原理说明、发布日期，以及重复发布同一主题，都会被拒绝。相似的英文说明会提示人工复核。CI 对照上一次提交或 PR 基线，阻止修改、删除已发布的日期、主题和知识点对应关系。修正说明或绘图时保留原日期和链接。初次迁移尚无旧发布记录可供对照。

Automated text checks cannot prove semantic novelty or mathematical correctness. Reviewers must reject renamed duplicates and superficial variants. There is no year-count limit on the ledger; sustaining five or ten years requires continuing to prepare distinct, correct topics. The site does not automatically scrape articles or publish unreviewed AI-generated mathematics.

文本检查不能证明知识点在语义上全新，也不能替代数学正确性核对。审核时应排除改名重复和表面变体。发布记录不设年数上限；持续五年、十年仍需不断补充准确且有实质差异的内容。网站不会自动抓取文章或发布未经核对的 AI 数学内容。

## Monthly preparation / 按月准备

At the beginning of each month, run `npm run math:check` and inspect the final scheduled date. Prepare the next month's distinct topics as drafts, review the formulas, sources, translations and diagrams, and only then use `npm run math:plan` to append the reviewed batch. Aim to finish replenishment while at least 14 future days remain; the existing daily Actions warning is the fallback reminder, not a substitute for review. Check the queue again after scheduling and retain the original dates, IDs and concepts of published entries.

每月初运行 `npm run math:check`，查看排期的最后日期。提前为下个月准备不重复的草稿，核对公式、来源、翻译和图示后，再用 `npm run math:plan` 追加已审核批次。目标是在未来库存仍有至少 14 天时完成补充；已有的每日 Actions 库存警告用于兜底提醒，不能代替审核。排期后再次检查库存，保留已发布条目的日期、标识与知识点。

## Daily operation / 每日运行

Each prepared date has a small same-origin JSON asset containing its pre-rendered MathJax/AMS equations and explanation. The homepage contains only one topic and requests only the current date; it does not download the entire library. It switches at Shanghai midnight, checks again after sleep, and retains at most four topic payloads in memory. A dated archive page stays on its original topic. Without JavaScript, the homepage shows the build-date entry with its actual date.

每个准备好的日期都有一个同源 JSON 文件，包含预先排版的 MathJax/AMS 公式与说明。首页只内嵌一个主题、按需获取当天文件，不下载整个题库。北京时间零点切换，休眠恢复时再次检查，内存最多保留四项。日期归档页面始终展示原主题。禁用 JavaScript 时，首页展示构建当天的内容并标明真实日期。

The mathematical backgrounds on the homepage and dated entries animate continuously at up to 30 fps while the scene and browser tab are visible. They pause offscreen or in a hidden tab and resume on return. There are no Appearance settings: system reduced-motion preferences and previously saved motion choices do not change this behavior. Animation requires no browser storage. Without JavaScript, equations, explanations and reading links remain available.

首页和日期页面的数学背景在画面与浏览器标签页可见时持续播放，帧率不超过每秒 30 帧。画面离屏或标签页隐藏时暂停，返回后恢复。页面没有“显示设置”；系统的减少动态效果偏好和以前保存的动态选择均不影响此行为。动画不依赖浏览器存储。禁用 JavaScript 时，公式、说明和延伸阅读链接仍可使用。

If no entry exists for today or loading fails, the homepage removes the old topic and offers the archive. It never labels yesterday's mathematics as today's. Transient failures are retried on return to the tab and hourly. No visitor identifiers or analytics are required.

当天没有内容或加载失败时，首页移除旧主题并显示往期入口，不会把昨天的内容标成今天。临时故障会在返回页面或每小时检查时重试。整个功能不需要访客标识或访问统计。

The [daily deployment workflow](../operations/deployment.md) refreshes static pages, the archive and sitemap, and checks inventory. Fewer than 14 future days produces an Actions warning and job summary; no entry today fails the inventory step. There is no separate notification service. Prepared topic changes at midnight do not depend on this workflow finishing on time.

现有[每日部署工作流](../operations/deployment.md)刷新静态页面、归档与站点地图，并检查库存。未来不足 14 天时产生 Actions 警告和任务摘要；当天没有内容时库存检查失败。未接入额外通知服务。已准备好的主题在零点切换，不依赖工作流准时完成。

GitHub schedules can be delayed, and scheduled workflows in public repositories are disabled after 60 days without repository activity. Keep the repository active as content is replenished; manually run the workflow after re-enabling it if necessary. See [GitHub's schedule documentation](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#schedule).

GitHub 定时任务可能延迟，公开仓库连续 60 天没有活动时，定时工作流会被停用。持续补充内容也会保持仓库活跃；必要时重新启用并手动运行工作流。参见 [GitHub 定时工作流文档](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#schedule)。
