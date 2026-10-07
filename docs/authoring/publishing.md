# Articles and learning series / 文章与学习专题

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

## Publish an HTML article / 发布 HTML 文章

1. Copy `examples/article.html` to `content/posts/your-slug.html`. This is the single article directory, including imported visual essays. Subdirectories are supported; article filenames must be unique across the site. Published URLs remain `/blog/your-slug.html`.

   把 `examples/article.html` 复制到 `content/posts/your-slug.html`。所有文章（包括导入的交互长文）统一放在这里。支持子目录，文章文件名必须全站唯一。发布链接仍为 `/blog/your-slug.html`。

2. Set the title, article summary, card description, category, tags, and date. Use `YYYY-MM-DD` for `date`, recording the first push or sharing date. Keep it unchanged on later edits; add `updated` to record an update date.

   填写标题、正文摘要、卡片简介、分类、标签与日期。`date` 使用 `YYYY-MM-DD`，记录首次推送或分享日期；后续编辑保留该日期，如需注明更新则添加 `updated`。

3. Write both languages directly in the same HTML: English first with `data-lang="en"`, immediately followed by Chinese with `data-lang="zh"`. Use `class="parallel-text"` on paired headings or paragraphs, as in the example below. Equations, code, SVG, Canvas, and experiments can sit outside language markers to be shared.

   直接在同一个 HTML 中写两种语言：英文使用 `data-lang="en"`，紧接的中文使用 `data-lang="zh"`。配对的标题或段落使用 `class="parallel-text"`，具体写法见下方示例。公式、代码、SVG、Canvas 与实验可放在语言标记之外，由两种语言共用。

4. Build locally or push to the main branch. New articles automatically appear in the notebook, categories, timeline, search, RSS, and sitemap. Headings `h2` and `h3` receive anchors automatically; `h2` headings also populate the reading contents.

   在本地构建或推送到主分支。新文章会自动进入知识库、分类、时间线、搜索、RSS 与站点地图。`h2` 和 `h3` 自动获得锚点，`h2` 同时进入阅读目录。

Keep published slugs and URLs stable when editing a title. They identify [discussion threads](../operations/discussions.md) and [shared article links](../operations/sharing.md); a URL migration needs an explicit redirect and a review of those associations.

修改标题时保留已发布的 slug 与 URL。它们关联[讨论话题](../operations/discussions.md)和[文章分享链接](../operations/sharing.md)；迁移 URL 时应明确配置跳转，并检查这些关联。

Place complete metadata in `<script type="application/json" id="article-metadata">`. The builder also accepts `title`, `meta[name=description]`, `meta[name=category]`, and `meta[name=article:published_time]`. All metadata lives in the article itself; no separate catalog or translation files are needed. Missing required metadata, invalid dates, or a missing language stop the build.

完整元信息写在 `<script type="application/json" id="article-metadata">` 中。构建器也支持 `title`、`meta[name=description]`、`meta[name=category]` 和 `meta[name=article:published_time]`。所有元信息都放在文章自身，无需单独的目录配置或译文文件。必填元信息缺失、日期无效或缺少一种语言都会阻止构建。

Use `cardDescriptionEn` and `cardDescription` for the short English and Chinese introductions in article listings. Keep the fuller `descriptionEn` and `description` for the article opening and metadata; a missing card description falls back to that summary. Follow the [writing guide](writing.md#decide-what-the-reader-will-gain--先想清读者能得到什么) for each field's purpose and examples.

`cardDescriptionEn` 与 `cardDescription` 用于文章列表中的英文与中文短简介。正文开头和页面元信息使用较完整的 `descriptionEn` 与 `description`；未填写卡片简介时回退到正文摘要。各处文字的职责与示例见[写作指南](writing.md#decide-what-the-reader-will-gain--先想清读者能得到什么)。

```html
<h2 class="parallel-text">
  <span lang="en" data-lang="en">The idea</span>
  <span lang="zh-CN" data-lang="zh">核心思路</span>
</h2>
<p class="parallel-text">
  <span lang="en" data-lang="en">English explanation.</span>
  <span lang="zh-CN" data-lang="zh">对应的中文解释。</span>
</p>
<div data-math="y = Ax"></div>
```

| Category ID<br>分类 ID | Category<br>分类 |
| --- | --- |
| `math` | Mathematics & Algorithms<br>数学与算法 |
| `physics` | Physics & Models<br>物理与模型 |
| `systems` | Systems & Performance<br>系统与性能 |
| `benchmark` | Benchmarks<br>基准测试 |
| `biology` | Life & Neuroscience<br>生命与神经科学 |
| `travel` | Travel & Observation<br>旅行与观察 |

`draft: true` hides a draft. `archiveOnly: true` includes an article in the archive and search, but excludes it from the homepage, notebook and category listings. A numeric `featured` is a secondary ordering key for articles with the same publication date; it does not select or order the homepage's recommendations. See [recommendation behavior](#recommendation-behavior--推荐的实际行为) below. Choose `spike`, `rust`, `gpu`, `matrix`, `math`, `wave`, `benchmark`, or `vision` for `art` to generate local SVG illustrations.

`draft: true` 隐藏草稿。`archiveOnly: true` 让文章进入归档与搜索，不进入首页、知识库和分类列表。数字形式的 `featured` 是同一发布日期文章的次级排序依据，不决定首页推荐的入选与顺序。具体见下方[推荐规则](#recommendation-behavior--推荐的实际行为)。`art` 可选 `spike`、`rust`、`gpu`、`matrix`、`math`、`wave`、`benchmark` 或 `vision`，用于生成本地 SVG 配图。

## Learning series / 学习专题

To add an article to a learning series, include an optional `series` object in its article metadata. Use the same `id`, `titleEn`, and `title` across installments and a unique positive `part` number. The build generates `/series/<id>/`, lists it in the notebook, adds it to the sitemap, and connects published installments with previous/next navigation. The homepage does not automatically list every new series; review its reading routes separately in `scripts/templates.mjs`. Drafts are excluded; planned articles should remain plain text until published. Duplicate part numbers and inconsistent series titles fail the build.

若要把文章加入学习专题，在元信息中填写可选的 `series` 对象。各期使用相同的 `id`、`titleEn` 和 `title`，并填写不重复的正整数 `part`。构建自动生成 `/series/<id>/`，在知识库添加入口、写入站点地图，并为已发布文章生成前后期导航。首页不会自动列出所有新增专题，需要另行核对 `scripts/templates.mjs` 中的阅读路径。草稿不计入；未发布规划应保留为普通文字。期数重复或专题名称不一致会阻止构建。

```json
{
  "series": {
    "id": "pytorch-internals",
    "titleEn": "Inside PyTorch",
    "title": "PyTorch 源码之旅",
    "part": 2
  }
}
```

Every installment automatically includes the same compact, collapsible series contents after its body, before related reading, with publication dates and the current installment highlighted. No preview markup is required. To list future parts, declare `series.roadmap` once in any installment’s metadata, for example `[{"part":3,"titleEn":"Operator dispatch","title":"算子调度"}]`. All installments share that plan; published metadata takes precedence as each part appears. New published parts are always added automatically.

每一期正文之后、相关阅读之前都会自动出现同一份紧凑、可折叠的专题目录，显示发布日期并突出当前期，无需添加预览标记。若需列出后续规划，在任一期元信息的 `series.roadmap` 中声明一次即可，例如 `[{"part":3,"titleEn":"Operator dispatch","title":"算子调度"}]`。各期共用这份规划；对应章节发布后，自动采用实际文章标题与日期。新增的已发布章节始终自动加入。

An optional `<div data-series-preview>` can embed the list elsewhere; set its attribute to a series ID to reference another series. Child entries with `data-series-part` and bilingual titles can add local planned topics. Avoid duplicating the automatically generated directory in ordinary installments.

如需在其他位置引用，可选用 `<div data-series-preview>`，属性值可指定其他专题 ID。带 `data-series-part` 与双语标题的子元素可补充局部规划。普通专题文章无需再重复插入自动生成的目录。

The series landing page can add an editorial reading guide from `content/editorial.json`. Each entry in `series` names an existing series `id`, gives bilingual `audience`, `prerequisites`, and `outcomes`, and lists short `paths` with a stable `id`, bilingual `title` and `description`, and published `parts` numbers. These guides organize already published material; they must not imply an unimplemented course, tool, or credential.

专题首页可从 `content/editorial.json` 补充阅读指引。`series` 中的条目指定已有专题 `id`，以双语填写 `audience`（适合谁）、`prerequisites`（先修知识）和 `outcomes`（读后能理解什么）；`paths` 使用稳定的 `id`、双语 `title` 和 `description`，以及已发布的 `parts` 期数。指引用于组织现有内容，不应暗示尚未完成的课程、工具或作者资历。

## Site copy and reading entry points / 站点文案与阅读入口

`content/editorial.json` contains bilingual homepage mission and about-page copy, series reading guides and curated article connections. `scripts/editorial-views.mjs` renders the guides and connections; `scripts/templates.mjs` composes the homepage and about page. README states the broader mission; `scripts/config.mjs` and `scripts/seo.mjs` supply site and page descriptions. When wording changes, review these surfaces for consistency without copying the full about page into each one. Content choices follow the [writing guide](writing.md).

`content/editorial.json` 保存双语首页使命、关于页文案、专题阅读指引和人工文章联系；`scripts/editorial-views.mjs` 渲染指引与联系，`scripts/templates.mjs` 组织首页和关于页。README 说明整体使命，`scripts/config.mjs` 与 `scripts/seo.mjs` 提供站点和页面描述。修改文案时检查这些位置是否一致，无需到处复制整段关于页。内容取舍遵循[写作指南](writing.md)。

### Recommendation behavior / 推荐的实际行为

These lists are generated at build time. They do not vary by visitor or require analytics.

以下列表在构建时生成，不按访客个人变化，也不依赖统计服务。

| Entry point / 入口 | Selection / 取文逻辑 | Source / 修改位置 |
| --- | --- | --- |
| Homepage “Three places to begin” / 首页“从这里读起” | Three fixed slugs: PyTorch part 1, the visual system and the Chaoshan essay. Missing or archive-only entries are omitted; no automatic replacement.<br>固定选取 PyTorch 第一期、人类视觉与潮汕游记；条目缺失或仅归档时省略，不自动补位。 | `home()` → `choices`, `scripts/templates.mjs` |
| Homepage “Recent notes” / 首页“最近更新” | Up to four visible articles ordered by publication `date`, excluding the fixed choices. An `updated` value does not move an article to the front.<br>按发布 `date` 排序，排除固定精选后取最多四篇；`updated` 不改变排序。 | `scripts/content.mjs`, `home()` in `scripts/templates.mjs` |
| About-page reading links / 关于页推荐 | Fixed links to the PyTorch series, visual-system article and Chaoshan essay. New articles do not replace them.<br>固定链接到 PyTorch 专题、人类视觉与潮汕游记，新文章不会自动替换。 | `about()` in `scripts/templates.mjs` |
| Article “Continue exploring” / 文末“继续探索” | Curated directed links first, then relevant automatic matches, up to three distinct articles.<br>人工指定的有向联系优先，再用自动相关结果补充，最多三篇且不重复。 | `content/editorial.json` → `related`; [related-reading rules / 相关阅读规则](#related-reading--相关阅读) |
| Daily-math further reading / 每日数学延伸阅读 | Optional explicit links for that topic, preserved on its dated page.<br>按知识点填写可选的明确链接，日期页面保留对应关系。 | `content/daily-math/topics.json` → `related`; [daily mathematics / 每日数学](daily-mathematics.md) |

### Related reading / 相关阅读

For a deliberate connection between articles, add a directed entry to `content/editorial.json` → `related`: `source` and `target` are article slugs; `reason` contains `en` and `zh`. Explain the conceptual step, not a generic “read more”. Keep each connection relevant to the source article and avoid self-links or duplicate targets.

要指定明确的跨文联系，在 `content/editorial.json` 的 `related` 数组中添加有方向的条目：`source` 和 `target` 填文章 slug，`reason` 含 `en` 与 `zh`。理由要说明知识如何衔接，避免泛泛地写“了解更多”。每条联系应贴合原文章，不指向自身、不重复目标。

Automatic matches use bilingual titles, tags, summaries, headings and body text to calculate TF-IDF cosine similarity; shared tags and nearby installments in the same series receive extra weight. Category is a small secondary signal, and publication date does not determine relevance. Code blocks, controls and roadmap links are excluded from body text. Weak matches are omitted rather than used to fill every slot. Add a bilingual HTML article as usual and rebuild to update recommendations.

自动匹配使用中英文标题、标签、摘要、章节标题与正文计算 TF-IDF 余弦相似度，共同标签和同专题相邻期数额外加权。分类只作为较弱的辅助信号，发布日期不决定相关性。正文分析排除代码块、控件和专题预览链接。相关性不足时不强行补满；照常新增双语 HTML 并构建即可更新推荐。

## Evidence and observation / 证据与观察

Apply the [writing guide](writing.md) to articles and site copy, including its standards for [technical evidence](writing.md#explain-mechanisms-without-skipping-the-hard-step--不跳过最难懂的一步) and [photographs, recollections and corrected travel locations](writing.md#write-observation-before-declaring-emotion--先写观察再让情绪发生).

文章与站点文案均遵循[写作指南](writing.md)，包括[技术论述的证据](writing.md#explain-mechanisms-without-skipping-the-hard-step--不跳过最难懂的一步)，以及[照片、回忆与游记地点修正](writing.md#write-observation-before-declaring-emotion--先写观察再让情绪发生)的要求。

Start from the [article template](../../examples/article.html); follow the [language rules](language.md) and [presentation standards](presentation.md).

从[文章模板](../../examples/article.html)开始，按[双语规范](language.md)和[排版规范](presentation.md)编写。
