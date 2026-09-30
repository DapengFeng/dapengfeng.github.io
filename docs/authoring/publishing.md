# Articles and learning series / 文章与学习专题

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

## Publish an HTML article / 发布 HTML 文章

1. Copy `examples/article.html` to `content/posts/your-slug.html`. This is the single article directory, including imported visual essays. Subdirectories are supported; article filenames must be unique across the site. Published URLs remain `/blog/your-slug.html`.

   把 `examples/article.html` 复制到 `content/posts/your-slug.html`。所有文章（包括导入的交互长文）统一放在这里。支持子目录，文章文件名必须全站唯一。发布链接仍为 `/blog/your-slug.html`。

2. Set the title, summary, category, tags, and date. Use `YYYY-MM-DD` for `date`, recording the first push or sharing date. Keep it unchanged on later edits; add `updated` to record an update date.

   填写标题、摘要、分类、标签与日期。`date` 使用 `YYYY-MM-DD`，记录首次推送或分享日期；后续编辑保留该日期，如需注明更新则添加 `updated`。

3. Write both languages directly in the same HTML: English first with `data-lang="en"`, immediately followed by Chinese with `data-lang="zh"`. Use `class="parallel-text"` on paired headings or paragraphs, as in the example below. Equations, code, SVG, Canvas, and experiments can sit outside language markers to be shared.

   直接在同一个 HTML 中写两种语言：英文使用 `data-lang="en"`，紧接的中文使用 `data-lang="zh"`。配对的标题或段落使用 `class="parallel-text"`，具体写法见下方示例。公式、代码、SVG、Canvas 与实验可放在语言标记之外，由两种语言共用。

4. Build locally or push to the main branch. New articles automatically appear in the notebook, categories, timeline, search, RSS, and sitemap. Headings `h2` and `h3` receive anchors automatically; `h2` headings also populate the reading contents.

   在本地构建或推送到主分支。新文章会自动进入知识库、分类、时间线、搜索、RSS 与站点地图。`h2` 和 `h3` 自动获得锚点，`h2` 同时进入阅读目录。

Place complete metadata in `<script type="application/json" id="article-metadata">`. The builder also accepts `title`, `meta[name=description]`, `meta[name=category]`, and `meta[name=article:published_time]`. All metadata lives in the article itself; no separate catalog or translation files are needed. Missing required metadata, invalid dates, or a missing language stop the build.

完整元信息写在 `<script type="application/json" id="article-metadata">` 中。构建器也支持 `title`、`meta[name=description]`、`meta[name=category]` 和 `meta[name=article:published_time]`。所有元信息都放在文章自身，无需单独的目录配置或译文文件。必填元信息缺失、日期无效或缺少一种语言都会阻止构建。

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

`draft: true` hides a draft. `archiveOnly: true` includes an article in the archive and search, but excludes it from the homepage notebook. A numeric `featured` sets its order among featured notes. Choose `spike`, `rust`, `gpu`, `matrix`, `math`, `wave`, `benchmark`, or `vision` for `art` to generate local SVG illustrations.

`draft: true` 隐藏草稿。`archiveOnly: true` 让文章仅进入归档与搜索，不进入首页知识库。数字形式的 `featured` 指定首页精选顺序。`art` 可选 `spike`、`rust`、`gpu`、`matrix`、`math`、`wave`、`benchmark` 或 `vision`，用于生成本地 SVG 配图。

## Learning series / 学习专题

To add an article to a learning series, include an optional `series` object in its article metadata. Use the same `id`, `titleEn`, and `title` across installments and a unique positive `part` number. The build generates `/series/<id>/`, links it from the homepage and notebook, adds it to the sitemap, and connects published installments with previous/next navigation. Drafts are excluded; planned articles should remain plain text until published. Duplicate part numbers and inconsistent series titles fail the build.

若要把文章加入学习专题，在元信息中填写可选的 `series` 对象。各期使用相同的 `id`、`titleEn` 和 `title`，并填写不重复的正整数 `part`。构建自动生成 `/series/<id>/`，在首页和知识库添加入口、写入站点地图，并为已发布文章生成前后期导航。草稿不计入；未发布规划应保留为普通文字。期数重复或专题名称不一致会阻止构建。

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

## Related reading / 相关阅读

“Continue exploring” is recomputed during each build from published articles. Bilingual titles, tags, summaries, headings and body text contribute to TF-IDF cosine similarity; shared tags and nearby installments in the same series receive extra weight. Category is a small secondary signal, and publication date does not determine relevance. Code blocks, controls and roadmap links are excluded from body text. Up to three related articles are shown; weak matches are omitted. No external service or visitor tracking is required. Add a bilingual HTML article as usual and rebuild to update recommendations.

“继续探索”在每次构建时根据已发布文章重新计算。中英文标题、标签、摘要、章节标题与正文参与 TF-IDF 余弦相似度计算，共同标签和同专题相邻期数额外加权。分类只作为较弱的辅助信号，发布日期不决定相关性。正文分析排除代码块、控件和专题预览链接。最多显示三篇，相关性不足时不强行补满。不依赖外部服务或访客追踪；照常新增双语 HTML 并构建即可更新推荐。

Start from the [article template](../../examples/article.html); follow the [language rules](language.md) and [presentation standards](presentation.md).

从[文章模板](../../examples/article.html)开始，按[双语规范](language.md)和[排版规范](presentation.md)编写。
