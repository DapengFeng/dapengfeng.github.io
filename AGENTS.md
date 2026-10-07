# Agent entry point / Agent 工作入口

Start with the [mission](README.md) and use the [documentation index](docs/README.md) to read only the guides relevant to the task. Each rule has one owning guide; link to it instead of copying it. A later explicit user instruction takes precedence. Distinguish current behavior from a proposal and update the owning guide when behavior changes.

先了解[网站使命](README.md)，再从[文档索引](docs/README.md)按任务阅读相关指南。每项规则只在一份归属文档维护，其他位置使用链接，不复制正文。用户后续的明确指示优先；区分现有行为与建议，行为变化时同步更新归属文档。

| Task / 任务 | Guide / 规范归属 |
| --- | --- |
| Mission, content structure, prose and evidence / 使命、内容组织、文笔与依据 | [Writing / 内容与写作](docs/authoring/writing.md) |
| Page layout, diagrams, controls, equations and code / 页面、图示、控件、公式与代码 | [Presentation / 设计与排版](docs/authoring/presentation.md) |
| Article metadata, dates, series and recommendation selection / 文章元信息、日期、专题与推荐取文 | [Publishing / 发布与推荐](docs/authoring/publishing.md) |
| Language behavior and paired markup / 语言行为与双语标记 | [Language / 双语规范](docs/authoring/language.md) |
| Daily topics, scheduling and background playback / 每日主题、排期与背景播放 | [Daily mathematics / 每日数学](docs/authoring/daily-mathematics.md) |
| Sources, preview, loading and checks / 源码、预览、加载与检查 | [Local workflow / 本地工作流](docs/development/local-workflow.md), [testing / 测试](docs/development/testing.md) |
| Deployment, search, analytics, discussions, sharing and support / 部署、搜索、统计、讨论、分享与赞赏 | [Service guides / 服务文档](docs/README.md#3-deployment-and-services--部署与服务) |

Inspect `git status` and the relevant diff before editing; preserve existing work. Edit sources, not generated `dist/`. Commit or push when requested. Keep ignored private research, local environment settings and account details out of public documentation and commits.

修改前检查 `git status` 和相关差异，保留已有工作。修改源码，不编辑生成的 `dist/`；按用户请求提交或推送。已忽略的私人研究、本地环境配置和账号资料不进入公开文档或提交。

Verify the current preview port and build state. A local preview, a commit, a push and a Pages deployment are different outcomes. Run checks appropriate to the change and inspect screenshots for visual work. Report the actual changes, checks and unresolved limits; do not claim that mocks prove a live publication, login or payment, or perform such actions as a test side effect.

核对当前预览端口与构建状态。本地预览、提交、推送和 Pages 部署是不同结果。按影响范围执行检查，视觉修改实际查看截图。交代真实修改、检查结果与剩余限制，不把模拟检查说成真实发布、登录或付款，也不把这些操作作为测试副作用。
