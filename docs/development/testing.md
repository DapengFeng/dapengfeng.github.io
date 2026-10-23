# Testing and validation / 测试与验证

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

## Site checks / 站点检查

```bash
npm run build
npm test
npx playwright install chromium webkit
npm run test:browsers
```

Browser checks use the browser versions pinned by Playwright, matching CI by default. Install Chromium and WebKit with `npx playwright install chromium webkit`; on Linux, add `--with-deps` when system libraries are missing. `PLAYWRIGHT_CHROMIUM_EXECUTABLE` can select an existing Chromium executable. Run `npm run test:browsers` for interaction, language, sorting, layout, formula, and accessibility checks. Every configured suite runs even if an earlier one fails; each has a five-minute limit. Per-suite logs and the result summary are saved under `test-results/browser/`; any failure makes the command exit nonzero. The compiler tests use isolated API responses; add `--live` to `node tests/compiler.mjs` to verify Godbolt itself. Static output is written to `dist/` and requires no server-side application.

浏览器检查默认使用 Playwright 锁定的浏览器版本，与 CI 保持一致。执行 `npx playwright install chromium webkit` 安装 Chromium 与 WebKit；Linux 缺少系统库时增加 `--with-deps`。可通过 `PLAYWRIGHT_CHROMIUM_EXECUTABLE` 指定已安装的 Chromium。`npm run test:browsers` 检查交互、语言、排序、排版、公式和无障碍。前面的测试失败后仍会继续检查后续项目；每项限时五分钟，逐项日志和结果汇总保存在 `test-results/browser/`，任意一项失败都会使命令返回非零状态。编译器测试使用隔离的接口响应；运行 `node tests/compiler.mjs --live` 可验证 Godbolt 本身。静态产物输出到 `dist/`，无需服务器端程序。

`npm run test:cross-browser` runs a small shared-UI smoke suite in Chromium and WebKit: desktop/mobile layouts, navigation, native dialogs, sharing, discussions, and the visual-system figure fallback. It saves fixed-view screenshots of the homepage, technical article and travel essay at 1440px and 390px widths under `test-results/visual/`. These screenshots support visual review; the suite does not compare pixels with a baseline or certify Safari on real devices.

`npm run test:cross-browser` 在 Chromium 与 WebKit 中执行少量共用界面冒烟检查，覆盖桌面／手机排版、导航、原生对话框、分享、讨论以及视觉系统图示的回退显示。首页、技术文章与游记在 1440px 和 390px 宽度下的固定视角截图保存在 `test-results/visual/`。截图用于人工复核，不进行像素基线比对，也不等于真实设备上的 Safari 验证。

Catalog checks derive category membership from the published search index and compare article URLs, so new posts do not require changing fixed totals. Search checks target known relevant and unrelated articles rather than assuming there is exactly one result. For asynchronous search, wait for the current query's expected state; an existing result from the previous query does not prove the update has finished. Keep fixed counts only for deliberate UI limits or article-specific fixtures.

文章清单检查从已构建的搜索索引获取分类成员，并核对文章链接，新增文章无需修改写死的总数。搜索检查针对已知相关与无关文章，不假定结果恰好只有一条。异步搜索应等待当前关键词对应的预期状态，上一轮仍显示的结果不能证明更新已完成。仅对明确的界面数量限制或特定文章的固定结构保留数量断言。

## Recorded performance experiments / 已记录的性能实验

The dense/sparse parallel-library experiment has its own [reproduction and measurement guide](../../benchmarks/parallel-libraries/README.md). Its C++, Rust and Python workloads are run manually, with independent correctness checks, separate timing and process-tree PSS passes, and a bounded sustained-load profile. Website builds publish the recorded artifacts; they do not run stress tests or enforce performance rankings on CI hardware.

稠密／稀疏并行库实验的命令与测量边界见[复现说明](../../benchmarks/parallel-libraries/README.md)。C++、Rust 与 Python 负载手动运行，包含独立正确性校验、分开的计时与进程树 PSS 测量，以及有边界的持续压测。网站构建只发布记录好的资源，不在 CI 机器上运行压测或按性能排名设置门槛。

## Open-model deployment economics / 开源模型部署经济性

The investment article uses `/blog/llm-deployment-roi.html` and has no fixed provider example. Its former DeepSeek URL redirects to this page. The currency-neutral [finance engine](../../content/assets/llm-economics/finance.js) evaluates explicit monthly cash flows; the [token model](../../content/assets/llm-economics/token-model.js) creates those flows from billing and measured capacity. The public calculator starts with neutral inputs. [Charts](../../content/assets/llm-economics/charts.js) recompute price, first-month cost, cumulative cash and demand sensitivity from the same parameters. The dated legacy example remains a regression fixture, not the article's default or a performance claim.

投资文章使用 `/blog/llm-deployment-roi.html`，不再采用固定供应商案例；原 DeepSeek 链接跳转到新页面。[通用财务引擎](../../content/assets/llm-economics/finance.js)接收同币种的逐月现金流，[token 模型](../../content/assets/llm-economics/token-model.js)由计费与实测容量生成现金流。公开计算器从中性参数开始，[动态图](../../content/assets/llm-economics/charts.js)以同一组参数重算单价、首月成本、累计现金与需求敏感性。带日期的旧案例仅保留为回归样本，不是文章默认值或性能保证。

The [token comparisons](../../content/assets/llm-economics/token-examples.js) hold total tokens fixed across 90/10, 50/50 and 10/90 input/output mixes, with 0% and 80% **input** cache hits. The live view uses the entered monthly total (or an explicitly labeled 100M reference), plus a 10M/100M/1,000M uncached volume matrix. All bills reuse `billTokens`; missing cache prices leave the 80% case unavailable. These are API-only controlled comparisons, not equal-task or hardware-capacity claims. The static example uses fictional CNY rates and remains readable without JavaScript; regenerate it with `node scripts/render-token-examples.mjs` and verify it with `--check`. The pure calculation and view tests cover conservation, weighted rates, shared chart scales, missing prices, escaping and export metadata.

[Token 对照](../../content/assets/llm-economics/token-examples.js)固定总 token 数，比较 90/10、50/50、10/90 输入／输出配比及 0%、80% **输入**缓存命中率。实时图采用所填月总量，未填时明确标出 100M 参考量；另有 10M／100M／1,000M 无缓存用量矩阵。账单统一复用 `billTokens`，缺缓存报价时不计算 80% 场景。这些只是 API 计费对照，不代表相同任务或硬件容量。静态算例采用虚构人民币单价，禁用 JavaScript 仍可阅读；用 `node scripts/render-token-examples.mjs` 重建，加 `--check` 校验。纯计算及视图测试覆盖守恒、加权单价、共用图表尺度、缺价、转义与导出来源信息。

`node --test tests/deepseek-roi.test.mjs` checks generic cash flows, replacement outlays, negative terminal values, independent token quantities, bottlenecks, fallback, caches, growth, conservation, investment ceilings and reversals. ROI is net cash gain divided by initial outlay, not an annualized return; payback excludes terminal proceeds. The archived calculation can still be checked with `node scripts/calculate-deepseek-roi.mjs --check`. Regenerate the generic cover with `node scripts/render-llm-economics-cover.mjs`.

`node --test tests/deepseek-roi.test.mjs` 检查通用现金流、换机支出、负终值、独立 token 用量、瓶颈、回退、缓存、增长、守恒、投资上限与反转。ROI 是净现金收益除以初投，不是年化收益；回本排除终值。旧计算样本仍可用 `node scripts/calculate-deepseek-roi.mjs --check` 校验；通用封面通过 `node scripts/render-llm-economics-cover.mjs` 重建。

`npm run test:economics` checks provider filtering, model search, imported prices and provenance, million-token and original-currency inputs, edits, currency changes, invalid inputs, reset, downloaded cash flows and chart updates on the built article. It covers desktop/mobile widths and language modes, plus catalog failure and no-JavaScript reading. Run parser/fallback tests and follow the source-refresh rules in [model pricing](../operations/model-pricing.md). Daily publications include the economics browser check because prices can change without source changes. Tests validate calculation and UI behavior, not real throughput, service acceptance or a live deployment.

`npm run test:economics` 检查构建后的供应商筛选、模型搜索、价格及来源、百万 token 与原币金额输入、修改、币种切换、非法输入、重置、下载现金流与图表更新，覆盖桌面／手机、语言模式、目录失败与禁用 JavaScript 阅读。解析／回退测试及刷新规则见[模型报价](../operations/model-pricing.md)。每日发布包含经济性浏览器检查，因为价格可能在源码不变时更新。测试验证运算与界面行为，不证明真实吞吐、服务验收或线上部署。

Inspect chart screenshots, not just bounding boxes. Labels and legends must remain legible at 320px and in English, Chinese and bilingual modes. Changing a model clears stale capacity; invalid inputs clear stale results. Rate plots can display without a workload, while cash plots need demand. Cumulative plots distinguish operating recovery from terminal resale, and sensitivity holds costs and capacity fixed. Monetary units must match every input and export.

实际查看图表截图，不能只看边界框。320px 下及中、英、双语模式的标签与图例都须清晰。换模型清除旧容量，非法输入清除旧结果；无用量可看单价，有需求才展示现金图。累计图区分运营回本与期末出售，敏感性固定成本和容量。输入与导出的全部金额须统一币种。

## Visual acceptance / 视觉验收

The DGX Spark / Mac Studio article's memory-budget and bandwidth figures are generated by [scripts/render-local-llm-budget.py](../../scripts/render-local-llm-budget.py), using Python, Matplotlib and Noto Sans CJK SC. Run `python3 scripts/render-local-llm-budget.py` to regenerate the bilingual desktop/mobile SVGs, cover and [calculation data](../../content/assets/local-llm/data.json). Its arithmetic checks validate a stated 70B memory example and conditional weight-streaming ceilings, not hardware benchmarks. The deployment commands are checked against upstream documentation, not executed on the target machines. Review current product generations, prices, availability, runtime versions and model support before updating the article; inspect both language editions, narrow-screen figures, stacked comparison tables and code blocks after rebuilding.

DGX Spark／Mac Studio 文章的内存预算与带宽图由 [scripts/render-local-llm-budget.py](../../scripts/render-local-llm-budget.py) 生成，依赖 Python、Matplotlib 与 Noto Sans CJK SC。运行 `python3 scripts/render-local-llm-budget.py` 可重建中英文桌面／手机 SVG、封面与[计算数据](../../content/assets/local-llm/data.json)。算术检查验证指定的 70B 内存算例及条件性权重读取上界，不是硬件跑分；部署命令核对上游文档，未在目标机器执行。更新文章前重新核对产品代际、价格、供货、运行时版本和模型支持；构建后查看两种语言、窄屏图、纵向展开的比较表与代码块。

The Roadster article's launch-threshold and fan-power figures are generated by [scripts/render-roadster-models.py](../../scripts/render-roadster-models.py). Run it manually with Python, NumPy and Matplotlib to check the numerical balances and regenerate the SVGs. The script records the illustrative assumptions; its checks validate those calculations, not Tesla performance. Website builds use the checked-in SVGs and do not install plotting dependencies.

Roadster 文章的一秒起步门槛图与风扇功耗图由 [scripts/render-roadster-models.py](../../scripts/render-roadster-models.py) 生成。使用装有 NumPy 和 Matplotlib 的 Python 手动运行，核对数值平衡并重新生成 SVG。脚本记录示例假设，校验验证的是模型计算，并非 Tesla 性能。网站构建使用已提交的 SVG，不安装绘图依赖。

Part 2's envelopes and decision-boundary figures are generated by [scripts/render-roadster-power.py](../../scripts/render-roadster-power.py), using the same Python plotting dependencies. It writes English/Chinese SVGs, narrow-screen versions of the target and boundary plots, a raster cover and [model.json](../../content/assets/roadster-power/model.json). The data separates the constant-acceleration target, equal-cell wiring comparison, controlled upgrades, gear sensitivity, downforce cost allowances and variable-acceleration one-second requirements. Checks cover power conservation, transition continuity, gear speed/torque relations and model speed bounds; analytic launch times and boundary points are checked against independent numerical integration. The original downforce comparison omits device costs; the cost map restores added mass and auxiliary draw from the shared DC budget while holding downforce fixed. It does not establish a realizable fan design. The one-second map varies both force and power capacities rather than retaining the baseline motor and grip limits. Fixed efficiencies and static envelopes remain illustrative; validation does not establish real vehicle performance or a battery safety rating.

第二期的包络图与决策边界图由 [scripts/render-roadster-power.py](../../scripts/render-roadster-power.py) 使用相同 Python 绘图依赖生成，输出中英文 SVG、目标与边界图的窄屏版本、位图封面和 [model.json](../../content/assets/roadster-power/model.json)。数据分别记录匀加速目标、相同电芯接线对照、受控升级、齿比敏感性、下压力代价预算，以及变加速的一秒需求。检查涵盖功率守恒、转折连续性、传动比对应的转速／转矩关系和模型转速边界，解析用时与边界点以独立数值积分核对。原始下压力对照省略装置代价，代价图则在固定下压力的前提下补回新增质量，并从共用直流预算扣除辅助用电；它不证明风扇设计可实现。一秒可行域同时改变力与功率能力，不沿用基线电机与抓地上限。固定效率与静态包络仍是示例假设；校验不证明真实车辆性能或电池安全额定值。

Part 3's tyre-load and force-allocation figures are generated by [scripts/render-roadster-contact.py](../../scripts/render-roadster-contact.py), with assumptions and calculated cases recorded in [model.json](../../content/assets/roadster-contact/model.json). The renderer checks the balance of all four normal loads, the concavity of the illustrative tyre load-sensitivity law and a self-consistent steady-cornering limit set by the first saturated axle with lateral load transfer. It also checks a specified minimum wheel-load reserve, sensitivity to the tyre exponent and controlled changes to balance, centre-of-mass height and transfer distribution. The four-wheel braking-and-turning calculation checks the total longitudinal/lateral forces, net yaw moment and individual tyre utilisations, then locates the first constraint reached as braking demand or the downforce state changes. Recompute loads and the stated allocation at each candidate state; the boundaries belong to that allocation rule, not a globally optimal controller.

第三期的轮胎载荷与受力分配图由 [scripts/render-roadster-contact.py](../../scripts/render-roadster-contact.py) 生成，假设和算例记录在 [model.json](../../content/assets/roadster-contact/model.json)。脚本检查四个车轮的法向载荷平衡、示例轮胎载荷敏感性规律的凹性，以及包含横向载荷转移、由首先饱和车轴决定的自洽稳态转弯极限；同时检查指定最小轮载余量、轮胎指数敏感性，以及对下压力分配、质心高度和载荷转移分配的受控改动。四轮制动转弯计算核对纵横向总力、净横摆力矩与每轮利用率，再定位制动需求或下压力状态改变时首先到达的约束。每个候选状态都重新计算轮载并应用指定分配；边界属于该规则，不代表全局最优控制器。

These are quasi-static, idealized comparisons, not Tesla measurements, transient simulations or full predictions of vehicle handling. The 500 N reserve is an illustrative constraint rather than a certified safety threshold; tyre-exponent scenarios are not a measured confidence interval. Pressure-loss comparisons do not establish actuator response or recovery time. The article's measurement plan remains proposed work. Regenerate with the same Python plotting dependencies and inspect both language editions, narrow-screen figures, locally scrollable equations and the validation table. At narrow widths, the table presents each check vertically with its condition and decision labels; verify that all three columns remain readable without horizontal scrolling. Figure captions and assumption notes use the article's quieter 16px type, below the 18px body copy. Part 3 keeps one reading column below 1024px, with a collapsed native contents panel when JavaScript is available; test both sides of that breakpoint, keyboard expansion and the expanded no-JavaScript fallback. The responsive picture sources reserve their own aspect ratios; delay image responses to confirm the surrounding copy stays in place.

这些是准静态理想化对照，不是 Tesla 实测、瞬态仿真或完整的操稳预测。500 N 余量是示例约束，不是经过认证的安全门槛；轮胎指数对照也不是实测置信区间。失压前后的状态对照不能证明执行器响应或恢复时间。文章中的测量计划仍是后续工作。使用相同 Python 绘图依赖重新生成，检查两种语言、窄屏图、可局部滚动的公式与验证计划表。窄屏时，表格按验证项纵向展开，并显示条件与判断标签；确认三列内容均可直接阅读，无需横向滚动。图注与假设说明使用较安静的 16px 字号，低于 18px 正文层级。第三期在 1024px 以下保持单栏阅读，启用 JavaScript 时原生目录面板默认收起；检查断点两侧、键盘展开及无 JavaScript 时保持展开的回退。响应式图片来源分别预留自身宽高比；延迟图片响应，确认前后正文不会因图片加载而移动。

The fourth Roadster installment uses [scripts/render-roadster-cold-gas.py](../../scripts/render-roadster-cold-gas.py) to generate bilingual desktop/mobile thrust, tank-state and storage-volume figures, plus [calculation data](../../content/assets/roadster-cold-gas/model.json). Run it with Python, NumPy, Matplotlib and Noto Sans CJK JP. The renderer checks its closed-form blowdown histories against independent mass/internal-energy integration, integrates discharged mass, enthalpy and impulse, and checks the ideal refill-work integral. These checks validate the stated ideal-gas calculation, not Tesla specifications, actual high-pressure air, regulators or vehicle performance. Both architectures use the same tank and initial thrust but different nozzle areas; the regulated trace stops at loss of regulation, not at zero thrust. The refill calculation holds the tank at 300 K throughout filling. Review the in-page flow schematic, both figure languages, responsive image aspect ratios, long formulas and stacked validation table at 320px through desktop widths. Part 4 uses the same tablet-width reading-column and native contents behavior as Part 3; keep no-JavaScript reading intact.

Roadster 第四期使用 [scripts/render-roadster-cold-gas.py](../../scripts/render-roadster-cold-gas.py) 生成双语桌面／手机推力、罐内状态与储气容积图，以及[计算数据](../../content/assets/roadster-cold-gas/model.json)。运行依赖 Python、NumPy、Matplotlib 与 Noto Sans CJK JP。脚本用独立的质量／内能积分核对放气闭式解，再核对排气质量、焓与冲量积分，以及理想补气功积分。这些检查验证所声明的理想气体计算，不验证 Tesla 规格、实际高压空气、稳压器或整车性能。两种供气方式保持相同储罐和初始推力，但喷嘴面积不同；稳压曲线结束代表无法继续稳压，不代表推力归零。补气计算在整个充气过程中保持罐温 300 K。检查页面流路示意、两种语言图示、响应式图片宽高比、长公式与纵向展开的验证表，覆盖 320px 至桌面宽度。第四期沿用第三期的平板单栏与原生目录行为，并保留无 JavaScript 阅读。

Part 5 uses [scripts/render-roadster-control.py](../../scripts/render-roadster-control.py) with Python, NumPy, Matplotlib and Noto Sans CJK JP to reproduce the bilingual force-budget, attainable-set and delay figures, the geometry schematic, cover and [calculation data](../../content/assets/roadster-control/model.json). Run `python3 scripts/render-roadster-control.py`. Its frozen planar example allows independently selected signed left/right tyre forces; actual powertrain topology, tyre dynamics and lateral-force changes are outside that relaxed set. Checks compare force ceilings with transformed-box edge intersections, inverse allocations with the region inequality, and angular responses with independent integration; they also check saturation before reallocation and directional reserve. Delay plots integrate one uncompensated moment, not a closed-loop recovery or safety result. Inspect the geometry, diagram language, mobile image aspect ratios, formulas and native contents at 320px through desktop widths. Website builds consume the checked-in assets without plotting dependencies.

第五期通过 [scripts/render-roadster-control.py](../../scripts/render-roadster-control.py)，使用 Python、NumPy、Matplotlib 与 Noto Sans CJK JP，复现双语受力预算、可行域、延迟图，以及几何示意、封面和[计算数据](../../content/assets/roadster-control/model.json)。运行 `python3 scripts/render-roadster-control.py`。冻结的平面算例允许独立选择有符号的左右轮胎力；真实动力拓扑、轮胎动态与横向力变化不包含在这一放宽集合内。检查将受力上限与变换后边界求交对比，将逆分配与区域不等式对比，再以独立积分核对角运动，并检查重分配之前的物理饱和与方向性余量。延迟图只积分单独的未抵消力矩，不代表闭环恢复或安全验证。检查 320px 至桌面宽度下的几何示意、图示语言、手机图宽高比、公式与原生目录。网站构建直接使用已提交资产，不安装绘图依赖。

Open and inspect the screenshots for a visual change; a passing overflow assertion is not sufficient. When shared styles or templates change, cover the homepage, a technical article, a travel essay and the about page in desktop and mobile layouts. Check English, Chinese and bilingual modes, and include a short viewport when the first screen is part of the requirement. For a local change, focus this review on the affected component and its shared neighbors.

视觉改动后要实际打开截图查看，通过溢出断言并不足够。共用样式或模板变化时，覆盖首页、技术文章、游记和关于页的桌面与手机布局；检查英文、中文和双语模式，涉及首屏要求时增加矮屏。局部修改则重点复核受影响组件及相邻的共用界面。

Look for orphaned title characters, long English strings, oversized decorative marks, distorted or unrelated photo crops, detached captions, clipped formulas and overlapping floating controls. Inherited CSS can be the cause: for example, a travel table of contents should explicitly reset a shared multi-column rule when changing to horizontal layout. Do not repair one viewport with unscoped selectors or hard-coded line breaks that damage another language.

注意标题孤字、英文长词撑宽、过大的装饰图形、图片变形或裁掉主体、图注与图片脱节、公式裁切、浮动控件遮挡。原因可能来自 CSS 继承：例如游记目录改成横向布局时，要明确重置共用的多栏规则。不使用无作用域的选择器或破坏其他语言的硬换行，只修好一个尺寸。

The cross-browser suite freezes a known published date and advances a controlled clock for stable screenshots. Motion behavior is tested separately by `npm run test:daily-math`, including visible progression, the 30 fps limit, ignored legacy static settings, offscreen/hidden-tab pause and resume, midnight changes and fixed dated entries. Stabilize tests with test tools; never reintroduce a production display setting for screenshots.

跨浏览器套件固定一个已发布日期，以受控时钟推进画面，保证截图稳定。播放行为单独由 `npm run test:daily-math` 验证，包括可见时变化、30 fps 上限、忽略旧静态设置、离屏与隐藏标签页的暂停恢复、跨日切换及固定日期条目。用测试工具稳定画面，不为截图重新引入线上显示设置。

## CI selection and local reproduction / CI 选择与本地复现

```bash
# Documentation checks need only Node.js / 文档检查只需 Node.js
node scripts/check-docs.mjs

# Inspect the complete three-group plan / 查看全量检查的三个分组
CI_CHECK_MODE=full node scripts/ci-plan.mjs

# Inspect daily-publication coverage / 查看每日发布的检查范围
CI_CHECK_MODE=daily node scripts/ci-plan.mjs

# Reproduce one group's suite list from the Actions summary / 复现 Actions 摘要中的某组测试
npm run test:browsers -- --suites browser,language,reading,analytics
```

Documentation-only changes to README, root `AGENTS.md` or Markdown under `docs/` use the lightweight documentation job. It checks local link targets, UTF-8, conflict markers and fenced blocks, but not the truth of the prose or remote URLs. Compare behavior claims against source and distinguish existing functionality from future suggestions. With only documentation edits, this check and `git diff --check` are sufficient; no site rebuild or browser sweep is needed.

只修改 README、根目录 `AGENTS.md` 或 `docs/` 下的 Markdown 时，使用轻量文档任务。检查覆盖本地链接目标、UTF-8、冲突标记和代码围栏，不判断文案事实或远程链接有效性。功能描述仍需与源码比对，并区分已有实现和未来建议。仅修改文档时，运行此检查与 `git diff --check` 即可，无需重建站点或执行浏览器扫描。

The default `npm run test:browsers` still runs every browser suite registered in `package.json`. `--suites` takes exact comma-separated names without the `test:` prefix; unknown, empty, or duplicate names fail immediately. The selector uses the complete push range or the pull request's merge base, including both sides of renames. Missing comparison history and unmapped files trigger all suites. Multiple kinds of targeted changes combine their test selections.

默认的 `npm run test:browsers` 仍执行 `package.json` 中登记的全部浏览器测试。`--suites` 使用逗号分隔的精确名称，不带 `test:` 前缀；未知、空白或重复名称会立即报错。CI 比较整次推送的改动范围，PR 则从共同祖先比较，文件改名的前后路径均计入。无法获取比较历史或出现未映射文件时，会触发全量测试。多类改动会合并各自的检查范围。

Each matrix job runs its assigned suites sequentially; the jobs run in parallel and do not cancel siblings on failure. Reports are available as `test-results-browser-1`, `test-results-browser-2`, and `test-results-browser-3` artifacts, including `visual/` screenshots when the group runs `cross-browser`. WebKit is installed only for those groups; the smaller daily-publication selection remains Chromium-only. All browser jobs use the output of the same build, and deployment waits for the final gate. See [deployment](../operations/deployment.md) for the selection table, schedules, and manual full regression.

每个矩阵任务内部顺序执行分配到的套件，任务之间并行，某一组失败不会取消其他组。报告分别保存在 `test-results-browser-1`、`test-results-browser-2` 和 `test-results-browser-3` 附件中；运行 `cross-browser` 的分组还包含 `visual/` 截图。仅这些分组安装 WebKit，每日发布保留较小的 Chromium 检查范围。各组验证同一次构建的产物，部署等待最后的汇总检查通过。选择规则、定时计划及手动完整回归见[部署说明](../operations/deployment.md)。

## PyTorch series checks / PyTorch 专题检查

`npm run test:pytorch` checks the first installment's C++ model, both animations, language modes, mobile layout, and no-JavaScript reading. Its Python/PyTorch examples require a separate local PyTorch 2.10.0 CPU environment; the website build does not install PyTorch. CPU outputs were checked against that version; CUDA diagrams are source-based illustrations rather than GPU measurements.

`npm run test:pytorch` 检查第一期的 C++ 模型、两组动画、语言模式、手机排版与无 JavaScript 阅读。文中的 Python／PyTorch 示例需要独立的本地 PyTorch 2.10.0 CPU 环境，网站构建不会安装 PyTorch。CPU 输出已按该版本核对；CUDA 图示依据源码，不是 GPU 测量结果。

`npm run test:tensor` checks the second installment’s C++ address model, five tensor layouts, storage aliasing, scan controls, bilingual labels, responsive layout, and no-JavaScript fallback. Its Python outputs were verified with PyTorch 2.10.0 CPU; the standalone C++ program was also compiled and run with GCC 14.2 on Godbolt.

`npm run test:tensor` 检查第二期的 C++ 寻址模型、五种张量布局、存储别名、扫描控件、双语标签、响应式排版和无 JavaScript 回退。Python 输出已在 PyTorch 2.10.0 CPU 下核验，独立 C++ 程序也已通过 Godbolt 的 GCC 14.2 编译并运行。

`npm run test:dispatch` checks part 3’s C++ dispatch model, five execution conditions, inference-mode bypass, finite playback, language switching, responsive layout, and the static diagram without JavaScript. The Python examples and custom-operator checks were verified separately with PyTorch 2.10.0 CPU.

`npm run test:dispatch` 检查第三期的 C++ 调度模型、五种执行条件、inference 模式的绕行路径、有限播放、语言切换、响应式排版与无 JavaScript 静态图。Python 示例和自定义算子检查已在独立的 PyTorch 2.10.0 CPU 环境核验。


`npm run test:autograd` checks part 4’s C++ dependency scheduler against finite differences, both branch orders, shared-node readiness, duplicate leaf edges, finite animation, language switching, mobile layout, and the no-JavaScript graph. Its six Python examples were verified with PyTorch 2.10.0 CPU.

`npm run test:autograd` 用有限差分核对第四期的 C++ 依赖调度模型，并检查两种分支顺序、共享节点就绪条件、重复叶子边、有限动画、语言切换、手机排版与无 JavaScript 图示。六个 Python 示例已在 PyTorch 2.10.0 CPU 中核验。

`npm run test:cuda-streams` checks part 5’s C++ ordering model, same-stream and cross-stream dependencies, the missing-wait hazard, host/GPU completion boundaries, finite playback, language modes, mobile layout and the static fallback. CUDA Python examples require a separate GPU environment; local checks cover syntax and the no-CUDA paths, not GPU execution or performance.

`npm run test:cuda-streams` 检查第五期的 C++ 顺序模型、同流与跨流依赖、缺少等待的危险时序、主机／GPU 完成边界、有限播放、语言模式、手机排版与静态回退。CUDA Python 示例需要独立的 GPU 环境；本地检查覆盖语法和无 CUDA 路径，不代表已验证 GPU 执行或性能。

`npm run test:compile` checks part 6’s C++ fusion/access model, compiler and cache-hit routes, variant reuse, finite playback, language modes, mobile layout and the static graph. The five complete Python examples were verified separately with PyTorch 2.10.0 CPU, including actual Inductor forward/backward compilation, graph capture, generated CPU code, graph breaks and the CPU benchmark. Website CI checks their syntax without installing PyTorch.

`npm run test:compile` 检查第六期的 C++ 融合／访问模型、编译与缓存命中路径、版本复用、有限播放、语言模式、手机排版及静态图。五个完整 Python 示例已在独立的 PyTorch 2.10.0 CPU 环境验证，覆盖实际 Inductor 前向／反向编译、图捕获、生成的 CPU 代码、图中断和 CPU 基准测试。网站 CI 检查其语法，不安装 PyTorch。

Part 7 uses static contract and generation diagrams. Its three Python programs were executed separately with PyTorch 2.10.0 CPU, checking add-family storage behavior, schema inspection, a broadcasting failure and a custom CPU/fake operator. The custom operator's four `opcheck` checks passed for inputs without gradients; this does not establish derivative support or CUDA execution. General article checks cover its layout, language modes and discovery; series tests cover its date and neighboring links. Website CI does not install PyTorch.

第七期使用静态契约与生成流程图。三个 Python 程序已在独立的 PyTorch 2.10.0 CPU 环境执行，核对 add 家族的存储行为、schema 内省、广播失败以及自定义 CPU／fake 算子。自定义算子的四项 `opcheck` 在无梯度输入下通过，不代表已支持求导或 CUDA 执行。通用文章检查覆盖排版、语言模式和内容发现，专题测试覆盖日期与前后期链接。网站 CI 不安装 PyTorch。

## Discussion integration / 讨论功能

`npm run test:discussions` checks article mapping, floating and inline views, draft preservation, language changes, OAuth callback handling, logout, message-source validation, offline behavior, and mobile layout with a mock widget. It does not post to GitHub or verify a real account's login and submission. See [Discussions setup](../operations/discussions.md) for the live integration.

`npm run test:discussions` 使用模拟组件检查文章关联、浮动与文末视图、草稿保留、语言切换、OAuth 回调、退出登录、消息来源校验、离线情况和手机排版。它不会向 GitHub 发帖，也不验证真实账号的登录和发表。实际接入见 [Discussions 配置](../operations/discussions.md)。
