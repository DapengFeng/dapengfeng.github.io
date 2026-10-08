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

## Visual acceptance / 视觉验收

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
