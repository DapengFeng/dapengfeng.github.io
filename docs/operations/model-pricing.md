# Model pricing catalogue / 模型价格目录

[Documentation index / 文档索引](../README.md) · [Deployment / 部署](deployment.md) · [Calculation checks / 测算检查](../development/testing.md#open-model-deployment-economics--开源模型部署经济性)

## Upstream and units / 上游与单位

The calculator uses **the same upstream as CC Switch: models.dev**, a community-maintained catalogue. [CC Switch's source](https://github.com/farion1231/cc-switch/blob/main/src/lib/modelsDev.ts) identifies `https://models.dev/api.json`; [its pricing adapter](https://github.com/farion1231/cc-switch/blob/main/src/lib/modelsDevPricing.ts) imports that catalogue. This website downloads the public feed directly; it does not read a private CC Switch installation, account or configuration. Label the source “CC Switch upstream · models.dev”, never “official provider prices”.

计算器使用 **CC Switch 的同一上游 models.dev**，即社区维护的模型目录。[CC Switch 源码](https://github.com/farion1231/cc-switch/blob/main/src/lib/modelsDev.ts)明确使用 `https://models.dev/api.json`，[价格适配器](https://github.com/farion1231/cc-switch/blob/main/src/lib/modelsDevPricing.ts)从该目录导入。本网站直接获取公开数据，不读取私人 CC Switch 安装、账号或配置。来源标记为“CC Switch 同源 · models.dev”，不得称为“供应商官方定价”。

The [models.dev schema](https://github.com/anomalyco/models.dev#schema) defines `cost.input`, `output` and `cache_read` in **USD per one million tokens**. Preserve those units; do not multiply them by a million or relabel them as CNY. [model-pricing.mjs](../../scripts/lib/model-pricing.mjs) adapts the feed into the same-origin [model-prices.json](../../content/assets/llm-economics/model-prices.json). The reader's workload and budget remain in the browser.

[models.dev 数据定义](https://github.com/anomalyco/models.dev#schema)将 `cost.input`、`output`、`cache_read` 统一定义为**每百万 token 的美元价格**。保持该单位，不再乘一百万，也不换标签冒充人民币。[model-pricing.mjs](../../scripts/lib/model-pricing.mjs)转换数据，浏览器读取同源 [model-prices.json](../../content/assets/llm-economics/model-prices.json)。读者的负载和预算留在浏览器中。

## Supported records and boundaries / 收录规则与边界

Import entries tagged `open_weights:true`, not deprecated, with text-only output modalities and valid input/output prices. Keep each **provider/model pair** separate: two services hosting the same weights can charge differently. Preserve provider documentation, API identifier, canonical model identifier when supplied, context limit, upstream update date and successful retrieval time. The open-weight tag is upstream metadata, not this site's licence audit; an API identifier is not an immutable checkpoint revision.

收录 `open_weights:true`、未弃用、输出模态仅含文本且输入／输出价格有效的条目。每个**供应商／模型组合**独立保留：同一权重由不同服务托管时，报价可能不同。保留供应商文档、API 标识、上游提供的 canonical 标识、上下文上限、上游更新日期与成功获取时间。开放权重标签属于上游元数据，不是本站的许可证审计；API 标识也不等于不可变的 checkpoint revision。

The current worksheet expresses ordinary input, cached input and output. Unknown cache pricing uses the ordinary input rate and `cacheSupported:false`; it is never treated as free caching. Exclude entries containing paid cache writes, context tiers, separate audio rates, a reasoning rate different from the output rate, or unknown cost dimensions. Do not silently omit these charges. The catalogue's `coverage.excludedCounts` records exclusions. Unlisted models can use manually supplied prices and costs.

当前工作表表达普通输入、缓存命中输入和输出三类价格。缓存价格缺失时使用普通输入价并设置 `cacheSupported:false`，不解释为免费缓存。另收缓存写入费、按上下文分档、音频独立计费、推理单价不同于输出价或含未知计费维度的条目暂不收录，不得忽略这些成本后导入。目录中的 `coverage.excludedCounts` 记录排除原因。未收录模型可以手动填写报价和成本。

The upstream's base rate is applied to both time buckets because it does not encode a peak/off-peak calendar. `billing:"flat"` describes this calculator assumption, not a guarantee of all-day flat provider billing. This matters for services such as DeepSeek. Contract discounts, access fees, rate limits, tools and other unrepresented charges need separate verification. Shared source-level `tierDefault`, `conditionsDefault` and `conditionsNoCache` avoid repeating long explanations in thousands of records; the UI derives a selected quote's display text from them.

上游单档价暂用于两个时段，因为它未编码峰谷日历。`billing:"flat"` 描述的是计算器假设，不保证供应商全天同价；DeepSeek 等服务尤其需要留意。合同优惠、接入费、限额、工具及其他未表示的费用需要另外核对。来源级共用 `tierDefault`、`conditionsDefault`、`conditionsNoCache`，避免数千条记录重复长说明；界面据此生成所选报价的提示文字。

## Refresh and failure behavior / 更新与失败行为

`checks.yml` refreshes the catalogue before building, including the daily publication run. Main/master can restore and save a validated catalogue cache; pull requests do not update that shared cache. The checked-in snapshot is the fallback when no cache exists. There is no automatic Git commit. `npm run build` itself uses the available local snapshot without a network request.

`checks.yml` 在构建前更新目录，现有每日发布任务也会运行。main/master 可以恢复并保存有效目录缓存；PR 不更新该共享缓存。没有缓存时回退到仓库中的快照，不自动提交 Git。`npm run build` 本身只使用已有本地快照，不发网络请求。

The updater uses one fixed public URL, a 20-second timeout, a 16 MB response limit, strict JSON/price validation and atomic writes. A network error, malformed response, or unexpected loss of more than half the supported entries preserves the previous validated snapshot. Its successful `fetchedAt` remains unchanged; `lastAttemptAt` and `status:"retained"` record the failed attempt. With no valid fallback, synchronization fails. Retrieval time is neither the upstream edit date nor the provider's price-effective date.

更新器使用一个固定公开网址、20 秒超时、16 MB 响应限制，严格校验 JSON 与价格，并原子写入。网络失败、响应格式异常或支持的条目突然减少一半以上时，保留上次有效快照；原来的成功获取时间 `fetchedAt` 不变，以 `lastAttemptAt` 和 `status:"retained"` 记录失败尝试。没有有效回退时同步失败。获取时间既不是上游编辑日期，也不是供应商价格生效日期。

Readers can filter by provider and search model names or identifiers without changing the current calculation. The UI shows retained snapshots and warns after three days. Loading a refreshed catalogue does not overwrite a reader's current scenario. Selecting a quote imports its rates; changing identity or rates marks it modified. Changing models clears measured throughput; changing currency clears monetary amounts instead of relabeling them. Exports preserve source and edits. Verify the provider's actual tariff and service terms before committing capital.

读者可以按供应商筛选、按模型名称或标识搜索，筛选本身不改变当前计算。界面显示沿用快照的状态，超过三天时提示过期。载入更新后的目录不覆盖读者当前场景；选择报价才导入单价，修改身份或价格后标记已修改。换模型清除实测吞吐，换币种清除金额而非只换标签；导出保留来源与修改信息。投入资金前核对供应商实际价目和服务条款。

## Local commands / 本地命令

```bash
# Network refresh / 联网更新
node scripts/sync-model-pricing.mjs

# Offline catalogue validation / 离线目录校验
node scripts/sync-model-pricing.mjs --check

# Deterministic adapter and fallback tests / 确定性适配与回退测试
node --test tests/model-pricing.test.mjs
```

Fixtures test unit mapping, provider separation, exclusions and fallback behavior. Only an actual successful fetch establishes that the public feed was reachable at that time; fixture tests do not verify live prices.

测试样本验证单位映射、供应商隔离、排除规则及失败回退。只有实际成功抓取才能说明公开数据源当时可达，样本测试不证明实时价格正确。
