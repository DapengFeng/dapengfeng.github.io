# Writing with clarity and attention / 清楚而细致的写作

[Documentation index / 文档索引](../README.md) · [Website mission / 网站使命](../../README.md)

Writing serves two parts of the mission: help readers understand something difficult, and preserve something they might otherwise overlook. Describe why the author writes and how the work helps readers; specific disciplines, tools and destinations belong in articles and reading paths, so future subjects fit without redefining the author. Keep the voice personal and concrete, without inventing credentials or promising unpublished coverage.

写作服务于使命的两个方面：帮助读者理解原本难懂的事物，留下容易被忽略的细节。作者表达重点放在为什么写、如何帮助读者；具体学科、工具与目的地放在文章和阅读路径中，使未来增加内容时不必重新定义作者身份。语气保持个人化、具体，不虚构资历，也不承诺尚未发布的内容。

These standards guide new work and revision; they do not claim every existing page already meets them or require the same headings, length or voice. Accuracy and faithfulness to the supplied material are requirements. Rhythm, openings and structure should fit the subject. See [publishing](publishing.md) for metadata and reading-entry implementation, and [presentation](presentation.md) for visual organization.

以下标准用于新增内容与修改，不代表每个现有页面都已达到，也不要求相同的标题、篇幅或声音。准确和忠于材料是底线；节奏、开头与结构应随内容变化。元信息和阅读入口的实现见[文章发布](publishing.md)，视觉组织见[排版与交互规范](presentation.md)。

## Decide what the reader will gain / 先想清读者能得到什么

Before drafting, write a short working note: who is reading, what they already know, what question or experience holds the piece together, and what they should understand or notice afterward. Identify the evidence, photograph or example that can carry that promise. This is a drafting aid, not a block that must appear on the page. Narrow an oversized promise instead of padding the article to match it.

动笔前用几句工作笔记说明：读者是谁、已有何种基础、什么问题或体验贯穿全文、读完能理解或注意到什么，再找出支撑这个承诺的证据、照片或例子。这是写作工具，不必作为说明框展示在页面上。承诺过大时收窄主题，不为凑成“全面介绍”填充内容。

Treat the title, listing summary and opening as different parts of that promise. The title identifies the subject and angle. The short summary states the useful relationship or distinctive observation. The opening gives readers a concrete reason to continue. Do not repeat the same sentence in all three or advertise animations and tools in place of the subject.

标题、列表摘要和开头共同兑现一个承诺，但各有职责：标题指出对象与角度，短摘要交代有用的关系或独特观察，开头用具体内容让读者愿意继续。不在三处重复同一句话，也不用动画与工具的介绍代替内容本身。

For example, a title such as “How light becomes a neural signal” promises a bounded mechanism; “Everything about human vision” promises far more. A matrix-multiplication summary can say that changing loop order changes data reuse and cache access, rather than merely announce an “in-depth performance analysis”. Follow the existing [summary metadata](publishing.md), without adding new fields for this planning step.

例如，“光如何变成神经信号”承诺解释一段具体机制，“一篇读懂人类视觉的一切”则承诺过多。矩阵乘法的摘要可以说“循环顺序改变数据复用和缓存访问”，不必写成“深入浅出地全面分析性能优化”。沿用现有[摘要元信息](publishing.md)，无需为这一步规划新增字段。

A technical card should state one concrete relationship, in roughly 20–35 English words or one short Chinese sentence; a travel card should name the places and observations. Introduce detailed terminology where the body explains its role. For example, a compiler card can say: “Compiled code can be reused while its input assumptions hold; changing those assumptions may trigger compilation again.”

技术卡片用一句话说明具体关系，英文约 20–35 词，中文尽量简短；游记卡片交代地点与观察。具体术语留到正文解释其职责时再展开。例如，编译文章的卡片可写：“输入仍满足缓存版本的条件时，编译代码才能复用；条件变化可能触发重新编译。”

## Organize the content before filling the page / 先组织内容，再填页面

Arrange sections by the dependencies of understanding or the progression of an experience, not by how much material was collected. Give prerequisites before they are used and connect each new level to the previous one. In observations, let place, time or a recurring image guide movement without inventing chronology. A section should add a step, evidence or perspective rather than repeat the introduction.

章节顺序依据理解上的依赖或体验的推进，不依据收集了多少材料。先交代必要前提，并把新的层次接回前一层；观察中可由地点、时间或反复出现的意象引路，但不虚构行程。一节应增加步骤、证据或视角，不能重复开头。

Keep a readable main path, with optional derivations, full programs or secondary comparisons nearby when useful. Do not hide a necessary premise in a collapsed panel or require an interaction or external site to obtain the main explanation. Each article should be understandable at its stated prerequisite level; a series can link to earlier foundations instead of repeating them in full.

保留清楚的主阅读线，必要时在近旁提供可选推导、完整程序或次要对比；理解结论必需的前提不能藏在折叠区，也不应要求读者操作控件或跳到外站才能获得主要解释。每篇文章在所需基础范围内应能独立读懂，专题可链接已有基础，避免整段重复。

Make navigation describe the published content it leads to. Categories group subjects, series establish a learning sequence, and related reading connects questions. Do not add a category for a single keyword or label unrelated articles as a course. Representative reading should demonstrate the mission and offer a useful first step; keep evidence, stable references and meaningful next reading available for readers who want to revisit or continue.

导航应准确描述它指向的已发布内容：分类组织主题，专题建立学习顺序，相关阅读连接问题。不为一个关键词就新增分类，也不把无关文章包装成课程。代表作应体现使命、提供合适的阅读起点；保留依据、稳定引用与有联系的延伸阅读，方便读者核对、回看或继续探索。

Depth should follow the question. A short observation need not acquire an experiment, and a technical derivation need not acquire a travel-style ending. Choose a form that provides enough evidence and room to think, then stop when the promise is fulfilled. Article count, word count, animation count and repeated engagement prompts are not content goals.

深度由问题决定。短篇观察不必强加实验，技术推导也不必加上游记式结尾。选择能给读者充分依据与思考空间的形式，兑现承诺后即可收束。不以文章数量、字数、动画数量或反复互动提示作为内容目标。

## Build paragraphs that move the explanation / 让段落推进理解

Give each paragraph one main job: establish an observation, explain a cause, interpret evidence or make a transition. Begin from something the preceding paragraph made available, add the missing step, and leave the reader ready for the next one. A short sentence can settle a point; a longer one can carry a relationship. Split a sentence when its nested qualifications obscure who does what, but avoid reducing connected prose to a stack of slogans.

每段承担一个主要任务：交代观察、解释原因、解读证据或完成过渡。从上一段已经建立的认识出发，补上缺的一步，再为下一段留下落点。短句可收住判断，长句可展开关系；层层限定使主语和动作模糊时就拆开，但不要把连贯文章拆成口号堆叠。

Prefer named objects and precise verbs to abstract praise. Replace “greatly improves the experience” with what changes for the reader or system. Replace “as is well known” with the explanation or source the claim needs. Remove sentences that only say “this is important”, “we will explore” or “it is worth noting” when the following sentence can state the point directly. Use lists for genuinely parallel items and prose for causal or narrative development.

用明确对象与准确动词替代抽象赞美。“大幅提升体验”应改成读者或系统具体发生什么变化；“众所周知”应改成论述需要的解释或来源。如果下一句已经能直接说明问题，就删掉只起铺垫作用的“这一点非常重要”“接下来深入探讨”“值得注意的是”。并列信息适合列表，因果与叙事展开适合连续文字。

Make transitions explain the dependency. “We know where the image forms; the next question is how a cell changes its electrical state” connects two problems. “Next, we introduce phototransduction” only announces a section. Headings should help readers recover the argument when scanning, without turning every paragraph into a numbered subsection.

过渡要交代前后依赖。“知道像落在哪里之后，还需要解释细胞怎样改变电状态”连接了两个问题；“接下来介绍光转导”只是宣布章节。小标题应让读者扫读时找回思路，不必把每段都切成编号小节。

## Explain mechanisms without skipping the hard step / 不跳过最难懂的一步

Carry one question and a stable example through the explanation. Start with a phenomenon or observable behavior, explain its physical, biological or computational mechanism, then show what a model preserves and leaves out. Introduce a term when the reader needs it, explain its role in the example, and use it consistently. Return to the opening question using the evidence already shown.

围绕一个问题，用稳定的例子贯穿讲解。从现象或可观察行为开始，解释物理、生物或计算机制，再说明模型保留和省略了什么。读者需要某个术语时再引入，说明它在例子里的作用，并保持用词一致。最后用已经展示的证据回答开头问题。

For example, the visual-system article follows a red cup from optical imaging to cell responses, neural signals and engineering operations. If light reduces transmitter release and two cell types respond in opposite ways, naming molecules alone does not explain the reversal: show the sign of each change and the responsible connection. Align signal traces in time, and show a center–surround filter's center value, surround value and subtraction before its formula or code.

例如，视觉系统一文以红杯为线索，从光学成像进入细胞响应、神经信号与工程运算。光照减少递质释放、两类细胞却作出相反响应时，单列分子名称无法解释反号：要展示每一步变化的方向和产生差异的连接。信号曲线在时间上对齐；中心—周围滤波先展示中心值、周围值和相减过程，再给公式或代码。

Give equations a reading path: what the symbols refer to, the assumptions or units that matter, what relation the equation expresses, and a small consequence or example. Derivations should explain why the next step follows. Code should implement an already understandable operation and show meaningful input and output. Keep notation, figure labels and variable names aligned. An analogy is useful only if its limits are clear; biological responses, physical models and engineering approximations must not silently become identical claims.

为公式提供阅读路径：符号指什么、哪些假设或单位重要、公式表达什么关系，以及一个小推论或例子。推导应解释下一步为何成立，代码应落实已能理解的操作，并展示有意义的输入与输出。符号、图中标签与变量名保持对应。类比需要交代边界，不能把生物响应、物理模型和工程近似悄悄写成同一件事。

Separate what is observed, derived, checked in source code, executed, measured, inferred or still uncertain. Cite a source beside the claim it supports; use primary sources for technical mechanisms and state relevant versions or conditions. A benchmark needs its workload, environment and measurement boundary. Do not turn an illustrative animation into measured timing or a local test into a universal result. Evidence should support the sentence, not merely decorate a reference list.

分清观察、推导、源码核对、实际执行、测量、推断和尚未确定的部分。来源靠近其支撑的论述，技术机制优先使用原始资料，并交代相关版本或条件；性能结论要说明工作负载、环境与计时边界。不把示意动画当作实测时序，也不把局部测试扩大为普遍结论。证据要支撑具体句子，不只是装饰参考文献列表。

## Write observation before declaring emotion / 先写观察，再让情绪发生

In travel and personal essays, choose the details that distinguish this place or moment: the lettering of a sign, the spacing of boats, the color of water against stone. Let a change of scene, scale or pace carry the reader. Vary longer descriptive sentences with shorter ones that let an image settle. Do not attach a metaphor to every object; a precise detail often needs no embellishment.

游记与随笔选择能区别这个地方、这个片刻的细节：招牌上的字、船只之间的距离、水色与石头的关系。用场景、观察尺度或节奏的变化带读者往前走。较长的描写之后，可以用短句让画面停住。不为每样东西都配比喻，准确的细节常常已经足够。

A photograph supports visible detail, not a claim about an unheard sound, a taste or the author's feelings. Use supplied recollections for those experiences; do not invent weather, conversations, food locations or memories to complete a narrative. If a needed fact is missing, omit it or ask for it when essential. Let emotion emerge from the selection and arrangement of details; avoid a technical report's scaffolding and stock endings about “healing”, “the poetry of distance” or what every journey supposedly teaches.

照片能够支撑可见细节，不能证明未听见的声音、味道或作者当时的心情；这些体验需以提供的回忆为依据，不为补齐叙事虚构天气、对话、美食地点或经历。缺少事实时可以省略，确实影响内容时再询问。情绪通过细节的取舍与安排自然形成，避免技术报告式结构，也不套用“治愈了心灵”“诗与远方”或“每次旅行都教会我们”的结尾。

An ending can return to an earlier object, pause at a specific scene or answer the piece's question. It need not elevate a small experience into a general life lesson. Images should support the adjacent passage; short captions can identify the scene without repeating the paragraph. See [presentation](presentation.md) for figure and layout choices.

结尾可以回到前文的一个物件，停在具体场景，或回应贯穿全文的问题，不必把一次小经历拔高成人生道理。图片支持紧邻段落，简短图注交代场景而不重复正文。图示与布局选择见[排版与交互规范](presentation.md)。

Preserve the owner's location corrections in existing work: in the Chaoshan essay, “给阿嫲的情书” belongs to Chaozhou's Paifang Street, and the food photographs were taken in Chaozhou. Its Chaozhou section describes sightseeing before food, with the food heading “碟中有余味”. Retain these facts when rearranging photographs or prose; visual similarity is not a reason to move a scene to another city.

维护旧文章时保留作者纠正过的地点：潮汕游记中，“给阿嫲的情书”位于潮州牌坊街，美食照片也拍于潮州。潮州部分先写游览、后写美食，美食小节标题为“碟中有余味”。调整图文顺序时保留这些事实，不因画面相似就把场景移到其他城市。

## Revision examples / 修改示例

These are editing examples, not new factual claims or memories to insert into articles. Check the evidence before using a specific detail.

以下用于示范修改方法，不是可以直接加入文章的新事实或回忆。采用具体细节前仍须核对材料。

| Aim / 目的 | Weak draft / 初稿问题 | Revision direction / 修改方向 |
| --- | --- | --- |
| Explain a mechanism / 解释机制 | “Caching significantly improves efficiency.”<br>“缓存显著提高效率。” | “When the next operation reuses data already in the cache, it can avoid another main-memory access.” State when reuse actually occurs.<br>“下一次运算复用仍在缓存中的数据时，可以减少一次主存访问。”继续说明复用何时发生。 |
| Replace general praise with observation / 用观察替换赞美 | “The old street is full of history and charm.”<br>“老街充满历史底蕴，别有韵味。” | “The shop signs hang out from the arcades; looking up, the characters overlap along the street.” Use only if the photograph supports it.<br>“招牌从骑楼下伸出来，抬头看，字一层叠着一层。”仅在照片支持时采用。 |
| Interpret rather than repeat a figure / 解读图示 | “The figure shows two curves.”<br>“如图所示，有两条曲线。” | Name what is held constant and where the curves diverge, then explain what that difference means.<br>说明什么条件保持不变、曲线从哪里开始分开，以及差异意味着什么。 |

## Revise in passes / 分层修改

1. **Purpose and truth:** does the draft deliver its promise, and is each factual or personal claim supported? Remove digressions, repair missing reasoning and narrow claims before polishing sentences.

   **主旨与事实**：正文是否兑现开头承诺，事实与亲身经历是否有依据？先删离题内容、补推理缺口、收窄结论，再润色句子。

2. **Sequence and connection:** can a reader follow the example without guessing a missing term, step or location? Check that each figure and section advances the same explanation or experience.

   **顺序与关联**：读者是否需要猜测术语、中间步骤或地点才能读下去？检查每张图、每一节是否推进同一段理解或体验。

3. **Language and rhythm:** remove empty praise, repeated introductions and unnecessary jargon; vary sentence length and read the paragraph aloud. Preserve a distinctive observation or judgment instead of smoothing every sentence into generic prose.

   **语言与节奏**：删空泛赞美、重复铺垫和非必要术语，调整长短句，并通读段落。保留有辨识度的观察与判断，不把每一句磨成通用文案。

4. **Bilingual reading:** maintain the same claims, examples, emphasis and limits in both languages, while writing natural sentences in each. Recheck titles, captions and interface labels as well as body text; do not translate Chinese literary expressions word for word. Follow the [language rules](language.md).

   **双语阅读**：两种语言保持相同的论述、例子、重点与边界，同时各自成句自然。除正文外，还要核对标题、图注和界面标签；中文文学表达不逐字硬译。遵循[双语规范](language.md)。

5. **The rendered page:** read text beside its figures at mobile and desktop widths. A sentence, caption or transition that works in a source file may become unclear when separated by a large image. Use [presentation standards](presentation.md) and [visual validation](../development/testing.md#visual-acceptance--视觉验收) for the final reading pass.

   **成品页面**：在手机与桌面尺寸中对照图文阅读。源码里顺畅的句子、图注或过渡，隔着大幅图片可能就不再清楚。最后按[排版与交互规范](presentation.md)与[视觉验收](../development/testing.md#visual-acceptance--视觉验收)检查真实阅读体验。
