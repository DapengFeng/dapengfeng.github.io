# FENG self-hosted typefaces / 自托管字体

These WOFF2 files are renamed, subsetted derivatives of the open-source fonts
below. They are licensed under the SIL Open Font License 1.1, **not** the
website's MIT license. Original copyright and license records remain embedded
in the fonts; the accompanying unmodified OFL files must travel with them.

这些 WOFF2 文件是下列开源字体的改名、子集化版本，遵循 SIL Open Font
License 1.1，**不适用网站的 MIT 许可证**。字体内保留原始版权与许可记录，
分发时须同时保留对应的原始 OFL 文件。

| Web family / 网页字体名 | Upstream / 原字体 | License / 许可证 |
| --- | --- | --- |
| Feng WenKai | [LXGW WenKai Regular / 霞鹜文楷](https://github.com/lxgw/LxgwWenKai) | [wenkai-OFL.txt](wenkai-OFL.txt) |
| Feng Song | [Noto Serif SC / 思源宋体的 Noto 简体中文版本](https://github.com/google/fonts/tree/main/ofl/notoserifsc) | [song-OFL.txt](song-OFL.txt) |
| Feng Garamond | [EB Garamond](https://github.com/google/fonts/tree/main/ofl/ebgaramond) | [garamond-OFL.txt](garamond-OFL.txt) |
| Feng Serif | [Source Serif 4](https://github.com/google/fonts/tree/main/ofl/sourceserif4) | [serif-OFL.txt](serif-OFL.txt) |

## Pinned inputs / 固定来源

- LXGW WenKai: `da0fdf92a0adc5d83c82ecde06d6665c97a2a45e`.
- Google Fonts: `86aced91de05cd4f96a06838124cca45be337b2b`.
- [sources.json](sources.json) records immutable download URLs, source sizes and
  SHA-256 digests. [manifest.json](manifest.json) records embedded upstream
  versions, output digests, byte sizes, weight ranges and character coverage.

[sources.json](sources.json)记录固定下载地址、源文件大小及 SHA-256 摘要；
[manifest.json](manifest.json)记录字体自带的上游版本、产物摘要、体积、字重范围
和实际字符覆盖。

## Modifications and coverage / 修改与覆盖

- Family, full and PostScript names use the new `Feng` names. Original author,
  copyright and OFL records are retained. No endorsement by the original font
  designers is implied.
- WenKai retains its regular weight. Noto Serif SC, EB Garamond and Source Serif
  4 retain variable weights 400–700. EB Garamond and Source Serif 4 include real
  italic faces; Source Serif 4's optical size is fixed at 20 for reading text.
- Chinese coverage includes every supported Hanzi from GB2312 (6,763 standard
  simplified/traditional characters), additional Hanzi used in the site's
  article/topic/interface sources at regeneration, and Chinese/fullwidth
  punctuation. It is **not limited to currently published text** and is not a
  full Unicode CJK collection. Uncovered rare characters use system fallbacks.
- The 1,600 most frequently used Hanzi plus punctuation form a reusable common
  package; the remaining characters are divided into groups of up to 512.
  CSS `unicode-range` loads only packages needed by the rendered page.
- Latin faces cover supported basic and extended Latin, Greek, reading
  punctuation and selected symbols. Hinting is removed, glyph outlines are
  retained, and all assets are compressed as WOFF2. Mathematical formula SVGs
  keep their own MathJax glyphs.

字体内部的家族名、全名和 PostScript 名均改为 `Feng` 名称，保留作者、版权与
OFL 信息，不表示原作者为网站背书。文楷使用常规字重；宋体及两款西文字体保留
400–700 可变字重。西文均提供真正的斜体；Source Serif 4 的光学尺寸固定为 20。

中文覆盖以 GB2312 的 6,763 个标准汉字为基础，补充重新生成时文章、每日数学
题库及界面中使用的其他汉字，以及中文和全角标点。**覆盖不限于现有文章用字**，
但不包含整个 Unicode 汉字集合；生僻字由系统后备字体显示。高频 1,600 字与
标点放在常用包，其他字每包最多 512 个，通过 `unicode-range` 按需加载。

## Regeneration / 重新生成

The checked-in assets are copied by the ordinary Node build. Website visitors
make no requests to Google Fonts or another font CDN; CI needs no font download
or Python installation. Regeneration is a separate maintenance operation:

日常 Node 构建只复制已提交资源，浏览器不请求 Google Fonts 或其他字体 CDN，
CI 无需下载字体或安装 Python。需要调整覆盖时另行执行：

```sh
python3 -m pip install --target /tmp/feng-font-tools fonttools==4.68.0 brotli==1.2.1
PYTHONPATH=/tmp/feng-font-tools python3 scripts/prepare-fonts.py --cache /tmp/feng-font-sources
```

The script verifies each pinned source checksum before writing output. Keep the
generated CSS, WOFF2 files, source manifest, output manifest and licenses
together. CSS uses `font-display: swap`; do not preload every CJK package.

脚本写入产物前会核对每个固定来源的摘要。字体、CSS、来源与产物清单、许可证
须一起维护。CSS 使用 `font-display: swap`，不应预加载全部汉字分包。
