#!/usr/bin/env python3
"""Render the DeepSeek self-hosting article's illustrative financial figures.

Requires Python, NumPy, Matplotlib, Noto Sans CJK SC, and Node.js.
Run from the repository root: python3 scripts/render-deepseek-roi.py

All financial results are evaluated by content/assets/deepseek-roi/model.js.
This script only samples and plots that model; it does not implement financial
formulas, measure inference performance, or assert a hardware throughput rate.
SVGs have English/Chinese and desktop/mobile variants. Inspection PNGs are
written outside the repository; the only shipped PNG is the abstract cover.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "deepseek-roi-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MultipleLocator
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PAPER, INK, MUTED = "#111e2b", "#f2f1e8", "#bdcbd4"
GRID, BLUE, MINT, AMBER, SLATE = "#435868", "#87cce9", "#b2dfc0", "#e9b980", "#8ca3b4"
PREVIEW = Path(tempfile.gettempdir()) / "deepseek-roi-figures"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Noto Sans CJK SC", "DejaVu Sans"],
    "svg.fonttype": "path", "svg.hashsalt": "feng-deepseek-roi",
    "axes.unicode_minus": False, "lines.solid_capstyle": "round",
})


def samples():
    """Use Node's authoritative evaluator, including all scenario assumptions."""
    source = r"""
import fs from 'node:fs';
import {defaults,evaluate} from './content/assets/deepseek-roi/model.js';
const scenarioData=JSON.parse(fs.readFileSync('./content/assets/deepseek-roi/scenarios.json','utf8'));
for (const scenario of scenarioData.scenarios) {
  const current=evaluate(scenario.overrides);
  if (Math.abs(current.NPV-scenario.result.NPV)>1e-6) {
    throw new Error('scenarios.json is stale; run scripts/calculate-deepseek-roi.mjs before rendering');
  }
}
const base=evaluate();
const capacityMultipliers=Array.from({length:201},(_,i)=>i/80);
const ceiling=[[24e9,6e9],[48e9,12e9],[64e9,16e9]].map(([monthlyInputTokens,monthlyOutputTokens])=>({
  monthlyInputTokens,monthlyOutputTokens,
  atDefault:evaluate({monthlyInputTokens,monthlyOutputTokens}).maxInitialInvestment,
  values:capacityMultipliers.map(multiplier=>evaluate({
    monthlyInputTokens,monthlyOutputTokens,
    inputTokensPerSecond:defaults.inputTokensPerSecond*multiplier,
    outputTokensPerSecond:defaults.outputTokensPerSecond*multiplier,
  }).maxInitialInvestment),
}));
console.log(JSON.stringify({defaults,scenarioData,base,ceiling,capacityMultipliers}));
"""
    return json.loads(subprocess.check_output(["node", "--input-type=module", "-e", source], cwd=ROOT))


def canvas(width, height):
    fig = plt.figure(figsize=(width / 72, height / 72), facecolor=PAPER)
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(height, 0))
    ax.set_axis_off()
    return fig, ax


def text(ax, x, y, label, size=24, color=INK, **kwargs):
    return ax.text(x, y, label, fontsize=size, color=color, va="top", linespacing=1.28, **kwargs)


def segment(ax, x, y, label, color, size, style="-", width=42):
    ax.plot([x, x + width], [y + size * .65] * 2, color=color, lw=3.1, ls=style)
    text(ax, x + width + 12, y, label, size, color)


def axes(fig, bounds, width, height, mobile):
    left, top, right, bottom = bounds
    ax = fig.add_axes([left / width, 1 - bottom / height, (right-left) / width, (bottom-top) / height], facecolor=PAPER)
    ax.tick_params(axis="both", colors=MUTED, labelsize=25 if mobile else 23, pad=8, length=0)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    for side in ["bottom", "left"]:
        ax.spines[side].set_color(GRID)
    ax.grid(axis="y", color=GRID, alpha=.48, lw=.8)
    ax.set_axisbelow(True)
    return ax


def save(fig, out, stem, lang, mobile, description):
    name = f"{stem}{'-mobile' if mobile else ''}-{lang}"
    target = out / f"{name}.svg"
    fig.canvas.draw()
    width, height = fig.get_size_inches() * fig.dpi
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text().strip():
            extent = artist.get_window_extent(fig.canvas.get_renderer())
            if extent.x0 < -1 or extent.y0 < -1 or extent.x1 > width + 1 or extent.y1 > height + 1:
                raise ValueError(f"Text outside {name}: {artist.get_text()!r}, {extent.bounds}")
    fig.savefig(target, format="svg", dpi=72, metadata={"Date": None, "Description": description})
    # Add screen-reader metadata while keeping deterministic geometry and paths.
    tree = ET.parse(target)
    root = tree.getroot()
    ns = "{http://www.w3.org/2000/svg}"
    title = ET.Element(ns + "title", {"id": name + "-title"})
    title.text = description
    root.insert(0, title)
    root.set("role", "img")
    root.set("aria-labelledby", name + "-title")
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    tree.write(target, encoding="utf-8", xml_declaration=True)
    # Matplotlib's whitespace-only SVG lines otherwise retain trailing spaces.
    svg_lines = target.read_text(encoding="utf-8").splitlines()
    target.write_text("\n".join(line.rstrip(" \t") for line in svg_lines) + "\n", encoding="utf-8")
    PREVIEW.mkdir(parents=True, exist_ok=True)
    fig.savefig(PREVIEW / f"{name}.png", dpi=72)
    plt.close(fig)


def bridge_figure(data, out, lang, mobile):
    """Show which part of the avoided API bill can recover initial investment."""
    zh = lang == "zh"
    w, h = (432, 1040) if mobile else (1080, 760)
    fig, page = canvas(w, h)
    margin, fs = (24, 25) if mobile else (36, 24)
    row = data["base"]["firstMonth"]
    api, fallback = row["fullApiCost"], row["overflowApiCost"]
    avoided, opex, savings = row["avoidedApiCost"], data["defaults"]["monthlyOpex"], row["operatingSavings"]
    title = ("API 账单中有多少\n能用于收回初投？" if zh else "How much of the bill\ncan repay the purchase?") if mobile else ("从每月 API 账单，到可用于收回初投的金额" if zh else "From the monthly API bill to cash that can recover the investment")
    text(page, margin, 24, title, 29 if mobile else 31)
    subtitle = ("月输入 240 亿 · 输出 60 亿 token\n基线组合 · 单位：元/月" if zh else "24B input · 6B output tokens/mo\nB = billion tokens · RMB/month") if mobile else ("月输入 240 亿 token · 月输出 60 亿 token · 相同缓存组合 · 单位：元/月" if zh else "24B input + 6B output tokens/month · B = billion tokens · same cache mix · RMB/month")
    text(page, margin, 118 if mobile else 81, subtitle, fs if mobile else 23, MUTED)
    labels = (["全 API 账单", "停机回退 API", "可避免的 API 账单", "本地运营费用", "可用于收回初投"] if zh else ["All-API bill", "API fallback", "Avoided API bill", "Local operations", "Remaining savings"])
    colors = [BLUE, AMBER, BLUE, AMBER, MINT]
    low = [0, avoided, 0, savings, 0]
    high = [api, api, avoided, avoided, savings]
    amounts = [f"{api:,.0f}", f"−{fallback:,.2f}", f"{avoided:,.2f}", f"−{opex:,.0f}", f"{savings:,.2f}"]
    if mobile:
        bar_left, bar_width, y0, gap = 24, 384, 214, 124
        for i in range(5):
            y = y0 + i * gap
            text(page, 24, y, labels[i], 25, colors[i])
            text(page, 408, y + 38, amounts[i], 27, colors[i], ha="right")
            by = y + 78
            page.plot([bar_left, bar_left + bar_width], [by + 13] * 2, color=GRID, lw=.8, alpha=.5)
            x1, x2 = bar_left + bar_width * low[i] / api, bar_left + bar_width * high[i] / api
            page.add_patch(Rectangle((x1, by), x2 - x1, 26, facecolor=colors[i], edgecolor="none"))
        text(page, 24, 883, f"每月剩余 {savings:,.2f} 元收回初投。\n未来账单、运营成本与使用年限，\n决定最终能收回多少。" if zh else "The balance can repay the\nupfront investment. Future bills,\ncosts and useful life matter.", 25, MUTED)
    else:
        ax = axes(fig, (88, 204, 1040, 582), w, h, False)
        ax.set_xlim(-.55, 4.55)
        ax.set_ylim(0, 52500)
        ax.set_yticks([0, 15000, 30000, 45000])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: "0" if y == 0 else f"{y/1000:g}k"))
        ax.set_xticks([])
        ax.spines["left"].set_visible(False)
        for i in range(5):
            ax.bar(i, high[i] - low[i], bottom=low[i], width=.57, color=colors[i], zorder=3)
            if i < 4:
                remaining = [api, avoided, avoided, savings][i]
                ax.plot([i + .285, i + 1 - .285], [remaining, remaining], color=MUTED, lw=1.1, ls=(0, (3, 4)), alpha=.8)
            offset = 3600 if i == 1 else 900
            ax.text(i, high[i] + offset, amounts[i], color=colors[i], fontsize=24, ha="center", va="bottom")
            if i == 1:
                ax.plot([i, i], [high[i] + 500, high[i] + 3300], color=AMBER, lw=1.1)
            label = labels[i]
            if zh:
                label = ["全 API\n账单", "停机回退\nAPI", "可避免的\nAPI 账单", "本地\n运营费用", "用于收回\n初始投资"][i]
            else:
                label = ["All-API\nbill", "Downtime\nAPI fallback", "Avoided\nAPI bill", "Local\noperations", "To recover\nthe investment"][i]
            ax.text(i, -3200, label, color=colors[i], fontsize=23, ha="center", va="top", linespacing=1.3, clip_on=False)
        text(page, 36, 690, f"{savings:,.2f} 元/月用于收回初投；未来账单、运营成本与使用年限决定最终回收额。" if zh else f"RMB {savings:,.2f}/month remains. Future bills, operating costs and lifetime determine the return.", 21, MUTED)
    save(fig, out, "cash-bridge", lang, mobile,
         f"Illustrative baseline: 24 billion input and 6 billion output tokens per month. Monthly cash bridge: full API bill RMB {api:,.2f}, less fallback API RMB {fallback:,.2f}, leaves avoided API cost RMB {avoided:,.2f}. Less RMB {opex:,.2f} local operations leaves RMB {savings:,.2f} to recover initial investment. All amounts come from the shared model. Bar lengths use a common monetary scale; this is not measured performance or realized profit.")


def cashflow_figure(data, out, lang, mobile):
    zh = lang == "zh"
    w, h = (432, 1040) if mobile else (1080, 730)
    fig, page = canvas(w, h)
    margin, fs = (24, 25) if mobile else (36, 24)
    text(page, margin, 24, ("何时收回购机现金？" if zh else "When does the cash\ncome back?") if mobile else ("相对全 API：36 个月累计现金流" if zh else "Cumulative cash flow against the all-API alternative"), 29 if mobile else 31)
    text(page, margin, 119 if mobile else 80,
         "月 0 支出 160 万元" if zh else "Month 0: ¥1.60 million outlay", fs, MUTED)
    definitions = [
        ("base-24b-6b", BLUE, "-", "输入 240 亿 · 输出 60 亿/月", "24B input · 6B output / month"),
        ("high-64b-16b", MINT, "-", "输入 640 亿 · 输出 160 亿/月", "64B input · 16B output / month"),
        ("base-api-minus-30", AMBER, "--", "240 亿入 · 60 亿出/月\nAPI 价格每年降 30%", "24B in · 6B out / month\nAPI prices fall 30%/year"),
    ]
    by_id = {r["id"]: r for r in data["scenarioData"]["scenarios"]}
    text(page, margin, 162 if mobile else 121,
         "以下用量均为 token" if zh else "All volumes in tokens; B = billion", fs if mobile else 22, MUTED)
    for i, (sid, color, style, zlabel, elabel) in enumerate(definitions):
        segment(page, margin if mobile else 36 + (i % 2) * 525, (207 + i * 49) if mobile else 166 + (i // 2) * 43,
                zlabel if zh else elabel, color, 23 if mobile else 23, style, 24 if mobile else 36)
    ax = axes(fig, (85, 430, 408, 782) if mobile else (106, 323, 1040, 582), w, h, mobile)
    all_y = []
    for sid, color, style, _, _ in definitions:
        rows = by_id[sid]["result"]["months"]
        values = [r["cumulativeOperatingNet"] / 1e6 for r in rows]
        all_y.extend(values)
        ax.plot([r["month"] for r in rows], values, color=color, lw=3.2, ls=style)
        ax.plot(36, values[-1], "o", color=color, ms=6, clip_on=False)
    ymin = np.floor(min(all_y) * 2) / 2
    ymax = max(.5, np.ceil(max(all_y) * 2) / 2)
    ax.set_ylim(ymin - .1, ymax + .1)
    ax.set_xlim(0, 36)
    ax.set_xticks([0, 12, 24, 36])
    ax.yaxis.set_major_locator(MultipleLocator(1 if mobile else .5))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{y:g}"))
    ax.axhline(0, color=INK, lw=1.5, alpha=.95)
    payback = by_id["high-64b-16b"]["result"]["simplePaybackMonth"]
    crossing = by_id["high-64b-16b"]["result"]["months"][payback]["cumulativeOperatingNet"] / 1e6
    ax.plot(payback, crossing, "o", color=MINT, markeredgecolor=INK, markeredgewidth=1.5, ms=8, zorder=6)
    ax.annotate(f"{payback} 月回本" if zh else f"Month {payback}", (payback, crossing),
                xytext=(17 if mobile else 22, .40), color=MINT, fontsize=25 if mobile else 23,
                arrowprops={"arrowstyle": "-", "color": MINT, "lw": 1.2}, va="center",
                bbox={"facecolor": PAPER, "edgecolor": "none", "pad": 2})
    text(page, 85 if mobile else 106, 381 if mobile else 280, "百万元" if zh else "¥ million", fs, MUTED)
    text(page, 245 if mobile else 573, 832 if mobile else 627, "月份" if zh else "Month", fs, MUTED, ha="center")
    foot = ("先扣购机款，再累计运营节省。\n未折现，不计期末残值。\n白线为现金回本门槛。" if zh else "The curves are undiscounted.\nNo resale value is included.\nWhite line: cash payback.") if mobile else ("月 0 扣购机款，之后累计相对全 API 的运营节省；未折现，也不计期末残值。" if zh else "Purchase at month 0, then savings versus all API. Undiscounted; excludes resale.")
    text(page, margin, 903 if mobile else 678, foot, fs if mobile else 22, MUTED)
    save(fig, out, "cumulative-cashflow", lang, mobile,
         "Illustrative cumulative undiscounted cash flow relative to all API for three monthly token scenarios: 24 billion input plus 6 billion output, 64 billion input plus 16 billion output, and the 24B-input/6B-output workload with API prices falling 30 percent annually. Each starts at minus 1.6 million RMB. The 64B-input/16B-output scenario first repays the investment at month 32. Operating savings exclude resale value. These are assumptions, not measured hardware performance.")


def investment_figure(data, out, lang, mobile):
    """Plot the conditional break-even investment, read directly from the model."""
    zh = lang == "zh"
    w, h = (432, 1100) if mobile else (1080, 800)
    fig, page = canvas(w, h)
    margin, fs = (24, 25) if mobile else (36, 24)
    title = ("吞吐与需求，决定\n最多值得投入多少" if zh else "What investment can\nthis demand support?") if mobile else ("合格吞吐与需求，决定可承受的初始投资" if zh else "What initial investment can the demand support?")
    text(page, margin, 24, title, 29 if mobile else 31)
    subtitle = ("36 个月 · 年折现率 10%\n月运营 2.4 万元 · 残值 24 万元" if zh else "36 months · 10%/year discount\n¥24k/month OpEx · ¥240k resale") if mobile else ("36 个月 · 年折现率 10% · 运营费用固定为 2.4 万元/月 · 期末残值固定为 24 万元" if zh else "36 months · 10% annual discount · fixed ¥24,000/month OpEx · fixed ¥240,000 resale")
    text(page, margin, 118 if mobile else 81, subtitle, fs if mobile else 22, MUTED)
    colors = [BLUE, SLATE, MINT]
    styles = ["-", (0, (5, 3)), "-"]
    for i, (row, color, style) in enumerate(zip(data["ceiling"], colors, styles)):
        incoming, outgoing = row["monthlyInputTokens"], row["monthlyOutputTokens"]
        amount = row["atDefault"] / (1e4 if zh else 1e6)
        value = f"1× 时：{amount:.1f} 万元" if zh else f"At 1×: ¥{amount:.3f}m"
        label = (f"入 {incoming/1e8:g} 亿 · 出 {outgoing/1e8:g} 亿/月" if zh
                 else f"{incoming/1e9:g}B in · {outgoing/1e9:g}B out / mo") + "\n" + value
        segment(page, margin if mobile else 36 + i * 347, 202 + 63*i if mobile else 139,
                label, color, 25 if mobile else 22, style, 24 if mobile else 26)
    segment(page, margin, 399 if mobile else 224,
            "160 万元假设预算" if zh else "¥1.60m assumed budget", AMBER, fs if mobile else 23, "--", 24 if mobile else 26)
    text(page, 24 if mobile else 106, 452 if mobile else 275,
         "初投上限 / 百万元" if zh else "Investment limit / ¥m", fs if mobile else 23, MUTED)
    ax = axes(fig, (75, 504, 408, 770) if mobile else (106, 323, 1040, 600), w, h, mobile)
    x = np.array(data["capacityMultipliers"])
    ax.axhspan(-.75, 0, color="#192938", alpha=.65)
    for row, color, style in zip(data["ceiling"], colors, styles):
        ax.plot(x, np.array(row["values"])/1e6, color=color, lw=3.1, ls=style, zorder=3)
    ax.axvline(1, color=MUTED, lw=1.1, ls=(0, (2, 4)), alpha=.75)
    for row, color in zip(data["ceiling"], colors):
        ax.plot(1, row["atDefault"]/1e6, "o", ms=7, color=color, markeredgecolor=INK, markeredgewidth=1.2, zorder=5)
    ax.axhline(data["defaults"]["capex"]/1e6, color=AMBER, lw=1.8, ls=(0, (5, 4)), zorder=4)
    ax.axhline(0, color=INK, lw=1.0, alpha=.75)
    ax.set_xlim(0, 2.5)
    ax.set_ylim(-.75, 3.5)
    ax.set_xticks([0, .5, 1, 1.5, 2, 2.5])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))
    ax.set_yticks([0, 1, 2, 3])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{y:g}"))
    ax.text(1, 3.43, "1×", color=MUTED, ha="center", va="top", fontsize=25 if mobile else 22)
    plateau = data["ceiling"][0]["atDefault"]/1e6
    ax.annotate("基线平台" if zh else "Base plateau", (2.15, plateau), xytext=(1.65, .29),
                color=BLUE, fontsize=25 if mobile else 23, ha="center", va="center",
                arrowprops={"arrowstyle": "-", "color": BLUE, "lw": 1.1},
                bbox={"facecolor": PAPER, "edgecolor": "none", "pad": 1})
    text(page, 240 if mobile else 573, 817 if mobile else 647,
         "成对容量倍率" if zh else "Paired capacity multiplier", fs, MUTED, ha="center")
    basis = ("1×：32k 输入 token/s\n       8k 输出 token/s" if zh else "1×: 32k input tokens/s\n       8k output tokens/s") if mobile else ("1× = 32k 输入 token/s + 8k 输出 token/s；按同倍同步扩展" if zh else "1× = 32k input tokens/s + 8k output tokens/s; both scale together")
    text(page, margin if mobile else 36, 867 if mobile else 691, basis, fs if mobile else 22, MUTED)
    foot = ("两项容量须在同负载下同时达标。\n图例为每月输入与输出 token。\n基线需求满足后，上限趋平。" if zh else "Paired rates must meet the SLO\ntogether, on the same workload.\nVolumes: tokens; B = billion.") if mobile else ("两项容量须在同负载下同时达标；不是独立峰值。图例用量为 token。成本、残值与折现率固定。" if zh else "Both rates must meet the SLO together, not be independent peaks. Volumes in tokens; B = billion.")
    text(page, margin, 978 if mobile else 746, foot, fs if mobile else 21, MUTED)
    save(fig, out, "investment-ceiling", lang, mobile,
         "Conditional 36-month break-even initial investment at monthly token workloads of 24B input/6B output, 48B input/12B output and 64B input/16B output, where B means one billion tokens. The horizontal axis scales a paired SLO-qualified capacity: 1x is 32,000 input tokens/second plus 8,000 output tokens/second, attained together for the same workload; both rates scale by the same multiplier from 0x to 2.5x. At 1x the limits are approximately RMB 0.824m, 1.457m and 1.766m. The amber line is the assumed RMB 1.6m purchase budget. Costs, 10 percent annual discount and fixed RMB 240,000 terminal resale stay constant. The base curve plateaus when online demand is served. This is not a hardware ranking or performance benchmark.")


def cover(data, out):
    fig, page = canvas(1200, 630)
    text(page, 64, 42, "DeepSeek worked example", 20, MUTED)
    text(page, 64, 102, "OPEN MODEL ECONOMICS", 47)
    text(page, 64, 169, "ROI · PAYBACK · NPV", 47)
    ax = fig.add_axes([.055, .17, .89, .37], facecolor=PAPER)
    x = np.array(data["capacityMultipliers"])
    ax.fill_between(x, np.array(data["ceiling"][-1]["values"])/1e6, 0, color=MINT, alpha=.06)
    for row, color, style in zip(data["ceiling"], [BLUE, SLATE, MINT], ["-", "--", "-"]):
        ax.plot(x, np.array(row["values"])/1e6, color=color, lw=4.2, ls=style)
    ax.axhline(data["defaults"]["capex"]/1e6, color=AMBER, lw=2.5, ls=(0, (5, 4)))
    ax.set_xlim(0, 2.5)
    ax.set_ylim(-.75, 3.5)
    ax.set_axis_off()
    text(page, 64, 551, "INPUT + OUTPUT TOKENS", 22, BLUE)
    text(page, 478, 551, "PAIRED TOKEN CAPACITY", 22, MINT)
    text(page, 945, 551, "BUDGET", 22, AMBER)
    fig.savefig(out / "cover.png", dpi=72, metadata={"Software": "FENG · Matplotlib", "Description": "Open model economics: ROI, payback and NPV. DeepSeek is a worked example. Input/output token demand and paired SLO-qualified token capacity determine a conditional investment limit. Curves are illustrative calculations, not measured hardware performance."})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "content/assets/deepseek-roi")
    args = parser.parse_args()
    data = samples()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for lang in ["en", "zh"]:
        for mobile in [False, True]:
            bridge_figure(data, args.output_dir, lang, mobile)
            cashflow_figure(data, args.output_dir, lang, mobile)
            investment_figure(data, args.output_dir, lang, mobile)
    cover(data, args.output_dir)
    print(f"Rendered 12 SVG figures and cover.png to {args.output_dir}")
    print(f"Inspection PNGs: {PREVIEW}")


if __name__ == "__main__":
    main()
