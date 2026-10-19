#!/usr/bin/env python3
"""Reproduce the local-LLM article's illustrative capacity and bandwidth figures.

Run: python3 scripts/render-local-llm-budget.py [--output-dir DIRECTORY]
Requires Matplotlib and Noto Sans CJK SC. No measured benchmark is generated.

Units are deliberate: parameter storage and bandwidth use decimal GB;
resident-memory budgets use binary GiB. A 70-billion-parameter dense model
at exactly four bits per parameter is an arithmetic example, not a specific
checkpoint. Metadata (10%) and runtime + OS reserve (16 GiB) are assumptions.
The full-context FP16 KV example has 80 layers, 8 KV heads, and head size 128.
Real backends, quantization formats and model architectures need their own
measurements; installed unified memory is not an available-GPU-memory limit.

The decode roofline assumes every generated token streams all 35 GB of raw
weights at advertised peak memory bandwidth, ignoring all other traffic and
work. It is an idealized upper bound within that model, not expected speed,
and does not cover MoE, batching, speculative decoding or weight reuse.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "local-llm-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle

GIB = 2**30
PARAMETERS = 70_000_000_000
BITS = 4
METADATA_FRACTION = .1
LAYERS, KV_HEADS, HEAD_DIM, KV_BYTES = 80, 8, 128, 2
RESERVE_GIB = 16
PAPER, INK, MUTED = "#111d24", "#f3f1e8", "#bac8ce"
GRID, BLUE, GOLD, MINT, SLATE = "#354c57", "#91cbe6", "#efbd7f", "#a8d4bb", "#708995"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Noto Sans CJK SC", "DejaVu Sans"],
    "svg.fonttype": "path", "svg.hashsalt": "feng-local-llm-budget",
    "axes.unicode_minus": False,
})


def budget(tokens: int, concurrency: int) -> dict:
    if not isinstance(tokens, int) or not isinstance(concurrency, int) or tokens <= 0 or concurrency <= 0:
        raise ValueError("Positive integer token count and concurrency required")
    weights = PARAMETERS * BITS // 8
    metadata = int(weights * METADATA_FRACTION)
    kv = 2 * LAYERS * KV_HEADS * HEAD_DIM * KV_BYTES * tokens * concurrency
    reserve = RESERVE_GIB * GIB
    parts = {"weights": weights, "quantization_metadata": metadata, "kv_cache": kv, "runtime_and_os_reserve": reserve}
    return {"tokens_per_sequence": tokens, "concurrent_sequences": concurrency,
            "bytes": parts, "GiB": {key: val/GIB for key, val in parts.items()},
            "total_bytes": sum(parts.values()), "total_GiB": sum(parts.values()) / GIB}


def all_data() -> dict:
    rows = [budget(tokens, n) for tokens, n in [(8192, 1), (32768, 1), (8192, 4), (32768, 4)]]
    roofs = [{"memory_bandwidth_GB_s": b, "weights_read_GB_per_token": 35,
              "ideal_tokens_s": b/35} for b in [273, 1200]]
    return {
        "kind": "illustrative arithmetic, not a benchmark or measured hardware memory limit",
        "units": {"GB": 10**9, "GiB": GIB},
        "model_assumptions": {"parameters": PARAMETERS, "weight_bits": BITS,
            "metadata_fraction_of_raw_weights": METADATA_FRACTION,
            "layers": LAYERS, "kv_heads": KV_HEADS, "head_dimension": HEAD_DIM,
            "kv_element_bytes": KV_BYTES, "kv_bytes_per_token_per_sequence": 2*LAYERS*KV_HEADS*HEAD_DIM*KV_BYTES,
            "runtime_and_os_reserve_GiB": RESERVE_GIB},
        "memory_cases": rows,
        "decode_bandwidth_rooflines": roofs,
        "decode_assumptions": ["Dense model, one sequence, one token per decoding step",
            "Every step streams 35 GB of raw weights at the stated peak bandwidth",
            "Ignores metadata, KV traffic, compute, synchronization and operating-system traffic",
            "Not a measured or expected token rate; not applicable to MoE, batching or speculative decoding"],
    }


def canvas(width, height):
    fig = plt.figure(figsize=(width/72, height/72), facecolor=PAPER)
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(height, 0))
    ax.set_axis_off()
    return fig, ax


def text(ax, x, y, s, size=20, color=INK, **kwargs):
    return ax.text(x, y, s, fontsize=size, color=color, va="top", linespacing=1.32, **kwargs)


def line(ax, x1, y1, x2, y2, color=GRID, width=1, **kwargs):
    ax.plot([x1, x2], [y1, y2], color=color, linewidth=width, **kwargs)


def block(ax, x, y, width, height, color, **kwargs):
    ax.add_patch(Rectangle((x, y), width, height, facecolor=color, edgecolor=PAPER, linewidth=1, **kwargs))


def save(fig, out, stem, lang, mobile, desc):
    suffix = "-mobile" if mobile else ""
    name = f"{stem}{suffix}-{lang}"
    svg_path = out / f"{name}.svg"
    fig.savefig(svg_path, metadata={"Date": None, "Description": desc})
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8")
    preview = Path(tempfile.gettempdir()) / "local-llm-figure-review"
    preview.mkdir(exist_ok=True)
    fig.savefig(preview / f"{name}.png", dpi=110)
    plt.close(fig)


def memory_figure(out, lang, mobile=False):
    zh = lang == "zh"
    w, h = (432, 960) if mobile else (864, 560)
    fig, ax = canvas(w, h)
    fs = 22 if mobile else 18
    text(ax, 28 if mobile else 36, 24, "70B，究竟占多少内存？" if zh else ("How much memory\ndoes 70B need?" if mobile else "How much memory does a 70B model need?"), 27 if mobile else 29)
    sub_y = 74 if zh or not mobile else 105
    sub = ("4-bit 稠密权重 · FP16 KV\n预算算例 · 单位 GiB" if zh else "4-bit dense · FP16 KV\nIllustrative budget / GiB") if mobile else ("4-bit 稠密权重 · FP16 KV · 预算算例，不是实测" if zh else "4-bit dense weights · FP16 KV · illustrative budget, not measured")
    text(ax, 28 if mobile else 36, sub_y, sub, fs, MUTED)
    legend_y = 172 if mobile else 120
    legends = [
        ("权重 32.6" if zh else "Weights 32.6", BLUE),
        ("元数据 3.3" if zh else "Metadata 3.3", GOLD),
        ("KV 随负载变化" if zh else "KV varies", MINT),
        ("预留 16" if zh else "Reserve 16", SLATE),
    ]
    for i, (label, color) in enumerate(legends):
        x = (28 + (i%2)*202) if mobile else (36+(i%2)*365)
        y = legend_y + (i//2)*(40 if mobile else 32)
        block(ax, x, y+4, 13, 18, color)
        text(ax, x+22, y, label, fs, MUTED)
    if not mobile:
        text(ax, 824, legend_y+34, "GiB", fs, MUTED, ha="right")
    rows = all_data()["memory_cases"]
    x0, bw = (28, 376) if mobile else (218, 532)
    y0, gap, bh = (310, 132, 36) if mobile else (230, 66, 29)
    keys = ["weights", "quantization_metadata", "kv_cache", "runtime_and_os_reserve"]
    colors = [BLUE, GOLD, MINT, SLATE]
    for i, row in enumerate(rows):
        y = y0+i*gap
        tokens, n, total = row["tokens_per_sequence"], row["concurrent_sequences"], row["total_GiB"]
        if mobile:
            text(ax, 28, y-42, f"{tokens:,} × {n}", 24)
            text(ax, 404, y-42, f"{total:.1f} GiB", 25, INK, ha="right")
            text(ax, 28, y+46, (f"KV {row['GiB']['kv_cache']:g} GiB · {n} 路序列" if zh else f"KV {row['GiB']['kv_cache']:g} GiB · {n} sequence" + ("s" if n>1 else "")), fs, MUTED)
        else:
            text(ax, 36, y-3, f"{tokens:,} × {n}", 21)
            text(ax, 36, y+23, f"KV {row['GiB']['kv_cache']:g} GiB", 16, MUTED)
            text(ax, 824, y-1, f"{total:.1f}", 24, INK, ha="right")
        left = 0
        for key, color in zip(keys, colors):
            val = row["GiB"][key]
            block(ax, x0+left/100*bw, y, val/100*bw, bh, color)
            left += val
    axis_y = 853 if mobile else 471
    for tick in [0, 50, 100]:
        x = x0+tick/100*bw
        line(ax, x, y0-8, x, axis_y-8, width=.7, alpha=.45)
        text(ax, x, axis_y, str(tick), fs if mobile else 17, MUTED, ha="center")
    line(ax, x0, axis_y-8, x0+bw, axis_y-8)
    foot = ("每路缓存 token 数 × 并发数\n预留含运行时与操作系统" if zh else "Cached tokens × sequences\nReserve: runtime + OS") if mobile else ("每路缓存 token 数 × 并发数；预留含运行时与操作系统，不代表 GPU 可用内存上限。" if zh else "Cached tokens × concurrent sequences; reserve includes runtime and OS. Not a GPU memory limit.")
    text(ax, 28 if mobile else 36, 901 if mobile else 512, foot, fs if mobile else 16, MUTED)
    save(fig, out, "memory-budget", lang, mobile,
         "Illustrative 70B four-bit dense memory budget. Weights 32.5963 GiB, assumed metadata 3.2596 GiB, runtime and OS reserve 16 GiB. FP16 KV: 80 layers, 8 KV heads, 128 head dimensions. Cases: 8192x1 = 54.3559 GiB, 32768x1 = 61.8559 GiB, 8192x4 = 61.8559 GiB, 32768x4 = 91.8559 GiB. Not measured hardware usage or a GPU memory limit.")


def decode_figure(out, lang, mobile=False):
    zh = lang == "zh"
    w, h = (432, 850) if mobile else (864, 540)
    fig, ax = canvas(w, h)
    fs = 22 if mobile else 18
    text(ax, 28 if mobile else 36, 24, ("能装下，不等于生成快" if zh else "Fitting is not\nthe same as speed") if mobile else ("装得下之后，还要看每个 token 读多少数据" if zh else "After it fits: how much must each token read?"), 27 if mobile else 28)
    text(ax, 28 if mobile else 36, (104 if mobile and not zh else 75), ("稠密模型 · 单路 · 逐 token 解码" if zh else "Dense · one stream · decode"), fs, MUTED)
    box_y = 154 if mobile else 122
    ax.add_patch(FancyBboxPatch((28 if mobile else 36, box_y), w-56 if mobile else w-72, 110 if mobile else 80,
                                boxstyle="round,pad=0,rounding_size=10", facecolor="#1c3039", edgecolor=GRID))
    text(ax, 46 if mobile else 56, box_y+15, "每 token 读取 35 GB 权重" if zh else ("Read 35 GB / token" if mobile else "Read 35 GB of weights / token"), fs if mobile else 21)
    text(ax, 46 if mobile else 56, box_y+(58 if mobile else 45), "理想上界 = 带宽 ÷ 35 GB" if zh else ("Ceiling = bandwidth / 35 GB" if mobile else "Ideal ceiling = bandwidth / 35 GB"), fs if mobile else 19, MINT)
    x0, bw = (28, 376) if mobile else (245, 480)
    y0, gap, bh = (342, 132, 36) if mobile else (259, 96, 34)
    for i, roof in enumerate(all_data()["decode_bandwidth_rooflines"]):
        y = y0+i*gap
        b, rate = roof["memory_bandwidth_GB_s"], roof["ideal_tokens_s"]
        color = BLUE if i==0 else MINT
        if mobile:
            text(ax, 28, y-45, f"{b:,} GB/s", 25, color)
            text(ax, 404, y-45, f"{rate:.1f} token/s", 25, INK, ha="right")
        else:
            text(ax, 36, y-2, f"{b:,} GB/s", 25, color)
            text(ax, 824, y-2, f"{rate:.1f}", 28, INK, ha="right")
        block(ax, x0, y, rate/40*bw, bh, color)
    ay = 548 if mobile else 406
    for tick in [0,10,20,30,40]:
        x=x0+tick/40*bw
        line(ax,x,y0-7,x,ay-8,width=.7,alpha=.45)
        text(ax,x,ay,str(tick),fs if mobile else 17,MUTED,ha="center")
    line(ax,x0,ay-8,x0+bw,ay-8)
    text(ax,w-28 if mobile else 824,ay+35 if mobile else 407,"token/s",fs if mobile else 17,MUTED,ha="right")
    if mobile:
        text(ax, 28, 627, "带宽上界，不是速度实测" if zh else "Bandwidth ceilings,\nnot measured speeds", 23, GOLD)
        text(ax, 28, 702, "忽略 KV、计算及其他开销\n不适用于 MoE、批处理\n或推测解码" if zh else "Ignores KV and compute.\nExcludes MoE, batching,\nand speculative decoding.", fs, MUTED)
    else:
        text(ax, 36, 465, "理想带宽上界，不是速度实测；忽略 KV、计算及其他开销。" if zh else "Ideal bandwidth ceilings, not measured speeds; KV traffic and compute are omitted.", 18, GOLD)
        text(ax, 36, 500, "不适用于 MoE、批处理或推测解码。" if zh else "This single-stream dense model excludes MoE, batching and speculative decoding.", 17, MUTED)
    save(fig, out, "decode-bound", lang, mobile,
         "Illustrative ideal bandwidth ceilings: 273 GB/s divided by 35 GB/token = 7.8 token/s; 1200 GB/s divided by 35 GB/token = 34.2857 token/s. One-stream dense decoding, all raw weights streamed each token. Not a benchmark or expected speed. Ignores all other work and traffic; not applicable to MoE, batching or speculative decoding.")


def cover(out):
    fig, ax = canvas(1200, 630)
    text(ax, 70, 55, "LOCAL AI / SYSTEMS", 20, MUTED)
    text(ax, 70, 115, "FIT THE MODEL.", 49)
    text(ax, 70, 180, "FEED THE COMPUTE.", 49)
    text(ax, 70, 267, "CAPACITY  /  BANDWIDTH  /  SOFTWARE", 19, MUTED)
    # A memory array and compute tile: an abstract data path, not product hardware.
    for i in range(9):
        for j in range(3):
            block(ax, 70+i*38, 365+j*42, 30, 34, BLUE if i<6 else MINT, alpha=.95)
    for j in range(3):
        line(ax, 412, 382+j*42, 775, 382+j*42, MINT if j==1 else SLATE, 3)
        for i in range(4):
            x=465+i*75+(12 if j==1 else 0)
            ax.add_patch(Circle((x,382+j*42),4,facecolor=INK,edgecolor="none"))
    ax.add_patch(FancyBboxPatch((790,339),225,176,boxstyle="round,pad=0,rounding_size=18",facecolor="#1d303a",edgecolor=MINT,linewidth=2))
    for i in range(4):
        for j in range(3):
            block(ax, 815+i*46, 364+j*43, 33, 30, MINT if (i+j)%3 else BLUE)
    text(ax, 70, 543, "MEMORY", 21, MUTED)
    text(ax, 790, 543, "COMPUTE", 21, MUTED)
    fig.savefig(out / "cover.png", dpi=72, metadata={"Software":"FENG · Matplotlib","Description":"Abstract local-LLM capacity, data movement and compute; not a product rendering."})
    plt.close(fig)


def check(data):
    assert data["model_assumptions"]["kv_bytes_per_token_per_sequence"] == 327680
    rows=data["memory_cases"]
    assert [r["GiB"]["kv_cache"] for r in rows] == [2.5, 10, 10, 40]
    assert rows[1]["total_bytes"] == rows[2]["total_bytes"]
    assert all(r["total_bytes"] == sum(r["bytes"].values()) for r in rows)
    assert math.isclose(rows[0]["total_GiB"], 54.35591912269592)
    assert [r["ideal_tokens_s"] for r in data["decode_bandwidth_rooflines"]] == [7.8, 1200/35]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=Path("content/assets/local-llm"))
    args=parser.parse_args()
    data=all_data()
    check(data)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    for lang in ["en","zh"]:
        for mobile in [False,True]:
            memory_figure(args.output_dir,lang,mobile)
            decode_figure(args.output_dir,lang,mobile)
    cover(args.output_dir)
    (args.output_dir/"data.json").write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    for path in args.output_dir.glob("*.svg"):
        ET.parse(path)
    print(json.dumps({"status":"ok","files":len(list(args.output_dir.iterdir())),"totals_GiB":[r["total_GiB"] for r in data["memory_cases"]],"previews":"/tmp/local-llm-figure-review"},ensure_ascii=False))


if __name__ == "__main__":
    main()
