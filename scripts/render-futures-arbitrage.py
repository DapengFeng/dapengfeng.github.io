#!/usr/bin/env python3
"""Render hypothetical futures teaching figures with Matplotlib.

Run from any directory with Python 3, Matplotlib and Noto Sans CJK JP installed.
Site builds consume the exported assets and do not need plotting dependencies.
The standard-library calculator in the asset directory owns the numeric inputs.
"""
from pathlib import Path
import json
import os
import runpy

os.environ.setdefault('MPLCONFIGDIR', '/tmp/feng-futures-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'content/assets/futures-arbitrage'
DATA = runpy.run_path(str(OUT / 'calculate-examples.py'))['calculate']()
BG, INK, TEAL, RED = '#f4f0e7', '#1d3442', '#176c60', '#b64d3d'
MUTED, GRID = '#617078', '#d7d9d1'
plt.rcParams.update({
    'font.family': 'Noto Sans CJK JP', 'font.size': 18,
    'svg.fonttype': 'path', 'svg.hashsalt': 'feng-futures-arbitrage-2026',
    'axes.unicode_minus': False, 'figure.facecolor': BG,
    'axes.facecolor': BG, 'savefig.facecolor': BG,
})


def canvas(width, height):
    return plt.figure(figsize=(width / 100, height / 100), dpi=100)


def text(fig, x, y, value, size=24, color=INK, weight='normal', ha='left'):
    width, height = fig.get_size_inches() * 100
    return fig.text(x / width, 1 - y / height, value, fontsize=size * .72,
                    color=color, weight=weight, ha=ha, va='top', linespacing=1.45)


def line(fig, x0, x1, y, color=GRID, width=1):
    w, h = fig.get_size_inches() * 100
    fig.add_artist(Line2D([x0 / w, x1 / w], [1 - y / h] * 2,
                         transform=fig.transFigure, color=color, linewidth=width))


def axes(fig, x, y, width, height):
    fw, fh = fig.get_size_inches() * 100
    ax = fig.add_axes([x / fw, 1 - (y + height) / fh, width / fw, height / fh])
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.spines['bottom'].set_color(GRID)
    ax.tick_params(colors=MUTED, length=0, pad=10, labelsize=24 * .72)
    ax.grid(axis='y', color=GRID, linewidth=.8)
    ax.set_axisbelow(True)
    return ax


def legend(fig, x, y, label, color, dashed=False, size=24):
    w, h = fig.get_size_inches() * 100
    fig.add_artist(Line2D([x / w, (x + 38) / w], [1 - (y + 13) / h] * 2,
                         transform=fig.transFigure, color=color, linewidth=2.6,
                         linestyle=(0, (5, 3)) if dashed else '-'))
    text(fig, x + 52, y, label, size, color)


def save(fig, stem, png=False):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Fail early if a manual label extends outside its exported artwork.
    for t in list(fig.texts) + [t for a in fig.axes for t in a.texts]:
        if not t.get_visible():
            continue
        box = t.get_window_extent(renderer)
        assert box.x0 >= -1 and box.x1 <= fig.bbox.width + 1, (stem, t.get_text(), box)
        assert box.y0 >= -1 and box.y1 <= fig.bbox.height + 1, (stem, t.get_text(), box)
    path = OUT / f'{stem}.svg'
    fig.savefig(path, metadata={'Date': None})
    path.write_text('\n'.join(s.rstrip() for s in path.read_text().splitlines()) + '\n')
    if png:
        fig.savefig(OUT / f'{stem}.png', dpi=100)
    # PNGs for inspection live outside the public asset directory.
    fig.savefig(Path('/tmp') / f'futures-{stem}.png', dpi=100)
    plt.close(fig)


def liquidity(lang, mobile):
    zh = lang == 'zh'
    rows = DATA['cash_and_carry']['interim_path']
    fig = canvas(480, 1030) if mobile else canvas(1200, 800)
    if mobile:
        text(fig, 30, 30, '01 / ' + ('资金约束' if zh else 'FUNDING CONSTRAINT'), 20, TEAL, 'bold')
        text(fig, 30, 76, '库存有浮盈，\n账户仍会缺现金' if zh else 'Inventory gains,\ncash runs short', 34, INK, 'bold')
        text(fig, 30, 190, '金额：元 · 构造的结算路径' if zh else 'CNY · Assumed settlements', 23, MUTED)
        legend(fig, 30, 236, '账户权益' if zh else 'Account equity', TEAL)
        legend(fig, 30, 278, '所需保证金' if zh else 'Required margin', RED, True)
        ax = axes(fig, 82, 345, 338, 335)
    else:
        text(fig, 64, 35, '01 / ' + ('资金约束' if zh else 'FUNDING CONSTRAINT'), 20, TEAL, 'bold')
        text(fig, 64, 77, '库存有浮盈，账户仍会缺现金' if zh else 'Inventory gains; the futures account runs short', 38, INK, 'bold')
        text(fig, 64, 135, '金额：元 · 构造的结算路径' if zh else 'CNY · Assumed settlement path', 23, MUTED)
        legend(fig, 580, 143, '账户权益' if zh else 'Account equity', TEAL)
        legend(fig, 870, 143, '所需保证金' if zh else 'Required margin', RED, True)
        ax = axes(fig, 112, 220, 930, 302)

    xx = [0, 1, 2]
    equity = [r['futures_account_equity'] for r in rows]
    margin = [r['required_margin'] for r in rows]
    ax.plot(xx, equity, color=TEAL, lw=3, marker='o', markersize=6, zorder=4)
    ax.plot(xx, margin, color=RED, lw=2.5, linestyle=(0, (5, 3)), marker='s', markersize=5, zorder=4)
    ax.set_xlim(-.12, 2.17 if mobile else 2.28)
    ax.set_ylim(1400, 3300)
    ax.set_xticks(xx, ['T0', 'T1', 'T2'])
    ax.set_yticks([1600, 2200, 2800, 3200], ['1,600', '2,200', '2,800', '3,200'])
    crossing = 1 + 230 / (230 + 540)
    crossing_value = 2360 - 700 * (crossing - 1)
    ax.fill_between([crossing, 2], [crossing_value, 1660], [crossing_value, 2200], color=RED, alpha=.14, zorder=2)
    ax.vlines(2, 1660, 2200, color=RED, linewidth=1.6, zorder=3)
    # Values at points are offset away from the lines and shared y-axis.
    points = [(0, 3060, '3,060', TEAL, (9, 10)),
              (1, 2360, '2,360', TEAL, (-10, 13)),
              (2, 1660, '1,660', TEAL, (-12, -30)),
              (0, 2060, '2,060', RED, (9, -29)),
              (1, 2130, '2,130', RED, (-10, -31)),
              (2, 2200, '2,200', RED, (-12, 12))]
    for x, y, value, color, offset in points:
        ax.annotate(value, (x, y), xytext=offset, textcoords='offset pixels',
                    fontsize=(23 if mobile else 25) * .72, color=color,
                    ha='left' if x == 0 else 'center', va='bottom')

    if mobile:
        line(fig, 30, 450, 745)
        text(fig, 30, 766, '可用资金' if zh else 'Available funds', 24, MUTED)
        for x, value in zip([30, 176, 330], ['+1,000', '+230', '−540']):
            text(fig, x, 810, value, 31, RED if value == '−540' else INK, 'bold')
        text(fig, 30, 865, 'T0', 21, MUTED); text(fig, 176, 865, 'T1', 21, MUTED); text(fig, 330, 865, 'T2', 21, MUTED)
        line(fig, 30, 450, 910)
        text(fig, 30, 930, '现货浮盈 +1,400 元' if zh else 'Unrealized spot gain: +1,400', 26 if zh else 24, TEAL, 'bold')
        text(fig, 30, 975, '浮盈不能自动补入期货账户。' if zh else 'Not cash in the futures account.', 23, MUTED)
    else:
        line(fig, 64, 1136, 594)
        text(fig, 64, 624, '可用资金' if zh else 'Available funds', 25, MUTED)
        for x, t, value in [(343, 'T0', '+1,000'), (606, 'T1', '+230'), (893, 'T2', '−540')]:
            text(fig, x, 619, value, 37, RED if value == '−540' else INK, 'bold')
            text(fig, x, 673, t, 21, MUTED)
        line(fig, 64, 1136, 721)
        text(fig, 64, 744, '现货浮盈 +1,400 元' if zh else 'Unrealized spot gain: +1,400 CNY', 25, TEAL, 'bold')
        text(fig, 460 if zh else 570, 744, '浮盈不能自动补入期货账户。' if zh else 'Not cash in the futures account.', 25, MUTED)
    save(fig, f'liquidity-{lang}' + ('-mobile' if mobile else ''))


def calendar(lang, mobile):
    zh = lang == 'zh'
    fig = canvas(480, 1040) if mobile else canvas(1200, 800)
    if mobile:
        text(fig, 30, 30, '02 / ' + ('跨期价差' if zh else 'CALENDAR SPREAD'), 20, TEAL, 'bold')
        text(fig, 30, 76, '两腿都涨，\n价差仍可缩小' if zh else 'Both rise.\nThe spread narrows.', 33, INK, 'bold')
        text(fig, 30, 187, '价格：元/吨 · 每腿 10 吨' if zh else 'CNY/t · 10 tonnes per leg', 23, MUTED)
        legend(fig, 30, 235, '近月：买入' if zh else 'Near: long', TEAL)
        legend(fig, 30, 277, '远月：卖出' if zh else 'Far: short', INK, True)
        ax = axes(fig, 88, 346, 322, 325)
    else:
        text(fig, 64, 35, '02 / ' + ('跨期价差' if zh else 'CALENDAR SPREAD'), 20, TEAL, 'bold')
        text(fig, 64, 77, '两腿都涨，价差仍可缩小' if zh else 'Both rise. The spread narrows.', 38, INK, 'bold')
        text(fig, 64, 135, '价格：元/吨 · 每腿 10 吨' if zh else 'CNY/t · 10 tonnes per leg', 23, MUTED)
        legend(fig, 620, 143, '近月：买入' if zh else 'Near: long', TEAL)
        legend(fig, 900, 143, '远月：卖出' if zh else 'Far: short', INK, True)
        ax = axes(fig, 120, 217, 900, 315)
    ax.set_xlim(-.15, 1.15)
    ax.set_ylim(3980, 4140)
    ax.set_yticks([4000, 4050, 4100], ['4,000', '4,050', '4,100'])
    ax.set_xticks([0, 1], ['开仓' if zh else 'Entry', '平仓' if zh else 'Exit'])
    near = [4000, 4050]; far = [4100, 4110]
    ax.fill_between([0, 1], near, far, color=TEAL, alpha=.09, zorder=2)
    ax.plot([0, 1], near, color=TEAL, lw=3, marker='o', markersize=7, zorder=4)
    ax.plot([0, 1], far, color=INK, lw=2.7, linestyle=(0, (5, 3)), marker='s', markersize=6, zorder=4)
    for x, gap, middle in [(0, 100, 4050), (1, 60, 4080)]:
        ax.vlines(x, near[x], far[x], color=TEAL, linewidth=1.4, zorder=3)
        ax.text(x + (.085 if x == 0 else -.085), middle, f'D = {gap}',
                color=TEAL, fontsize=(24 if mobile else 28) * .72,
                ha='left' if x == 0 else 'right', va='center', weight='bold',
                bbox={'facecolor': BG, 'edgecolor': 'none', 'pad': 3, 'alpha': .94})
    for x, y, value, color, offset in [
        (0, 4100, '4,100', INK, (0, 11)), (1, 4110, '4,110', INK, (0, 11)),
        (0, 4000, '4,000', TEAL, (0, -30)), (1, 4050, '4,050', TEAL, (0, -30))]:
        ax.annotate(value, (x, y), xytext=offset, textcoords='offset pixels',
                    fontsize=(23 if mobile else 26) * .72, color=color, ha='center', va='bottom')
    if mobile:
        text(fig, 30, 738, '近月 +500' if zh else 'Near +500', 27, TEAL, 'bold')
        text(fig, 260, 738, '远月 −100' if zh else 'Far −100', 27, RED, 'bold')
        line(fig, 30, 450, 794)
        text(fig, 30, 820, '毛利' if zh else 'Gross', 24, MUTED)
        text(fig, 450, 810, '+400', 39, INK, 'bold', ha='right')
        text(fig, 30, 880, '假设总成本' if zh else 'Assumed total cost', 24, MUTED)
        text(fig, 450, 880, '−40', 27, RED, ha='right')
        line(fig, 30, 450, 933)
        text(fig, 30, 960, '净收益 · 元' if zh else 'Net · CNY', 24, MUTED)
        text(fig, 450, 950, '+360', 39, TEAL, 'bold', ha='right')
    else:
        line(fig, 64, 1136, 600)
        for x, label, value, color in [
            (64, '近月' if zh else 'Near leg', '+500', TEAL),
            (326, '远月' if zh else 'Far leg', '−100', RED),
            (588, '毛利' if zh else 'Gross', '+400', INK),
            (880, '净收益' if zh else 'Net', '+360', TEAL)]:
            text(fig, x, 623, label, 24, MUTED)
            text(fig, x, 668, value, 41, color, 'bold')
        text(fig, 64, 748, '两次观察；连线仅作示意。' if zh else 'Two observations; connecting lines are schematic.', 23, MUTED)
        text(fig, 1136, 748, '假设总成本 40 元' if zh else 'Assumed total cost: 40 CNY', 23, MUTED, ha='right')
    save(fig, f'calendar-{lang}' + ('-mobile' if mobile else ''))


def cover():
    fig = canvas(1200, 630)
    text(fig, 64, 46, 'FENG / FINANCE & DECISIONS', 20, TEAL, 'bold')
    text(fig, 64, 109, 'FUTURES &\nARBITRAGE', 66, INK, 'bold')
    text(fig, 64, 317, '期货交易规则与套利模型', 34, INK)
    text(fig, 64, 383, '价差如何兑现，现金如何周转', 25, MUTED)
    ax = fig.add_axes([.54, .2, .4, .57])
    ax.set_xlim(-.1, 1.1); ax.set_ylim(3980, 4140); ax.axis('off')
    ax.fill_between([0, 1], [4000, 4050], [4100, 4110], color=TEAL, alpha=.1)
    ax.plot([0, 1], [4000, 4050], color=TEAL, lw=4, marker='o', markersize=8)
    ax.plot([0, 1], [4100, 4110], color=INK, lw=3, linestyle=(0, (5, 3)), marker='s', markersize=7)
    ax.vlines([0, 1], [4000, 4050], [4100, 4110], color=TEAL, lw=2)
    ax.text(.06, 4050, '100', color=TEAL, fontsize=32, weight='bold', va='center')
    ax.text(.94, 4080, '60', color=TEAL, fontsize=32, weight='bold', va='center', ha='right')
    line(fig, 64, 1136, 534)
    text(fig, 64, 558, 'Rules · Cash flows · Spreads', 24, MUTED)
    text(fig, 1136, 558, 'Illustrative prices', 20, MUTED, ha='right')
    save(fig, 'cover', png=True)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'example.json').write_text(json.dumps(DATA, ensure_ascii=False, indent=2) + '\n')
    for language in ['en', 'zh']:
        for is_mobile in [False, True]:
            liquidity(language, is_mobile)
            calendar(language, is_mobile)
    cover()
    print('Validated both hypothetical examples; rendered 8 responsive SVG figures and the cover.')
