#!/usr/bin/env python3
"""Reproduce synthetic SMC teaching figures, not investment performance.

Requires matplotlib and Noto Sans CJK JP. Run python3 scripts/render-smc.py.
The checked-in OHLC file is the input; no market data or network is used.
"""
from pathlib import Path
import argparse
import json
import math
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/feng-smc-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'content/assets/smc'
DATA = json.loads((OUT / 'example.json').read_text())
BARS = DATA['bars']
BG, PANEL, GRID = '#101a21', '#15232d', '#2a3c48'
FG, MUTED, UP, DOWN, ZONE = '#edf3f4', '#a5b8c4', '#61d3b8', '#f09383', '#efc57a'
plt.rcParams.update({'font.family': 'Noto Sans CJK JP', 'font.size': 14,
                     'svg.fonttype': 'path', 'svg.hashsalt': 'feng-smc-2026',
                     'axes.unicode_minus': False, 'figure.facecolor': BG,
                     'axes.facecolor': BG, 'savefig.facecolor': BG})


def validate():
    assert len(BARS) == 12
    assert all(b['low'] <= min(b['open'], b['close']) <=
               max(b['open'], b['close']) <= b['high'] for b in BARS)
    k = DATA['pivot_k']
    def confirmed(p, field, high, at):
        if at < p + k:
            return False
        center = BARS[p - 1][field]
        neighbors = BARS[p - k - 1:p - 1] + BARS[p:p + k]
        return all(center > b[field] if high else center < b[field] for b in neighbors)
    assert not confirmed(3, 'low', False, 4)
    assert confirmed(3, 'low', False, 5)
    assert not confirmed(6, 'high', True, 7)
    assert confirmed(6, 'high', True, 8)
    assert BARS[8]['low'] < 98 < BARS[8]['close']
    assert all(b['close'] <= 107 for b in BARS[8:10])
    assert BARS[10]['close'] > 107
    assert BARS[10]['low'] == BARS[8]['high']  # strict FVG test fails at 11
    assert BARS[11]['low'] > BARS[9]['high']
    assert (BARS[9]['high'], BARS[11]['low']) == (102, 105)
    assert BARS[10]['open'] < 102 < 105 < BARS[10]['close']
    last_down = max(b['bar'] for b in BARS[:11] if b['close'] < b['open'])
    assert last_down == 8
    r = DATA['separate_risk_example']
    risk = r['entry'] - r['stop']
    win = (r['target'] - r['entry'] - r['all_in_cost_per_unit']) / risk
    loss = (r['stop'] - r['entry'] - r['all_in_cost_per_unit']) / risk
    breakeven = -loss / (win - loss)
    assert math.isclose(win, 1.95) and math.isclose(loss, -1.05)
    assert math.isclose(breakeven, .35)
    assert math.isclose(.4 * win + .6 * loss, .15)
    assert math.isclose((.5 * win + .5 * loss), .45)
    return {'ohlc_valid': True, 'pivot_confirmation_checked': True,
            'strict_fvg_checked': True, 'retrospective_ob_candidate_checked': True,
            'risk_arithmetic_checked': True, 'market_profitability_tested': False}


def label(fig, x, y, s, size=16, color=FG, weight='normal', ha='left', va='top'):
    return fig.text(x, y, s, fontsize=size, color=color, weight=weight,
                    ha=ha, va=va, linespacing=1.5)


def canvas(w, h):
    return plt.figure(figsize=(w / 100, h / 100), dpi=100)


def save(fig, stem):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    artists = list(fig.texts) + [t for ax in fig.axes for t in ax.texts]
    for text in artists:
        box = text.get_window_extent(renderer)
        assert 0 <= box.x0 and box.x1 <= fig.bbox.width, (stem, text.get_text(), box)
        assert 0 <= box.y0 and box.y1 <= fig.bbox.height, (stem, text.get_text(), box)
    for ext in ('svg', 'png'):
        destination = OUT / f'{stem}.{ext}'
        fig.savefig(destination, dpi=100, metadata={'Date': None} if ext == 'svg' else None)
        if ext == 'svg':
            destination.write_text(
                '\n'.join(line.rstrip() for line in destination.read_text(encoding='utf-8').splitlines()) + '\n',
                encoding='utf-8',
            )
    plt.close(fig)


def candles(ax, bars, width=.56):
    for b in bars:
        c = UP if b['close'] >= b['open'] else DOWN
        ax.vlines(b['bar'], b['low'], b['high'], color=c, linewidth=1.7, zorder=4)
        y = min(b['open'], b['close'])
        ax.add_patch(Rectangle((b['bar'] - width / 2, y), width,
                              max(abs(b['close'] - b['open']), .05),
                              color=c, zorder=5))


def axes_style(ax, ticksize=15):
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.spines['bottom'].set_color(GRID)
    ax.grid(axis='y', color=GRID, linewidth=.65, alpha=.7)
    ax.set_axisbelow(True)
    ax.tick_params(axis='both', colors=MUTED, labelsize=ticksize, length=0, pad=8)


def head(fig, number, title, sub, mobile):
    label(fig, .065, .968, number, 12 if mobile else 13, UP, 'bold')
    label(fig, .065, .925 if mobile else .912, title, 22 if mobile else 27, FG, 'bold')
    label(fig, .065, .853 if mobile else .843, sub, 14 if mobile else 16, MUTED)


def pxlabel(fig, x, y, s, size=19.5, color=FG, weight='normal'):
    width, height = fig.get_size_inches() * 100
    return label(fig, x / width, 1 - y / height, s, size, color, weight)


def pixel_axes(fig, x, y, width, height):
    fw, fh = fig.get_size_inches() * 100
    return fig.add_axes([x / fw, 1 - (y + height) / fh, width / fw, height / fh])


def mobile_head(fig, number, title, sub):
    pxlabel(fig, 32, 30, number, 15, UP, 'bold')
    pxlabel(fig, 32, 78, title, 22, FG, 'bold')
    pxlabel(fig, 32, 137, sub, 19.5, MUTED)


def chronology_mobile(lang):
    zh = lang == 'zh'
    fig = canvas(480, 1460)
    mobile_head(fig, '03 / ' + ('信息的时间顺序' if zh else 'INFORMATION TIMING'),
                '信号在何时才可知？' if zh else 'When is a signal known?',
                '12 根合成 K 线\n所有判定使用收盘信息。' if zh else
                '12 synthetic bars.\nDecisions use the bar close.')
    ax = pixel_axes(fig, 76, 260, 372, 365)
    axes_style(ax, 19)
    candles(ax, BARS, .5)
    ax.set_xlim(.5, 12.6); ax.set_ylim(96, 112)
    ax.set_xticks([1, 3, 5, 7, 9, 11]); ax.set_yticks([98, 102, 107, 110])
    ax.tick_params(axis='x', labelsize=19)
    ax.hlines(98, 5, 12.5, color=ZONE, linestyle=(0, (4, 4)), linewidth=1.4)
    ax.hlines(107, 8, 12.5, color=ZONE, linestyle=(0, (4, 4)), linewidth=1.4)
    ax.scatter([3, 6], [98, 107], s=50, facecolors=BG, edgecolors=ZONE, linewidths=1.5, zorder=7)
    ax.scatter([5, 8], [98, 107], s=20, color=ZONE, zorder=7)
    pxlabel(fig, 140, 688, 'K 线序号 · 1–12' if zh else 'Bar number · 1–12', 19, MUTED)
    rows = [
        ('05 / 08', '高低点确认' if zh else 'Levels confirmed',
         '低点 3 = 98，高点 6 = 107。\n各等待后面两根 K 线。' if zh else
         'Low 3 = 98; high 6 = 107.\nWait two later bars for each.'),
        ('09', '下探后收回' if zh else 'Sweep + reclaim',
         '最低价 97 < 98。\n收盘价 100.5 > 98。' if zh else 'Low 97 < 98.\nClose 100.5 > 98.'),
        ('11', '向上结构突破' if zh else 'Upward break',
         '收盘价 108.5 > 107。\n第 8 根此时成为 OB 候选。' if zh else
         'Close 108.5 > 107.\nBar 8 is now an OB candidate.'),
        ('12', 'FVG 确认' if zh else 'FVG confirmed',
         '第 12 根低价 > 第 10 根高价。\n区间：[102, 105]。' if zh else
         'Low 12 = 105 > high 10 = 102.\nZone: [102, 105].')]
    for i, (number, name, detail) in enumerate(rows):
        y = 755 + 145 * i
        pxlabel(fig, 32, y, number, 19.5, ZONE, 'bold')
        pxlabel(fig, 155, y, name, 20, FG, 'bold')
        pxlabel(fig, 32, y + 47, detail, 19.5, MUTED)
    pxlabel(fig, 32, 1360, '虚线从确认时刻开始。\n合成示例，不是收益记录。' if zh else
            'Levels begin when confirmed.\nSynthetic prices, not returns.', 19.5, MUTED)
    save(fig, f'chronology-{lang}-mobile')


def fvg_mobile(lang):
    zh = lang == 'zh'
    fig = canvas(480, 1220)
    mobile_head(fig, '04 / ' + ('区间的含义' if zh else 'WHAT THE ZONE MEANS'),
                'FVG 不是“没有成交”' if zh else 'An FVG is not “no trading”',
                '三根 K 线的几何条件，\n不是订单簿证据。' if zh else
                'Three-bar geometry.\nNot order-book evidence.')
    ax = pixel_axes(fig, 76, 245, 372, 390); axes_style(ax, 19)
    ax.set_xlim(9.35, 12.8); ax.set_ylim(98.5, 111.4)
    ax.set_xticks([10, 11, 12]); ax.set_yticks([100, 102, 105, 108, 110])
    ax.axhspan(102, 105, facecolor=ZONE, alpha=.12, zorder=1)
    ax.hlines([102, 105], 9.5, 12.65, color=ZONE, linewidth=1.2, linestyle=(0, (4, 4)))
    candles(ax, BARS[9:12], .5)
    ax.plot([12.48, 12.58, 12.58, 12.48], [102, 102, 105, 105], color=ZONE, linewidth=2)
    pxlabel(fig, 186, 705, 'K 线序号' if zh else 'Bar number', 19, MUTED)
    pxlabel(fig, 32, 765, '[102, 105]', 29, ZONE, 'bold')
    pxlabel(fig, 32, 830, '本例第 11 根实体：\n101.5 到 108.5，跨过 FVG。' if zh else
            'Here, bar 11 opens at 101.5\nand closes at 108.5.', 19.5, FG)
    pxlabel(fig, 32, 935, 'OHLC 不能证明每个价位\n都发生交易，也不能识别\n未成交的机构订单。' if zh else
            'OHLC cannot prove trades at\nevery tick, or the existence\nof unfilled institutional orders.', 19.5, MUTED)
    pxlabel(fig, 32, 1120, '合成价格；回补没有保证。' if zh else
            'Synthetic prices.\nA revisit is not guaranteed.', 19.5, MUTED)
    save(fig, f'fvg-{lang}-mobile')


def risk_mobile(lang):
    zh = lang == 'zh'
    fig = canvas(480, 1320)
    mobile_head(fig, '05 / ' + ('从形态到收益' if zh else 'FROM PATTERN TO PAYOFF'),
                '先算成本，再谈胜率' if zh else 'Costs change the odds',
                '独立算例；每单位\n往返总成本 0.10。' if zh else
                'Separate example. Round-trip\ncost: 0.10 per unit.')
    ax = pixel_axes(fig, 32, 245, 416, 320)
    ax.set_xlim(0, 1); ax.set_ylim(97, 105); ax.axis('off')
    ax.axhspan(100, 104, xmin=.1, xmax=.5, facecolor=UP, alpha=.1)
    ax.axhspan(98, 100, xmin=.1, xmax=.5, facecolor=DOWN, alpha=.1)
    for y, color, name in [(104, UP, '目标' if zh else 'Target'),
                           (100, FG, '入场' if zh else 'Entry'),
                           (98, DOWN, '止损' if zh else 'Stop')]:
        ax.hlines(y, .10, .50, color=color, linewidth=2)
        ax.text(.57, y, f'{name}  {y}', fontsize=20, color=color, va='center')
    ax.text(.13, 102, '+4', color=UP, fontsize=28, weight='bold', va='center')
    ax.text(.13, 99, '−2', color=DOWN, fontsize=28, weight='bold', va='center')
    pairs = [('+1.95R', '获利净值 (4 − 0.10) / 2' if zh else 'Net win (4 − 0.10) / 2', UP),
             ('−1.05R', '亏损净值 (−2 − 0.10) / 2' if zh else 'Net loss (−2 − 0.10) / 2', DOWN),
             ('35%', '盈亏平衡胜率\n1.05 / (1.95 + 1.05)' if zh else
              'Break-even win rate\n1.05 / (1.95 + 1.05)', ZONE)]
    for i, (value, caption, color) in enumerate(pairs):
        y = 630 + i * 150
        pxlabel(fig, 32, y, value, 29, color, 'bold')
        pxlabel(fig, 32, y + 55, caption, 19.5, MUTED)
    pxlabel(fig, 32, 1120, 'R = 2；假定成交价为\n100、98 或 104。' if zh else
            'R = 2. Assumed fills:\n100, 98 or 104.', 19.5, FG)
    pxlabel(fig, 32, 1230, '跳空、滑点、部分成交\n会改变实际结果。' if zh else
            'Gaps, slippage and partial fills\ncan change the outcome.', 18.5, MUTED)
    save(fig, f'risk-{lang}-mobile')


def chronology(lang, mobile):
    if mobile:
        return chronology_mobile(lang)
    zh = lang == 'zh'
    fig = canvas(1200, 820)
    head(fig, '03 / '+('信息的时间顺序' if zh else 'INFORMATION TIMING'),
         '信号在何时才可知？' if zh else 'When is a signal known?',
         '12 根合成 K 线 · 所有判定使用收盘信息' if zh else
         '12 synthetic bars · decisions at bar close', False)
    ax = fig.add_axes([.065, .42, .89, .36]); axes_style(ax)
    candles(ax, BARS)
    ax.set_xlim(.5, 12.6); ax.set_ylim(96, 112)
    ax.set_xticks(range(1, 13)); ax.set_yticks([98, 102, 107, 110])
    ax.hlines(98, 5, 12.5, color=ZONE, linestyle=(0, (4, 4)), linewidth=1.4, alpha=.85)
    ax.hlines(107, 8, 12.5, color=ZONE, linestyle=(0, (4, 4)), linewidth=1.4, alpha=.85)
    ax.scatter([3, 6], [98, 107], s=50, facecolors=BG, edgecolors=ZONE, linewidths=1.5, zorder=7)
    ax.scatter([5, 8], [98, 107], s=20, color=ZONE, zorder=7)
    ax.text(.5, -.20, 'K 线序号' if zh else 'Bar number', transform=ax.transAxes,
            ha='center', color=MUTED, fontsize=14)
    notes = [
        ('05 / 08', '高低点确认' if zh else 'Levels confirmed',
         '低点 3 = 98，高点 6 = 107\n均等待后面两根 K 线。' if zh else
         'Low 3 = 98; high 6 = 107.\nWait two bars to confirm.'),
        ('09', '下探后收回' if zh else 'Sweep + reclaim',
         '最低价 97 < 98\n收盘价 100.5 > 98。' if zh else 'Low 97 < 98.\nClose 100.5 > 98.'),
        ('11', '向上结构突破' if zh else 'Upward break',
         '收盘价 108.5 > 107\n此时才标记第 8 根为 OB 候选。' if zh else
         'Close 108.5 > 107.\nOB candidate: bar 8.'),
        ('12', 'FVG 得到确认' if zh else 'FVG confirmed',
         '105 > 102（第 10 根最高价）\n候选区间为 [102, 105]。' if zh else
         'Low 12 = 105; high 10 = 102.\nThe zone is [102, 105].')]
    for i, (n, t, body) in enumerate(notes):
        x = .065 + .235 * i
        label(fig, x, .296, n, 17, ZONE, 'bold')
        label(fig, x, .245, t, 16, FG, 'bold')
        label(fig, x, .196, body, 14, MUTED)
    label(fig, .065, .035, '虚线从确认时刻开始；合成示例，不是收益记录。' if zh else
          'Dashed levels begin when confirmed. Synthetic example, not a return history.', 13, MUTED)
    save(fig, f'chronology-{lang}')


def fvg(lang, mobile):
    if mobile:
        return fvg_mobile(lang)
    zh = lang == 'zh'
    fig = canvas(1200, 720)
    head(fig, '04 / '+('区间的含义' if zh else 'WHAT THE ZONE MEANS'),
         'FVG 不是“没有成交”' if zh else 'An FVG is not “no trading”',
         '三根 K 线的几何条件，不是订单簿证据。' if zh else
         'A three-bar price pattern, not an observation of the order book.', False)
    ax = fig.add_axes([.065, .20, .50, .57]); axes_style(ax)
    ax.set_xlim(9.35, 12.8); ax.set_ylim(98.5, 111.4)
    ax.set_xticks([10, 11, 12]); ax.set_yticks([100, 102, 105, 108, 110])
    ax.axhspan(102, 105, facecolor=ZONE, alpha=.12, zorder=1)
    ax.hlines([102, 105], 9.5, 12.65, color=ZONE, linewidth=1.2, linestyle=(0, (4, 4)))
    candles(ax, BARS[9:12], .5)
    ax.plot([12.48, 12.58, 12.58, 12.48], [102, 102, 105, 105], color=ZONE, linewidth=2)
    ax.text(.5, -.14, 'K 线序号' if zh else 'Bar number', color=MUTED,
            transform=ax.transAxes, ha='center', fontsize=14)
    label(fig, .625, .74, '[102, 105]', 34, ZONE, 'bold')
    label(fig, .625, .635, '第 12 根最低价 > 第 10 根最高价' if zh else
          'Low 12 > high 10', 18, FG, 'bold')
    label(fig, .625, .552, '本例中间一根实体：101.5 到 108.5。\n这个实体的价格范围跨过了 FVG。' if zh else
          'This example: body 101.5 to 108.5.\nIts range spans the FVG.', 17, FG)
    label(fig, .625, .365, '但 OHLC 不能证明每个价位均有成交，\n也不能识别谁在下单、多少订单未成交。' if zh else
          'OHLC cannot establish trades at every tick,\nwho placed orders, or unfilled volume.', 15, MUTED)
    label(fig, .065, .04, '合成示例；“回补”只是待检验的后续事件。' if zh else
          'Synthetic example. A later revisit is an outcome to test, not a promise.', 13, MUTED)
    save(fig, f'fvg-{lang}')


def risk(lang, mobile):
    if mobile:
        return risk_mobile(lang)
    zh = lang == 'zh'
    fig = canvas(1200, 740)
    head(fig, '05 / '+('从形态到收益' if zh else 'FROM PATTERN TO PAYOFF'),
         '先算成本，再谈胜率' if zh else 'Price the cost before the edge',
         '独立算术示例 · 每单位往返总成本 0.10' if zh else
         'Separate arithmetic example · round-trip cost 0.10 per unit', False)
    ax = fig.add_axes([.065, .21, .40, .55])
    ax.set_xlim(0, 1); ax.set_ylim(97, 105); ax.axis('off')
    ax.axhspan(100, 104, xmin=.1, xmax=.5, facecolor=UP, alpha=.1)
    ax.axhspan(98, 100, xmin=.1, xmax=.5, facecolor=DOWN, alpha=.1)
    for y, color, name in [(104, UP, '目标' if zh else 'Target'),
                           (100, FG, '入场' if zh else 'Entry'),
                           (98, DOWN, '止损' if zh else 'Stop')]:
        ax.hlines(y, .10, .50, color=color, linewidth=2)
        ax.text(.57, y, f'{name}  {y}', fontsize=21, color=color, va='center')
    ax.text(.13, 102, '+4', color=UP, fontsize=30, weight='bold', va='center')
    ax.text(.13, 99, '−2', color=DOWN, fontsize=30, weight='bold', va='center')
    pairs = [('+1.95R', '获利净值  (4 − 0.10) / 2' if zh else 'Net win  (4 − 0.10) / 2', UP),
             ('−1.05R', '亏损净值  (−2 − 0.10) / 2' if zh else 'Net loss  (−2 − 0.10) / 2', DOWN),
             ('35%', '盈亏平衡胜率  1.05 / 3' if zh else 'Break-even win rate  1.05 / 3', ZONE)]
    for i, (value, caption, color) in enumerate(pairs):
        y = .72 - i * .19
        label(fig, .57, y, value, 32, color, 'bold')
        label(fig, .57, y - .075, caption, 17, MUTED)
    label(fig, .065, .095, 'R = 2；假定成交价分别为 100、98 或 104。' if zh else
          'R = 2; assumed fills are 100, 98 or 104.', 16, FG)
    label(fig, .065, .047, '实际跳空、滑点和部分成交会改变这组结果。' if zh else
          'Gaps, slippage and partial fills can change these outcomes.', 14, MUTED)
    save(fig, f'risk-{lang}')


def cover():
    fig = canvas(1200, 630)
    label(fig, .065, .9, 'FENG / QUANTITATIVE RESEARCH', 13, UP, 'bold')
    label(fig, .065, .73, 'SMART\nMONEY\nCONCEPTS', 43, FG, 'bold')
    label(fig, .065, .23, 'Liquidity. Price formation.\nRisk and expected return.', 18, MUTED)
    ax = fig.add_axes([.50, .17, .455, .62]); axes_style(ax, 12)
    candles(ax, BARS)
    ax.set_xlim(.4, 12.8); ax.set_ylim(95, 113)
    ax.axhspan(102, 105, xmin=.71, xmax=.98, facecolor=ZONE, alpha=.12)
    ax.hlines(98, 5, 12.6, color=ZONE, linestyle=(0, (4, 5)), linewidth=1.2)
    ax.hlines(107, 8, 12.6, color=ZONE, linestyle=(0, (4, 5)), linewidth=1.2)
    ax.set_xticks([]); ax.set_yticks([98, 102, 107, 110]); ax.spines['bottom'].set_visible(False)
    label(fig, .955, .095, 'SYNTHETIC OHLC / NO PERFORMANCE CLAIM', 10, MUTED, ha='right')
    save(fig, 'cover')


def impact_paths():
    """Two assumed paths, not an estimated decomposition of observed prices."""
    time = [i / 10 for i in range(161)]
    pressure = [4 * math.exp(-t / 3.5) for t in time]
    information = [4.0 for _ in time]
    assert time[0] == 0 and time[-1] == 16
    assert pressure[0] == information[0] == 4
    assert all(a > b for a, b in zip(pressure, pressure[1:]))
    assert 0 < pressure[-1] < .05
    assert all(p == 4 for p in information)
    return time, pressure, information


def impact_plot(ax, time, values, color, mobile=False):
    axes_style(ax, 19 if mobile else 16)
    ax.set_xlim(-2, 16.6)
    ax.set_ylim(-.4, 5.2)
    ax.set_xticks([0, 8, 16])
    ax.set_yticks([0, 2, 4])
    ax.axhline(0, color=GRID, linewidth=1)
    ax.axvline(0, color=MUTED, alpha=.35, linestyle=(0, (2, 4)), linewidth=1)
    # The neutral segment and jump are identical in both hypothetical cases.
    ax.plot([-2, 0, 0], [0, 0, 4], color=FG, linewidth=2.2, zorder=4)
    ax.plot(time, values, color=color, linewidth=3.0, zorder=5)
    ax.scatter([0], [4], s=42, color=FG, edgecolor=BG, linewidth=1, zorder=6)
    ax.scatter([16], [values[-1]], s=34, color=color, zorder=6)
    ax.text(.6, 4.7, '+4', fontsize=19 if mobile else 18, color=FG, va='center')


def impact(lang, mobile):
    zh = lang == 'zh'
    time, pressure, information = impact_paths()
    if mobile:
        fig = canvas(480, 1280)
        pxlabel(fig, 32, 30, '02 / ' + ('价格形成' if zh else 'PRICE FORMATION'), 15, UP, 'bold')
        pxlabel(fig, 32, 78, '同一次上涨，\n可能有不同原因' if zh else
                'Same initial move.\nDifferent possible causes.', 22, FG, 'bold')
        pxlabel(fig, 32, 189, '纵轴：价格变化；横轴：时间。\n两者均使用任意单位。' if zh else
                'Vertical: price change.\nHorizontal: time. Arbitrary units.', 18.5, MUTED)
        panels = [(267, 335, pressure, ZONE, 'A  暂时性买入压力' if zh else
                   'A  Temporary buying pressure',
                   '假设：压力消退，价格回落。' if zh else
                   'Assume the pressure fades.'),
                  (676, 744, information, UP, 'B  信息导致重新定价' if zh else
                   'B  Repricing after information',
                   '假设：新的估值持续存在。' if zh else
                   'Assume the new value persists.')]
        for title_y, ax_y, values, color, title, note in panels:
            pxlabel(fig, 32, title_y, title, 19.5, color, 'bold')
            ax = pixel_axes(fig, 76, ax_y, 372, 235)
            impact_plot(ax, time, values, color, True)
            pxlabel(fig, 32, ax_y + 300, note, 19.5, MUTED)
        pxlabel(fig, 32, 1114, '初次上涨不能识别原因。\n真实交易可同时含有两种效应。\n示意路径，不是预测或拟合。' if zh else
                'The initial move cannot reveal cause.\nReal trades may combine both effects.\nAssumed paths, not a forecast.', 18, MUTED)
        save(fig, f'impact-{lang}-mobile')
        return
    fig = canvas(1200, 740)
    label(fig, .065, .956, '02 / ' + ('价格形成' if zh else 'PRICE FORMATION'), 13, UP, 'bold')
    label(fig, .065, .89, '同一次上涨，可能有不同原因' if zh else
          'The same initial move can have different causes', 27, FG, 'bold')
    label(fig, .065, .807, '纵轴为价格变化，横轴为时间，均用任意单位；两条路径由假设给定。' if zh else
          'Vertical: price change. Horizontal: time. Arbitrary units; both paths are assumed.', 16, MUTED)
    for x, values, color, title, note in [
            (.065, pressure, ZONE, 'A  暂时性买入压力' if zh else
             'A  Temporary buying pressure', '假设：压力消退，价格逐步回落。' if zh else
             'Assumption: the pressure fades.'),
            (.555, information, UP, 'B  信息导致重新定价' if zh else
             'B  Repricing after information', '假设：新的估值持续存在。' if zh else
             'Assumption: the new value persists.')]:
        label(fig, x, .695, title, 19, color, 'bold')
        ax = fig.add_axes([x + .035, .29, .37, .33])
        impact_plot(ax, time, values, color)
        label(fig, x + .035, .205, note, 15.5, MUTED)
    label(fig, .065, .104, '初次上涨不能识别原因；真实交易可能同时包含两种效应。' if zh else
          'The initial rise cannot identify the cause. Real trades may combine both effects.', 16, FG)
    label(fig, .065, .046, '示意路径，不是市场预测、经验拟合或收益证据。' if zh else
          'Illustrative paths, not a market forecast, empirical fit or evidence of trading returns.', 13, MUTED)
    save(fig, f'impact-{lang}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', choices=('finance',), help='Render only impact diagrams and the finance cover.')
    args = parser.parse_args()
    checks = validate()
    impact_paths()
    checks['assumed_impact_paths_checked'] = True
    if args.only == 'finance':
        for lang in ('en', 'zh'):
            for mobile in (False, True):
                impact(lang, mobile)
        cover()
        print('Impact-path arithmetic and figure text bounds passed. Rendered four impact diagrams and cover.')
        return
    for lang in ('en', 'zh'):
        for mobile in (False, True):
            for render in (chronology, fvg, risk, impact):
                render(lang, mobile)
    cover()
    checks['figure_text_bounds_checked'] = True
    (OUT / 'checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps(checks, indent=2))
    print('Rendered 16 responsive language figures and 1 cover, each SVG + PNG.')


if __name__ == '__main__':
    main()
