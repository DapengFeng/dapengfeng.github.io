#!/usr/bin/env python3
"""Reproduce Roadster 06's illustrative thermal and gas-recovery calculations.

Run with NumPy, Matplotlib and Noto Sans CJK JP installed:
  python3 scripts/render-roadster-repeatability.py [--output-dir DIRECTORY]

These are explicit design examples, not Tesla measurements or specifications.
The thermal node represents electrical-drive losses only. A 1 MW electrical
port, 3 s pulse and 90% port-to-wheel efficiency give 3 MJ input, 2.7 MJ wheel
work and 300 kJ loss. The loss is deposited into one 60 kJ/K node; constant
conductance 0.5 kW/K cools it to a 25 deg C reservoir. The impulse approximation
omits cooling during the 3 s heat pulse. Period means start-to-start time.
The gas example is separate: ideal air, 300 K cooled endpoints and constant
500 m/s effective exhaust velocity. Pressures are absolute. Real high-pressure
air, fill heating, aftercooling and compressor maps require a richer model.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "roadster-repeatability-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
import numpy as np

Q, C, G, AMBIENT, LIMIT = 300000., 60000., 500., 25., 45.
PULSE_S, INPUT_W, DRIVE_EFF = 3., 1000000., .90
IMPULSE_NS, VE, R, GAS_T = 5000., 500., 287., 300.
P_LOW, P_HIGH, P_ATM, COMP_EFF, COMP_W = 1e7, 3e7, 1e5, .65, 20000.
PAPER, INK, MUTED = "#101b24", "#f2eee3", "#b5c3cc"
GRID, MINT, BLUE, AMBER, RED = "#506778", "#bcecc0", "#8cd1f4", "#f4b780", "#e99d9d"


def peak_rise(n, period, capacity=C, conductance=G, heat=Q):
    if n < 1 or period <= 0 or capacity <= 0 or conductance < 0 or heat < 0:
        raise ValueError("Positive count, period, capacity; nonnegative G and Q required")
    if conductance == 0:
        return n * heat / capacity
    exponent = -period * conductance / capacity
    return heat / capacity * math.expm1(n * exponent) / math.expm1(exponent)


def steady_rise(period, capacity=C, conductance=G, heat=Q):
    if period <= 0 or capacity <= 0 or conductance <= 0 or heat < 0:
        raise ValueError("Positive period, C and G; nonnegative Q required")
    return heat / capacity / -math.expm1(-period * conductance / capacity)


def minimum_period(capacity=C, conductance=G, heat=Q, rise_limit=LIMIT-AMBIENT):
    if capacity <= 0 or conductance <= 0 or heat < 0 or rise_limit <= 0:
        raise ValueError("Invalid thermal parameters")
    if heat == 0:
        return 0.
    if heat / capacity >= rise_limit:
        return math.inf
    return -capacity / conductance * math.log1p(-heat / capacity / rise_limit)


def refill(volume=None, p_low=P_LOW, p_high=P_HIGH, temperature=GAS_T):
    if not P_ATM <= p_low < p_high or temperature <= 0:
        raise ValueError("Refill requires positive temperature and p_atm <= p_low < p_high")
    consumed = IMPULSE_NS / VE
    volume = consumed * R * temperature / (p_high-p_low) if volume is None else volume
    if volume <= 0:
        raise ValueError("Positive volume required")
    primitive = lambda p: p * math.log(p / P_ATM) - p
    reversible = volume * (primitive(p_high)-primitive(p_low))
    electrical = reversible / COMP_EFF
    mass_low, mass_high = p_low*volume/(R*temperature), p_high*volume/(R*temperature)
    adiabatic_temperature = temperature*(mass_low/mass_high)**.4
    return {"consumed_kg": consumed, "volume_m3": volume,
            "cooled_low_mass_kg": p_low*volume/(R*temperature),
            "cooled_high_mass_kg": p_high*volume/(R*temperature),
            "ideal_adiabatic_residual_temperature_K": adiabatic_temperature,
            "ideal_adiabatic_residual_pressure_Pa": mass_low*R*adiabatic_temperature/volume,
            "conditioning_time_s": None,
            "reversible_work_J": reversible, "electrical_input_J": electrical,
            "electrical_input_kWh": electrical/3.6e6,
            "refill_wait_s": electrical/COMP_W,
            "start_to_start_gas_period_s": PULSE_S+electrical/COMP_W}


def check_model():
    # Direct recurrence independently validates the geometric-series solution.
    largest = 0.
    for period in [5.,20.,60.,300.]:
        for capacity in [30000.,60000.,120000.]:
            for conductance in [0.,250.,500.,1000.]:
                state = 0.
                for n in range(1,31):
                    if n > 1:
                        state *= math.exp(-period*conductance/capacity)
                    state += Q/capacity
                    largest = max(largest,abs(state-peak_rise(n,period,capacity,conductance)))
    assert largest < 1e-10
    assert peak_rise(1,20) == Q/C
    assert peak_rise(10,20,conductance=0.) == 50.
    assert peak_rise(10,20,heat=0.) == 0.
    assert abs(steady_rise(10000)-Q/C) < 1e-12
    assert math.isinf(minimum_period(heat=(LIMIT-AMBIENT)*C))
    assert minimum_period(heat=0.) == 0.
    critical = minimum_period()
    assert math.isclose(steady_rise(critical),LIMIT-AMBIENT,abs_tol=1e-12)
    assert steady_rise(critical*.99) > LIMIT-AMBIENT
    assert steady_rise(critical*1.01) < LIMIT-AMBIENT
    assert math.isclose(INPUT_W*PULSE_S,INPUT_W*PULSE_S*DRIVE_EFF+Q,abs_tol=1e-8)
    # RK4 integrates C*dtheta/dt=P_heat-G*theta over finite 3 s heat pulses.
    # This deliberately differs from the instantaneous-heat approximation.
    finite_errors=[]
    finite_results=[]
    dt=.02
    for period in [20.,60.]:
        state=0.
        for n in range(1,11):
            for duration,power in [(PULSE_S,Q/PULSE_S),(period-PULSE_S,0.)]:
                for _ in range(round(duration/dt)):
                    f=lambda value:(power-G*value)/C
                    k1=f(state);k2=f(state+.5*dt*k1);k3=f(state+.5*dt*k2);k4=f(state+dt*k3)
                    state += dt*(k1+2*k2+2*k3+k4)/6
                if power:
                    effective=Q/PULSE_S/G*(-math.expm1(-PULSE_S*G/C))
                    analytic=peak_rise(n,period,heat=effective*C)
                    finite_errors.append(abs(state-analytic))
                    if n == 10:
                        finite_results.append({"period_s":period,"peak10_C":AMBIENT+state})
    assert max(finite_errors) < 1e-9
    gas=refill()
    assert math.isclose(gas["cooled_high_mass_kg"]-gas["cooled_low_mass_kg"],gas["consumed_kg"],abs_tol=1e-12)
    assert math.isclose(gas["consumed_kg"]*VE,IMPULSE_NS,abs_tol=1e-12)
    # Midpoint integration of R*T*ln(p/pa) dm checks reversible refill work.
    dm=gas["consumed_kg"]/100000
    masses=gas["cooled_low_mass_kg"]+(np.arange(100000)+.5)*dm
    pressure=masses*R*GAS_T/gas["volume_m3"]
    numeric=float(np.sum(R*GAS_T*np.log(pressure/P_ATM)*dm))
    relative_error=abs(numeric-gas["reversible_work_J"])/numeric
    assert relative_error < 1e-10
    assert gas["electrical_input_J"] >= gas["reversible_work_J"]
    effective=Q/PULSE_S/G*(-math.expm1(-PULSE_S*G/C))
    assert math.isclose(gas["ideal_adiabatic_residual_pressure_Pa"]/P_HIGH,
                        (gas["cooled_low_mass_kg"]/gas["cooled_high_mass_kg"])**1.4,rel_tol=1e-12)
    return {"recurrence_cases":1440,"maximum_recurrence_residual_K":largest,
            "finite_pulse_RK4_max_residual_K":max(finite_errors),
            "finite_pulse_peak_results":finite_results,
            "finite_pulse_effective_jump_K":effective,
            "finite_pulse_lower_peak_rise_percent":100*(1-effective/(Q/C)),
            "finite_pulse_minimum_period_s":minimum_period(heat=effective*C),
            "refill_work_midpoint_relative_error":relative_error,
            "checks":["energy partition","gas mass and impulse conservation","zero heat",
                      "zero cooling","isolated pulse limit","thermal feasibility boundary"]}


def configure(lang):
    family="DejaVu Sans" if lang=="en" else "Noto Sans CJK JP"
    findfont(FontProperties(family=family),fallback_to_default=False)
    plt.rcParams.update({"font.family":family,"font.size":20,"svg.fonttype":"path",
                        "svg.hashsalt":f"roadster-repeatability-v1-{lang}",
                        "axes.unicode_minus":False,"mathtext.fontset":"dejavusans",
                        "lines.solid_capstyle":"round"})


def canvas(lang,mobile,height=8.2,mheight=11.5):
    configure(lang)
    return plt.figure(figsize=(5.4,mheight) if mobile else (12,height),facecolor=PAPER)


def save(fig,output,stem,lang,mobile,description):
    name=f"{stem}-{'mobile-' if mobile else ''}{lang}"
    path=output/f"{name}.svg"
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for text in fig.texts:
        bounds=text.get_window_extent(renderer)
        assert bounds.x0 >= 0 and bounds.x1 <= fig.bbox.width, (name,text.get_text(),bounds.bounds)
        assert bounds.y0 >= 0 and bounds.y1 <= fig.bbox.height, (name,text.get_text(),bounds.bounds)
    fig.savefig(path,metadata={"Date":None,"Creator":"FENG · Matplotlib","Description":description})
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    ET.parse(path)
    fig.savefig(Path(tempfile.gettempdir())/f"roadster-repeatability-{name}.png",dpi=100)
    plt.close(fig)


def axes(fig,rect,mobile):
    ax=fig.add_axes(rect,facecolor=PAPER)
    ax.grid(axis="y",color=GRID,alpha=.35,lw=.7)
    ax.set_axisbelow(True)
    for side,spine in ax.spines.items():
        spine.set_visible(side in ("left","bottom"));spine.set_color(GRID)
    ax.tick_params(colors=MUTED,labelsize=20 if mobile else 19,length=0,pad=7)
    return ax


def heading(fig,title,sub,mobile):
    fig.text(.075,.955,title,color=INK,fontsize=25 if mobile else 30,va="top",linespacing=1.25)
    fig.text(.075,.865 if mobile else .858,sub,color=MUTED,fontsize=20,va="top",linespacing=1.45)


def render_cycles(output,lang,mobile):
    zh=lang=="zh";fig=canvas(lang,mobile,8.2,11.5)
    heading(fig,"重复功率脉冲的峰值温度" if zh else "Peak temperature under\nrepeated power pulses" if mobile else "Peak temperature under repeated power pulses",
            "每次热输入 300 kJ · 环境 25°C" if zh else "300 kJ per pulse · Ambient 25°C",mobile)
    for i,(color,label,marker) in enumerate([(AMBER,"20 s 周期" if zh else "20 s period","^"),(BLUE,"60 s 周期" if zh else "60 s period","o")]):
        x=.075+i*(.47 if mobile else .29)
        fig.add_artist(plt.Line2D([x,x+.04],[.79,.79],transform=fig.transFigure,color=color,lw=3,marker=marker,markersize=5))
        fig.text(x+.06,.79,label,color=color,fontsize=20,va="center")
    ax=axes(fig,[.15,.31,.79,.39] if mobile else [.09,.28,.85,.42],mobile)
    for period,color,marker in [(20.,AMBER,"^"),(60.,BLUE,"o")]:
        n=np.arange(1,11);ys=[AMBIENT+peak_rise(int(k),period) for k in n]
        # Above-limit demand is a counterfactual extension, not safe operation.
        if period == 20.:
            crossing=5+(LIMIT-ys[4])/(ys[5]-ys[4])
            ax.plot([*n[:5],crossing],[*ys[:5],LIMIT],color=color,lw=2.7)
            ax.plot([crossing,*n[5:]],[LIMIT,*ys[5:]],color=color,lw=2.7,ls=(0,(3,3)))
            ax.scatter(n,ys,color=color,marker=marker,s=30,zorder=4)
        else:
            ax.plot(n,ys,color=color,lw=2.7,marker=marker,markersize=5)
        ax.axhline(AMBIENT+steady_rise(period),color=color,lw=1.2,ls=(0,(4,5)),alpha=.8)
        ax.annotate(f"{ys[-1]:.1f}",(10,ys[-1]),xytext=(-6,10 if period==20 else -25),textcoords="offset points",ha="right",color=color,fontsize=22)
    ax.axhline(LIMIT,color=RED,lw=1.5,ls=(0,(2,4)))
    ax.text(1.2,46.2,"45°C 上限" if zh else "45°C limit",color=RED,fontsize=20)
    ax.scatter([6],[AMBIENT+peak_rise(6,20)],s=110,facecolor=PAPER,edgecolor=AMBER,lw=2,zorder=6)
    ax.set(xlim=(.65,10.45),ylim=(25,62),xticks=[1,3,5,7,10],yticks=[25,35,45,55])
    ax.text(0,1.10,"峰温 / °C" if zh else "Peak temperature / °C",transform=ax.transAxes,color=INK,fontsize=21)
    ax.set_xlabel("功率脉冲次数" if zh else "Pulse number",color=MUTED,fontsize=21,labelpad=10)
    footer=("20 s：第 6 次超过上限\n水平虚线：周期稳态 57.6 / 37.7°C\n超限曲线：未限功率的计算延拓" if zh else "20 s: first breach on pulse 6.\nSteady peaks: 57.6 / 37.7°C.\nAbove-limit curve assumes\nno power reduction.") if mobile else ("20 s 周期第 6 次超限；水平虚线对应周期稳态 57.6°C / 37.7°C。\n超过 45°C 的虚线曲线是未限功率的计算延拓。" if zh else "20 s exceeds 45°C on pulse 6; periodic peaks are 57.6°C / 37.7°C.\nThe dashed above-limit curve assumes no power reduction.")
    fig.text(.075,.045 if mobile else .06,footer,color=MUTED,fontsize=19 if mobile else 20,va="bottom",linespacing=1.5)
    save(fig,output,"thermal-cycles",lang,mobile,"Illustrative lumped drive-loss node: heat 300 kJ, capacity 60 kJ/K, conductance 0.5 kW/K, ambient 25 C. Instantaneous-pulse approximation, periods measured start to start. 20 s period exceeds 45 C at run 6 and reaches 51.418 C at run 10; periodic peak 57.569 C. 60 s reaches 37.622 C at run 10, periodic peak 37.707 C. No actual Tesla data.")


def render_design(output,lang,mobile):
    zh=lang=="zh";fig=canvas(lang,mobile,8.5,12.7)
    heading(fig,"热容与散热能力的不同作用" if zh else "Thermal capacity and\ncooling conductance" if mobile else "Thermal capacity and cooling conductance",
            "同为 20 s 周期 · 每次 300 kJ" if zh else "20 s period · 300 kJ per pulse",mobile)
    labels=["基准","热容 ×2","散热 ×2"] if zh else ["Base","C ×2","G ×2"]
    colors=[AMBER,BLUE,MINT];parameters=[(C,G),(2*C,G),(C,2*G)]
    for idx,steady in enumerate([False,True]):
        rect=[.16,.49-idx*.30,.77,.20] if mobile else [.09+idx*.47,.29,.38,.39]
        ax=axes(fig,rect,mobile)
        vals=[AMBIENT+(steady_rise(20,c,g) if steady else peak_rise(10,20,c,g)) for c,g in parameters]
        ax.bar(np.arange(3),np.array(vals)-AMBIENT,bottom=AMBIENT,width=.56,color=colors,zorder=3)
        for j,value in enumerate(vals):
            ax.text(j,value+1,f"{value:.1f}",ha="center",color=colors[j],fontsize=22)
        ax.axhline(LIMIT,color=RED,lw=1.2,ls=(0,(2,4)))
        ax.set(ylim=(25,64),xlim=(-.65,2.65),xticks=[0,1,2],xticklabels=labels,yticks=[25,35,45,55])
        ax.text(0,1.13,("周期稳态 / °C" if steady else "第 10 次 / °C") if zh else ("Periodic limit / °C" if steady else "Pulse 10 / °C"),transform=ax.transAxes,color=INK,fontsize=22)
    text=("热容翻倍可延缓升温，长期仍超限。\n散热翻倍使周期稳态降到 42.6°C。\n基准：C = 60 kJ/K，G = 0.5 kW/K" if zh else "Doubling C delays accumulation;\nits periodic peak still exceeds 45°C.\nDoubling G gives a 42.6°C limit.") if mobile else ("热容翻倍延缓升温，但周期稳态仍超过 45°C；散热翻倍使周期稳态降到 42.6°C。\n基准 C = 60 kJ/K，G = 0.5 kW/K；环境 25°C。" if zh else "Doubling C delays accumulation but still exceeds 45°C in periodic operation.\nDoubling G gives 42.6°C. Base: C = 60 kJ/K, G = 0.5 kW/K; ambient 25°C.")
    fig.text(.075,.035 if mobile else .065,text,color=MUTED,fontsize=18 if mobile else 19,va="bottom",linespacing=1.55)
    save(fig,output,"thermal-design",lang,mobile,"Illustrative 20 s repetition. Baseline C=60 kJ/K and G=0.5 kW/K has run10 / periodic peaks 51.418 / 57.569 C. Doubling capacity alone gives 42.679 / 56.267 C; doubling conductance alone gives 42.009 / 42.639 C. Every run deposits 300 kJ; ambient is 25 C. Increasing C mostly buys transient thermal buffer; it does not replace sufficient heat rejection.")


def render_refill(output,lang,mobile):
    zh=lang=="zh";fig=canvas(lang,mobile,8.4,12.0);gas=refill()
    heading(fig,"补气需要恢复质量与温度" if zh else "Refill restores mass\nat a defined temperature" if mobile else "Refill restores mass at a defined temperature",
            "定温状态 300 K · 理想气体 · 43.05 L" if zh else "300 K conditioned states · 43.05 L",mobile)
    ax=axes(fig,[.17,.44,.74,.26] if mobile else [.09,.28,.48,.40],mobile)
    masses=np.linspace(5,15,101);pressure=masses*R*GAS_T/gas["volume_m3"]/1e6
    ax.plot(masses,pressure,color=BLUE,lw=3)
    ax.scatter([5,15],[10,30],s=65,color=BLUE,zorder=5)
    ax.set(xlim=(3.5,16.5),ylim=(5,36),xticks=[5,10,15],yticks=[10,20,30])
    ax.set_xlabel("罐内气体质量 / kg" if zh else "Gas mass in tank / kg",color=MUTED,fontsize=20,labelpad=12)
    ax.text(0,1.12,"绝对压力 / MPa" if zh else "Absolute pressure / MPa",transform=ax.transAxes,color=INK,fontsize=20)
    ax.annotate("5 kg · 10 MPa",(5,10),xytext=(2,13),textcoords="offset points",color=BLUE,fontsize=19)
    ax.annotate("15 kg · 30 MPa",(15,30),xytext=(-4,12),textcoords="offset points",ha="right",color=BLUE,fontsize=19)
    x=.075 if mobile else .66
    ys=[.31,.23,.15] if mobile else [.62,.46,.30]
    entries=[("10 kg","5 kN·s / 500 m/s"),("1.93 kWh","补气 · 效率 65%" if zh else "Refill · 65% efficiency"),("347.9 s","20 kW 输入" if zh else "20 kW input")]
    for y,(value,label) in zip(ys,entries):
        fig.text(x,y,value,color=INK if value=="10 kg" else AMBER,fontsize=26 if mobile else 29,va="center")
        fig.text(x+(.36 if mobile else 0),y if mobile else y-.055,label,color=MUTED,fontsize=18 if mobile else 17,va="center")
    text=("起点：残余气体已调温至 300 K。\n调温所需时间未计入补气时间；\n曲线不是实际快速充气轨迹。" if zh else "Start: gas already at 300 K.\nConditioning time is additional.\nNot the actual fast-fill trajectory.")
    fig.text(.075,.028 if mobile else .055,text,color=MUTED,fontsize=19,va="bottom",linespacing=1.5)
    save(fig,output,"gas-refill",lang,mobile,"Separate ideal-gas example at cooled endpoints 300 K. 5 kNs impulse at constant effective exhaust velocity 500 m/s consumes 10 kg air. Tank 43.05 L holds 5 kg at 10 MPa and 15 kg at 30 MPa absolute. Ambient inlet 0.1 MPa; reversible isothermal refill work 4.523 MJ. Electrical input at assumed 65% overall efficiency is 6.958 MJ or 1.933 kWh, requiring 347.916 s at 20 kW. This is not a real fast-fill trajectory or a high-pressure vessel design; pressure alone cannot establish recovered gas mass or readiness.")


def render_recovery(output,lang,mobile):
    zh=lang=="zh";fig=canvas(lang,mobile,8.4,11.9);gas=refill()
    heading(fig,"重复周期由恢复瓶颈决定" if zh else "Recovery limits\nthe repetition period" if mobile else "Recovery limits the repetition period",
            "乐观参考：未计气体调温时间" if zh else "Optimistic reference:\nconditioning omitted" if mobile else "Optimistic reference: conditioning omitted",mobile)
    ax=axes(fig,[.17,.35,.75,.36] if mobile else [.09,.30,.59,.40],mobile)
    power=np.linspace(10,100,300);gasperiod=PULSE_S+gas["electrical_input_J"]/power/1000
    ax.plot(power,gasperiod,color=AMBER,lw=3)
    ax.axhline(minimum_period(),color=BLUE,lw=2,ls=(0,(5,4)))
    ax.scatter([20],[gas["start_to_start_gas_period_s"]],s=70,color=AMBER,zorder=5)
    ax.annotate("20 kW · 350.9 s",(20,gas["start_to_start_gas_period_s"]),xytext=(8,12),textcoords="offset points",color=AMBER,fontsize=20)
    ax.set(xlim=(8,103),ylim=(0,760),xticks=[20,40,60,80,100] if not mobile else [20,60,100],yticks=[0,200,400,600])
    ax.set_xlabel("压缩机输入功率 / kW" if zh else "Compressor electrical power / kW",color=MUTED,fontsize=19 if mobile else 20,labelpad=12)
    ax.text(0,1.12,"起动到起动周期 / s" if zh else "Start-to-start period / s",transform=ax.transAxes,color=INK,fontsize=20)
    if mobile:
        fig.text(.075,.235,"补气：3 + 347.9 = 350.9 s" if zh else "Gas: 3 + 347.9 = 350.9 s",color=AMBER,fontsize=20)
        fig.text(.075,.19,"热约束：34.5 s" if zh else "Thermal constraint: 34.5 s",color=BLUE,fontsize=20)
    else:
        fig.text(.74,.61,"补气约束" if zh else "Gas constraint",color=AMBER,fontsize=21)
        fig.text(.74,.55,"350.9 s",color=AMBER,fontsize=29)
        fig.text(.74,.44,"热约束" if zh else "Thermal constraint",color=BLUE,fontsize=21)
        fig.text(.74,.38,"34.5 s",color=BLUE,fontsize=29)
    text=("此处设气体调温时间为零。\n实际补气路径：3 s + 调温 + 347.9 s。\n独立并行时再与热约束取较大值。" if zh else "Here, conditioning time is set to zero.\nGas path: 3 s + conditioning + 347.9 s.\nIndependent parallel recovery assumed.") if mobile else ("此处设气体调温时间为零；实际补气路径为 3 s + 调温时间 + 347.9 s。\n仅在恢复资源独立且可并行时，与热约束取较大值；此图不证明整车已就绪。" if zh else "Gas recovery takes 3 s + conditioning + 347.9 s; conditioning is set to zero here.\nThe greater thermal/gas bound assumes independent parallel resources.\nThis reference does not establish vehicle readiness.")
    fig.text(.075,.035 if mobile else .06,text,color=MUTED,fontsize=17 if mobile else 19,va="bottom",linespacing=1.5)
    save(fig,output,"recovery-limits",lang,mobile,"Illustrative optimistic start-to-start reference with gas conditioning time set to zero. Impulsive thermal limit is 34.522 s. Defined gas recovery path: 3 s pulse, then condition residual gas to 300 K, then refill isothermally from 10 to 30 MPa. At 20 kW the refill takes 347.916 s; the conditional period is 350.916 s plus unknown conditioning time. The greater thermal/gas bound applies only with independent parallel resources; no vehicle readiness is proven.")


def render_cover(output):
    configure("en");fig=plt.figure(figsize=(12,6.3),facecolor=PAPER)
    fig.text(.075,.87,"ROADSTER 06",color=MUTED,fontsize=18,weight="bold")
    fig.text(.075,.74,"REPEATABLE PERFORMANCE",color=INK,fontsize=31,weight="bold")
    fig.text(.075,.56,"Thermal capacity. Cooling. Gas recovery.",color=MUTED,fontsize=20)
    ax=fig.add_axes([.54,.11,.40,.37],facecolor=PAPER);ax.axis("off")
    for period,color in [(20.,AMBER),(60.,BLUE)]:
        n=np.arange(1,11);ys=[AMBIENT+peak_rise(int(i),period) for i in n]
        ax.plot(n,ys,color=color,lw=3,marker="o",markersize=4)
    ax.axhline(LIMIT,color=RED,lw=1.2,ls=(0,(3,5)))
    ax.set(xlim=(.7,10.3),ylim=(25,55))
    fig.text(.075,.33,"34.5 s",color=BLUE,fontsize=35)
    fig.text(.075,.25,"thermal cycle bound",color=MUTED,fontsize=18)
    fig.text(.075,.13,"350.9 s",color=AMBER,fontsize=35)
    fig.text(.075,.07,"idealized gas-recovery reference",color=MUTED,fontsize=16)
    fig.text(.075,.025,"Illustrative model · conditioning omitted",color=MUTED,fontsize=13)
    fig.savefig(output/"cover.png",dpi=100,metadata={"Description":"Roadster 06 illustrative repeatability: thermal and gas recovery constraints. Model assumptions, not Tesla performance."})
    plt.close(fig)


def summary(checks):
    return {"boundary":"Illustrative design calculations, not Tesla specifications, measurements or a validated vehicle simulation.",
        "thermal":{"heat_J":Q,"heat_capacity_J_K":C,"conductance_W_K":G,"ambient_C":AMBIENT,
           "limit_C":LIMIT,"time_constant_s":C/G,"instantaneous_rise_K":Q/C,
           "drive_electrical_input_W":INPUT_W,"pulse_duration_s":PULSE_S,"port_to_wheel_efficiency":DRIVE_EFF,
           "drive_input_J":INPUT_W*PULSE_S,"wheel_work_J":INPUT_W*PULSE_S*DRIVE_EFF,
           "period_definition":"start to start; impulse approximation deposits heat at each 3 s pulse end; consecutive deposits remain one period apart",
           "minimum_period_s":minimum_period(),
           "cases":[{"period_s":p,"peak_temperatures_C":[AMBIENT+peak_rise(n,p) for n in range(1,11)],"periodic_peak_C":AMBIENT+steady_rise(p)} for p in [20.,60.]],
           "design_comparisons":[{"name":name,"capacity_J_K":c,"conductance_W_K":g,"period_s":20.,"peak10_C":AMBIENT+peak_rise(10,20,c,g),"periodic_peak_C":AMBIENT+steady_rise(20,c,g),"minimum_period_s":minimum_period(c,g)} for name,c,g in [("base",C,G),("double_C",2*C,G),("double_G",C,2*G)]]},
        "gas":{"impulse_N_s":IMPULSE_NS,"effective_exhaust_velocity_m_s":VE,"specific_gas_constant_J_kg_K":R,
           "cooled_endpoint_temperature_K":GAS_T,"low_absolute_pressure_Pa":P_LOW,"high_absolute_pressure_Pa":P_HIGH,
           "inlet_absolute_pressure_Pa":P_ATM,"overall_isothermal_efficiency":COMP_EFF,"compressor_electrical_power_W":COMP_W,**refill()},
        "equations":{"thermal_ODE":"C*dtheta/dt=P_heat-G*theta","pulse_recurrence":"theta_n=exp(-Delta/tau)*theta_(n-1)+Q/C, theta_0=0",
            "periodic_peak":"theta_inf=(Q/C)/(1-exp(-Delta/tau))","minimum_period":"-tau*ln(1-Q/(C*theta_limit)); Q/C < theta_limit",
            "gas_mass":"m=p*V/(R*T)","reversible_refill":"V*[p*ln(p/p_atm)-p] evaluated from p_low to p_high",
            "gas_period":"pulse_duration+conditioning_time+reversible_refill/(efficiency*compressor_power); conditioning_time unknown; plotted reference sets it to zero"},
        "assumptions":["Drive-loss heat node only; battery I^2R and compressor heat are outside this node",
            "Constant C,G,ambient; no internal gradients, active chiller, coolant transients, speed dependence or thermal derating",
            "Instantaneous heat approximation ignores cooling during the 3 s pulse; finite-pulse RK4 comparison is reported",
            "Thermal and gas examples use independent temperature assumptions: 25 deg C reservoir and 300 K cooled gas",
            "Ideal gas and fixed effective exhaust velocity; real air at 10-30 MPa needs real-gas properties",
            "Reversible isothermal compression into a tank maintained at 300 K; 65% efficiency refers to this minimum work, including electrical losses",
            "Defined reference path: after 3 s pulse, residual gas must first be conditioned to 300 K before isothermal refill. Conditioning time is unknown and set to zero only in the optimistic plotted reference",
            "Adiabatic ideal-gas discharge of 15 to 5 kg would leave about 193 K and 6.44 MPa, not 300 K and 10 MPa. This separate limiting model is not a real discharge prediction",
            "Maximum of recovery bounds valid only for concurrent independent resources; no vehicle-wide sufficiency claim",
            "Drive input already includes Q; compressor electrical input is not all battery heat and must not be double counted"],"checks":checks}


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parents[1]/"content/assets/roadster-repeatability")
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    checks=check_model()
    for lang in ["en","zh"]:
        for mobile in [False,True]:
            for render in [render_cycles,render_design,render_refill,render_recovery]:render(args.output_dir,lang,mobile)
    render_cover(args.output_dir)
    (args.output_dir/"model.json").write_text(json.dumps(summary(checks),ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"output":str(args.output_dir),"checks":checks,"minimum_thermal_period_s":minimum_period(),"gas":refill()},ensure_ascii=False,indent=2))


if __name__=="__main__":main()
