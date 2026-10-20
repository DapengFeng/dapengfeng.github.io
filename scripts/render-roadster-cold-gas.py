#!/usr/bin/env python3
"""Reproduce Roadster 04's illustrative cold-gas blowdown calculation and figures.

Requires NumPy, Matplotlib and Noto Sans CJK JP. Run:
    python3 scripts/render-roadster-cold-gas.py [--output-dir DIRECTORY]

This is NOT a Tesla specification, hardware design, measured performance or a
real-air prediction. Dry air is represented as a calorically perfect ideal gas
even at 30 MPa and very low final temperatures. Those conditions make real-gas
properties, heat exchange and component testing consequential. The isothermal
case is a different idealized heat boundary, not a certified real-system bound.

Both arrangements have a 100 L rigid tank, initially at 30 MPa absolute and
300 K. Each has a different, fixed-area, ideal converging nozzle sized to give
10 kN initially. One nozzle connects directly to the tank. The other is fed by
an ideal regulator holding 3 MPa, with zero required pressure headroom or line
loss. Plenum stagnation temperature equals tank temperature (ideal-gas,
adiabatic throttling); nozzle expansion is isentropic. Exit equals throat.
Pressure thrust is retained. Both histories stop when tank pressure reaches
3 MPa. A separate 8 kN criterion defines the direct system's service interval;
it is an arbitrary mission threshold, not a standard or a 10 kN rating.

The adiabatic well-mixed tank obeys d(m cv T)/dt = -mdot cp T. Closed-form
histories are checked against an independently integrated mass/internal-energy
ODE, plus numerical mass, energy and impulse integrals. No post-regulation tail
or subcritical nozzle regime is included. SVGs have EN/ZH and narrow versions;
inspection PNGs are written to /tmp, not shipped.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "roadster-cold-gas-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.patches import FancyBboxPatch, Polygon
import numpy as np

GAMMA, R, TI = 1.4, 287., 300.
PI, PR, PA, VOLUME, TARGET = 30e6, 3e6, 101325., .1, 10000.
CV, CP = R / (GAMMA - 1), GAMMA * R / (GAMMA - 1)
CRITICAL = (2 / (GAMMA + 1)) ** (GAMMA / (GAMMA - 1))
THRUST_FACTOR = (GAMMA + 1) * CRITICAL
FLOW_FACTOR = math.sqrt(GAMMA / R) * (2 / (GAMMA + 1)) ** ((GAMMA + 1) / (2 * (GAMMA - 1)))
M0 = PI * VOLUME / (R * TI)
B = (GAMMA - 1) / 2
PAPER, INK, MUTED = "#101b24", "#f2eee3", "#b5c3cc"
GRID, MINT, BLUE, AMBER = "#506778", "#bcecc0", "#8cd1f4", "#f4b780"


def nozzle_area(p0):
    if not math.isfinite(p0) or p0 <= PA / CRITICAL:
        raise ValueError("Choked upstream absolute pressure required")
    return TARGET / (THRUST_FACTOR * p0 - PA)


def nozzle(p0, t0, area):
    p0, t0 = np.asarray(p0), np.asarray(t0)
    if np.any(p0 <= PA / CRITICAL) or np.any(t0 <= 0) or area <= 0:
        raise ValueError("Positive temperature/area and choked pressure ratio required")
    mdot = area * p0 / np.sqrt(t0) * FLOW_FACTOR
    ve = np.sqrt(GAMMA * R * t0 * 2 / (GAMMA + 1))
    pressure_force = (p0 * CRITICAL - PA) * area
    return mdot, mdot * ve + pressure_force, ve, pressure_force


def x_end(thermal="adiabatic"):
    return (PR / PI) ** (1 / GAMMA) if thermal == "adiabatic" else PR / PI


def duration(x, regulated=True, thermal="adiabatic"):
    if thermal not in ("adiabatic", "isothermal"):
        raise ValueError("Unknown thermal boundary")
    x = np.asarray(x)
    if np.any(x <= 0) or np.any(x > 1):
        raise ValueError("Remaining mass fraction must lie in (0,1]")
    p0 = PR if regulated else PI
    flow0 = float(nozzle(p0, TI, nozzle_area(p0))[0])
    if thermal == "isothermal":
        return M0 / flow0 * ((1 - x) if regulated else -np.log(x))
    return M0 / flow0 * ((1 - x ** (1 + B)) / (1 + B) if regulated else (x ** (-B) - 1) / B)


def history(regulated=True, thermal="adiabatic", samples=1201, *, until_fraction=None):
    exponent = GAMMA if thermal == "adiabatic" else 1.
    x = np.linspace(1., x_end(thermal) if until_fraction is None else until_fraction, samples)
    temp = TI * x ** (exponent - 1)
    pressure = PI * x ** exponent
    plenum = np.full_like(x, PR) if regulated else pressure
    area = nozzle_area(PR if regulated else PI)
    flow, force, ve, pressure_force = nozzle(plenum, temp, area)
    return {"time": duration(x, regulated, thermal), "mass": M0 * x, "temperature": temp,
            "pressure": pressure, "flow": flow, "force": force, "ve": ve,
            "pressure_force": pressure_force, "x": x}


def force_fraction_to_mass(fraction):
    if not 0 < fraction <= 1:
        raise ValueError("Force fraction must lie in (0,1]")
    return ((fraction * TARGET / nozzle_area(PI) + PA) / (THRUST_FACTOR * PI)) ** (1 / GAMMA)


def impulse(x, regulated=True, thermal="adiabatic"):
    if regulated:
        return float(TARGET * duration(x, True, thermal))
    flow0 = float(nozzle(PI, TI, nozzle_area(PI))[0])
    if thermal == "isothermal":
        integral = THRUST_FACTOR * PI * (1 - x) + PA * math.log(x)
    else:
        q = (GAMMA + 1) / 2
        integral = THRUST_FACTOR * PI * (1 - x ** q) / q - PA * (1 - x ** (1 - q)) / (1 - q)
    return M0 / flow0 * nozzle_area(PI) * integral


def sample_at_time(t, regulated=True):
    p0 = PR if regulated else PI
    flow0 = float(nozzle(p0, TI, nozzle_area(p0))[0])
    end = float(duration(x_end(), regulated))
    if not 0 <= t <= end:
        raise ValueError("Time is outside the modeled interval")
    if regulated:
        x = (1 - (1 + B) * flow0 / M0 * t) ** (1 / (1 + B))
    else:
        x = (1 + B * flow0 / M0 * t) ** (-1 / B)
    temp, p = TI * x ** (GAMMA - 1), PI * x ** GAMMA
    f, thrust, ve, fp = nozzle(PR if regulated else p, temp, nozzle_area(p0))
    return {"time_s": t, "tank_pressure_Pa_abs": p, "tank_temperature_K": temp,
            "gas_remaining_kg": M0*x, "mass_flow_kg_s": float(f), "thrust_N": float(thrust),
            "jet_velocity_m_s": float(ve), "pressure_thrust_N": float(fp)}


def check_model():
    """Independent ODE and conservation integrals guard the closed-form curves."""
    errors = {}
    for regulated in [True, False]:
        name = "regulated" if regulated else "direct"
        area = nozzle_area(PR if regulated else PI)
        end = float(duration(x_end(), regulated))
        state = np.array([M0, M0 * CV * TI])
        def deriv(y):
            m, u = y
            temp = u / (m * CV)
            p = m * R * temp / VOLUME
            flow = float(nozzle(PR if regulated else p, temp, area)[0])
            return np.array([-flow, -flow * CP * temp])
        # Integrate two extensive states, never the closed-form x,T,p relations.
        steps = 6000
        dt = end / steps
        for _ in range(steps):
            k1 = deriv(state)
            k2 = deriv(state + dt*k1/2)
            k3 = deriv(state + dt*k2/2)
            k4 = deriv(state + dt*k3)
            state += dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        expected = np.array([M0*x_end(), M0*x_end()*CV*TI*x_end()**(GAMMA-1)])
        rel = float(np.max(np.abs(state-expected)/expected))
        assert rel < 1e-8, (name, state, expected)
        h = history(regulated, samples=20001)
        trap = np.trapezoid
        mass_error = abs(trap(h["flow"], h["time"]) - (M0-h["mass"][-1]))/M0
        energy_out = trap(h["flow"]*CP*h["temperature"],h["time"])
        energy_error = abs(energy_out-(M0*CV*TI-h["mass"][-1]*CV*h["temperature"][-1]))/(M0*CV*TI)
        j = impulse(x_end(),regulated)
        impulse_error = abs(trap(h["force"], h["time"])-j)/j
        assert max(mass_error, energy_error, impulse_error) < 1e-7
        # Independent momentum expression agrees with temperature-cancelled thrust.
        assert np.allclose(h["force"], area*(THRUST_FACTOR*(PR if regulated else h["pressure"])-PA),rtol=1e-12)
        assert np.all(h["pressure"] >= PR-1e-6)
        assert np.all((PR if regulated else h["pressure"]) * CRITICAL > PA)
        errors[name] = {"ode_max_relative_error": rel, "mass_integral_relative_error":mass_error,
                        "energy_integral_relative_error":energy_error,"impulse_integral_relative_error":impulse_error}
    assert math.isclose(float(sample_at_time(0)["thrust_N"]), TARGET, rel_tol=1e-12)
    assert math.isclose(float(sample_at_time(float(duration(force_fraction_to_mass(.8),False)),False)["thrust_N"]), .8*TARGET,rel_tol=1e-12)
    # A constant-temperature tank needs Qdot = mdot R T, not zero heat flow.
    iso = history(True,"isothermal",samples=20001)
    heat = np.trapezoid(iso["flow"]*R*TI,iso["time"])
    assert math.isclose(heat, (M0-iso["mass"][-1])*R*TI,rel_tol=1e-12)
    warm_pressure=M0*x_end()*R*TI/VOLUME
    pressure_grid=np.linspace(warm_pressure,PI,20001)
    refill_numeric=VOLUME*np.trapezoid(np.log(pressure_grid/PA),pressure_grid)
    primitive=lambda p:VOLUME*(p*math.log(p/PA)-p)
    refill_exact=primitive(PI)-primitive(warm_pressure)
    refill_error=abs(refill_numeric-refill_exact)/refill_exact
    assert refill_error < 1e-9
    errors["isothermal_refill_integral_relative_error"]=refill_error
    return errors


def configure(lang):
    family = "DejaVu Sans" if lang == "en" else "Noto Sans CJK JP"
    findfont(FontProperties(family=family), fallback_to_default=False)
    plt.rcParams.update({"font.family":family,"font.size":20,"svg.fonttype":"path",
        "svg.hashsalt":f"roadster-cold-gas-v1-{lang}","axes.unicode_minus":False,
        "mathtext.fontset":"dejavusans","lines.solid_capstyle":"round"})


def axis_style(ax, mobile=False):
    ax.set_facecolor(PAPER)
    ax.grid(axis="y",color=GRID,alpha=.3,lw=.7)
    ax.set_axisbelow(True)
    for k,s in ax.spines.items():
        s.set_visible(k in ("left","bottom"));s.set_color(GRID)
    ax.tick_params(colors=MUTED,labelsize=22 if mobile else 20,length=0,pad=9)


def save(fig,output,stem,lang,mobile,description):
    name=f"{stem}-{'mobile-' if mobile else ''}{lang}"
    path=output/f"{name}.svg"
    fig.savefig(path,metadata={"Date":None,"Creator":"FENG · Matplotlib","Description":description})
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    ET.parse(path)
    fig.savefig(Path(tempfile.gettempdir())/f"roadster-cold-gas-{name}.png",dpi=100)
    plt.close(fig)


def title(fig,zh,mobile,en,cn,subtitle_en,subtitle_cn):
    fig.text(.07 if mobile else .055,.965,cn if zh else en,color=INK,fontsize=26 if mobile else 29,va="top",linespacing=1.2)
    fig.text(.07 if mobile else .055,.875 if mobile else .865,subtitle_cn if zh else subtitle_en,color=MUTED,fontsize=21 if mobile else 20,va="top",linespacing=1.4)


def render_thrust(output,lang,mobile):
    configure(lang);zh=lang=="zh"
    fig=plt.figure(figsize=(6,10.7) if mobile else (12,8.4),facecolor=PAPER)
    title(fig,zh,mobile,"Same initial thrust.\nDifferent usable burst." if mobile else "Same initial thrust. Different usable burst.","峰值相同，持续时间不同",
          "100 L · 30 MPa · 300 K\nIdeal adiabatic tank" if mobile else "100 L · 30 MPa · 300 K · Ideal adiabatic tank",
          "100 L · 30 MPa · 300 K\n理想绝热气瓶" if mobile else "100 L · 30 MPa · 300 K · 理想绝热气瓶")
    ax=fig.add_axes([.16,.345,.79,.30] if mobile else [.09,.34,.855,.39]);axis_style(ax,mobile)
    for reg,c in [(False,AMBER),(True,BLUE)]:
        h=history(reg)
        ax.plot(h["time"],h["force"]/1000,color=c,lw=3)
        xcut=x_end() if reg else force_fraction_to_mass(.8)
        active=history(reg,until_fraction=xcut)
        ax.fill_between(active["time"],active["force"]/1000,0,color=c,alpha=.10 if reg else .16)
        ax.scatter(h["time"][-1],h["force"][-1]/1000,s=45,color=c,edgecolors=PAPER,zorder=5)
    end=float(duration(x_end()))
    direct80=float(duration(force_fraction_to_mass(.8),False))
    ax.axhline(8,color=MUTED,lw=1,ls=(0,(4,5)),alpha=.8)
    ax.scatter([direct80],[8],color=AMBER,s=40,edgecolor=PAPER,lw=.7,zorder=6)
    ax.text(3.98,8.15,"8 kN",ha="right",va="bottom",color=MUTED,fontsize=22 if mobile else 20)
    ax.text(1 if not mobile else 0,1.08,"推力 / kN" if zh else "Thrust / kN",transform=ax.transAxes,
            color=INK,fontsize=22,ha="right" if not mobile else "left")
    ax.set(xlim=(0,4),ylim=(0,11.8),xticks=[0,1,2,3,4],yticks=[0,4,8,12])
    ax.set_xlabel("时间 / s" if zh else "Time / s",color=MUTED,fontsize=22,labelpad=9)
    # Two short legends directly above the plot, no crossing annotation leaders.
    ly=.755 if mobile else .78
    for j,(c,label) in enumerate([(BLUE,"稳压至 3 MPa" if zh else "Regulated at 3 MPa"),(AMBER,"直接放气" if zh else "Direct blowdown")]):
        lx=.16 if mobile else .09+j*.41
        yy=ly-j*.047 if mobile else ly
        fig.add_artist(plt.Line2D([lx,lx+.07],[yy,yy],transform=fig.transFigure,color=c,lw=3))
        fig.text(lx+.09,yy,label,color=c,fontsize=21 if mobile else 20,va="center")
    # Service windows are explicit durations, separate from the low-force tail.
    by=.21 if mobile else .205
    if not mobile:
        fig.text(.09,by+.05,"任务区间" if zh else "Service interval",fontsize=23,color=INK)
    for j,(c,label) in enumerate([(BLUE,"10 kN：1.32 s" if zh else "10 kN: 1.32 s"),(AMBER,"≥ 8 kN：0.305 s" if zh else "≥ 8 kN: 0.305 s")]):
        lx=.07 if mobile else .09+j*.43
        yy=by-j*.055 if mobile else by-.017
        fig.text(lx,yy,label,color=c,fontsize=24 if mobile else 25,va="center")
    footer=(("终点：瓶压 3 MPa，仍有余气\n蓝线结束仅表示稳压终止\n阴影：各自任务区间内的冲量" if mobile else
             "终点：瓶压 3 MPa，仍有余气；蓝线结束仅表示稳压终止。\n阴影：各自任务区间内的冲量。") if zh else
            ("Endpoints: 3 MPa; gas remains.\nBlue line ends when regulation stops.\nShading: service-interval impulse." if mobile else
             "Endpoints: 3 MPa; gas remains. Blue line ends when regulation stops.\nShading: impulse during each service interval."))
    fig.text(.07 if mobile else .09,.025,footer,color=MUTED,fontsize=18 if mobile else 19,va="bottom",linespacing=1.4)
    save(fig,output,"thrust-history",lang,mobile,
         "Ideal-air illustrative model, not Tesla data. Different fixed nozzle areas give the same initial 10 kN. Regulated 3 MPa nozzle stays at 10 kN for 1.322 s; direct tank nozzle drops below 8 kN after 0.305 s. Both curves stop at 3 MPa tank pressure; the direct tail is not rated 10 kN service. Shaded areas integrate their different service intervals.")


def render_state(output,lang,mobile):
    configure(lang);zh=lang=="zh"
    fig=plt.figure(figsize=(6,13.5) if mobile else (12,7.4),facecolor=PAPER)
    title(fig,zh,mobile,"Constant thrust.\nRising mass flow." if mobile else "Constant thrust. Rising mass flow.","推力不变，耗气加快",
          "Regulated 3 MPa · 10 kN\nSame 1.32 s burst" if mobile else "Regulated 3 MPa · 10 kN · The same 1.32 s burst",
          "3 MPa 稳压 · 10 kN\n同一次 1.32 s 喷射" if mobile else "3 MPa 稳压 · 10 kN · 同一次 1.32 s 喷射")
    h=history();end=h["time"][-1]
    entries=[("pressure",1e6,AMBER,"气瓶压力 / MPa" if zh else "Tank pressure / MPa",[0,15,30],(0,33),"30 → 3"),
             ("temperature",1,BLUE,"气体温度 / K" if zh else "Gas temperature / K",[150,225,300],(135,320),"300 → 155"),
             ("flow",1,MINT,"质量流率 / kg/s" if zh else "Mass flow / kg/s",[15,20,25],(15,29),"18.9 → 26.3")]
    for i,(key,scale,c,label,ticks,limits,change) in enumerate(entries):
        rect=[.19,.60-i*.215,.75,.15] if mobile else [.075+i*.32,.28,.22,.405]
        ax=fig.add_axes(rect);axis_style(ax,mobile)
        ax.plot(h["time"],h[key]/scale,color=c,lw=3)
        ax.scatter([0,end],[h[key][0]/scale,h[key][-1]/scale],color=c,s=35,zorder=4)
        ax.set(xlim=(0,1.4),xticks=[0,.6,1.2],ylim=limits,yticks=ticks)
        ax.text(0,1.2,label,transform=ax.transAxes,color=c,fontsize=22 if mobile else 20,va="bottom")
        ax.text(.03 if key == "flow" else 1,.9,change,transform=ax.transAxes,
                ha="left" if key == "flow" else "right",color=INK,fontsize=22 if mobile else 20,va="top")
        if mobile and i < 2:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("时间 / s" if zh else "Time / s",color=MUTED,fontsize=22 if mobile else 20,labelpad=7)
    footer=(("155 K 是理想气体结果，\n需校核实际气体物性与传热。" if mobile else
             "降温降低喷速，稳压提高耗气量。\n155 K 是理想气体结果；实际需校核物性与传热。") if zh else
            ("155 K is ideal-gas output,\nnot a real-air forecast." if mobile else
             "Cooling lowers jet speed; flow rises.\n155 K is ideal-gas output, not a real-air forecast."))
    fig.text(.07 if mobile else .075,.035 if mobile else .06,footer,color=MUTED,fontsize=19 if mobile else 20,va="bottom",linespacing=1.5)
    save(fig,output,"regulated-state",lang,mobile,
         "Same regulated ideal-air burst. Absolute tank pressure falls from 30 to 3 MPa, tank/plenum temperature from 300 to 155.38 K, mass flow rises from 18.91 to 26.27 kg/s. Constant-pressure converging nozzle thrust remains 10 kN because momentum and pressure thrust are temperature independent under these assumptions. Real-gas properties and heat transfer are omitted.")


def render_volume(output,lang,mobile):
    configure(lang);zh=lang=="zh"
    fig=plt.figure(figsize=(6,10.1) if mobile else (12,7.8),facecolor=PAPER)
    title(fig,zh,mobile,"A longer burst\nneeds more storage." if mobile else "A longer burst needs more storage.","持续时间，要由气瓶体积支付",
          "10 kN · 30 MPa to 3 MPa\nSame regulator and nozzle" if mobile else "10 kN · 30 MPa to 3 MPa · Same regulator and nozzle",
          "10 kN · 30 MPa 降至 3 MPa\n稳压器和喷口不变" if mobile else "10 kN · 30 MPa 降至 3 MPa · 稳压器和喷口不变")
    ax=fig.add_axes([.18,.335,.76,.38] if mobile else [.09,.30,.56,.43]);axis_style(ax,mobile)
    times=np.linspace(0,3,100)
    for thermal,c,ls in [("isothermal",MINT,(0,(5,4))),("adiabatic",BLUE,"-")]:
        t100=float(duration(x_end(thermal),True,thermal))
        ax.plot(times,times/t100*100,color=c,lw=3,ls=ls)
        ax.scatter([1],[100/t100],color=c,s=45,zorder=5)
    ax.set(xlim=(0,3.1),ylim=(0,250),xticks=[0,1,2,3],yticks=[0,100,200])
    ax.set_xlabel("10 kN 持续时间 / s" if zh else "Duration at 10 kN / s",color=MUTED,fontsize=22 if mobile else 20,labelpad=10)
    ax.text(0,1.1,"气瓶内体积 / L" if zh else "Tank internal volume / L",color=INK,fontsize=22,transform=ax.transAxes)
    lx=.07 if mobile else .71
    yy=.22 if mobile else .60
    fig.text(lx,yy,"每持续 1 秒" if zh else "For each second",color=INK,fontsize=23 if mobile else 22)
    fig.text(lx,yy-(.052 if mobile else .09),"绝热 75.6 L" if zh else "Adiabatic 75.6 L",color=BLUE,fontsize=23 if mobile else 22)
    fig.text(lx,yy-(.102 if mobile else .18),"等温 60.3 L" if zh else "Isothermal 60.3 L",color=MINT,fontsize=23 if mobile else 22)
    footer=(("两种理想热边界，非实际性能界限。\n等温需要补热；体积未含外壳。" if mobile else
             "两种理想热边界，并非实际性能上下界。\n等温需要持续补热；均未包含气瓶外壳。") if zh else
            ("Ideal heat boundaries only.\nIsothermal requires heat input." if mobile else
             "Ideal heat boundaries, not certified limits.\nIsothermal needs heat input; shell volume excluded."))
    fig.text(.07 if mobile else .09,.025 if mobile else .055,footer,color=MUTED,fontsize=18 if mobile else 19,va="bottom",linespacing=1.5)
    save(fig,output,"volume-budget",lang,mobile,
         "Volume scales linearly with regulated 10 kN duration, under fixed initial 30 MPa absolute/300 K and 3 MPa ideal-regulator cutoff. Adiabatic model requires 75.63 L per second; isothermal 60.30 L per second requires continuous heat input. These are alternative idealized thermal boundary conditions, not bounds for real compressed air. Tank hardware and exterior volume excluded.")


def render_cover(output):
    configure("en")
    fig=plt.figure(figsize=(12,6.3),facecolor=PAPER)
    ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,120),ylim=(0,63));ax.axis("off")
    # Schematic storage and controlled nozzle, no claimed production hardware.
    tank=FancyBboxPatch((12,21),28,20,boxstyle="round,pad=0,rounding_size=9",facecolor="#172e3d",edgecolor=BLUE,lw=1.5)
    ax.add_patch(tank)
    for x in [20,32]:ax.plot([x,x],[22,40],color=GRID,lw=2,alpha=.8)
    ax.plot([40,51],[31,31],color=BLUE,lw=4)
    ax.add_patch(FancyBboxPatch((49,26),10,10,boxstyle="round,pad=0,rounding_size=2",facecolor="#203743",edgecolor=BLUE,lw=1.5))
    ax.plot([59,67],[31,31],color=BLUE,lw=4)
    ax.add_patch(Polygon([[66,25],[77,29],[77,33],[66,37]],closed=True,facecolor="#354e5c",edgecolor=INK,lw=1.4))
    for i in range(26):
        y=(i-12.5)/5
        xs=np.linspace(77,110,180)
        ys=31+y*(1+(xs-77)/15)+.3*np.sin((xs-77)*.4+i*.5)
        ax.plot(xs,ys,color=BLUE,alpha=.10+(.45 if abs(y)<.8 else .06),lw=.9)
    ax.text(10,55,"ROADSTER 04",color=MUTED,fontsize=19,weight="bold")
    ax.text(10,48,"IMPULSE HAS A CLOCK",color=INK,fontsize=30,weight="bold")
    ax.text(26,31,"30 MPa",color=BLUE,fontsize=22,ha="center",va="center")
    ax.text(10,10,"10 kN",color=INK,fontsize=32)
    ax.text(36,10,"1.32 s",color=BLUE,fontsize=32)
    ax.text(64,11,"IDEAL-AIR EXAMPLE",color=MUTED,fontsize=16)
    fig.savefig(output/"cover.png",dpi=100,metadata={"Description":"Schematic ideal-air tank and nozzle. 10 kN for 1.32 seconds is an illustrative calculation, not Tesla performance."})
    plt.close(fig)


def summary(checks):
    result={"boundary":"Illustrative ideal-air calculation, not Tesla specifications or real-air prediction.",
        "parameters":{"gamma":GAMMA,"R_J_kg_K":R,"initial_temperature_K":TI,"initial_pressure_Pa_abs":PI,
            "regulator_pressure_Pa_abs":PR,"ambient_pressure_Pa_abs":PA,"tank_internal_volume_m3":VOLUME,"initial_thrust_N":TARGET},
        "assumptions":["Calorically perfect ideal air even at 30 MPa and down to 155 K; real-air deviations not calculated",
            "Rigid, well-mixed tank; adiabatic primary model, isothermal alternative requires heat input",
            "Ideal regulator with no pressure headroom or line losses; no post-regulation tail",
            "Ideal converging nozzle, exit equals throat, choked isentropic expansion, pressure thrust retained",
            "Different nozzle throat areas chosen to give the same initial force; not a shared-hardware comparison",
            "8 kN is an illustrative service threshold; neither a standard nor a nominal 10 kN rating"],
        "initial_gas_mass_kg":M0,"critical_pressure_ratio":CRITICAL,"checks":checks,
        "regulated_one_second":sample_at_time(1),"cases":{},"volume_per_second_L":{}}
    for thermal in ["adiabatic","isothermal"]:
        result["volume_per_second_L"][thermal]=100/float(duration(x_end(thermal),True,thermal))
        for reg in [True,False]:
            name=f"{'regulated' if reg else 'direct'}_{thermal}"
            h=history(reg,thermal);area=nozzle_area(PR if reg else PI)
            result["cases"][name]={"throat_area_m2":area,"throat_diameter_mm":1000*math.sqrt(4*area/math.pi),
                "initial_mass_flow_kg_s":float(h["flow"][0]),"end_mass_flow_kg_s":float(h["flow"][-1]),
                "end_time_s":float(h["time"][-1]),"end_thrust_N":float(h["force"][-1]),
                "end_tank_temperature_K":float(h["temperature"][-1]),"gas_remaining_kg":float(h["mass"][-1]),
                "total_impulse_until_3MPa_N_s":impulse(x_end(thermal),reg,thermal)}
    x80=force_fraction_to_mass(.8)
    result["direct_adiabatic_80_percent_service"]={"threshold_N":.8*TARGET,"duration_s":float(duration(x80,False)),
        "impulse_N_s":impulse(x80,False),"tank_pressure_Pa_abs":PI*x80**GAMMA,"gas_remaining_kg":M0*x80}
    result["isothermal_regulated_heat_input_W"]=float(nozzle(PR,TI,nozzle_area(PR))[0])*R*TI
    result["ideal_isothermal_refill_after_warming"]={}
    for label,gas_remaining in [("full_regulated_burst",M0*x_end()),("one_second_regulated",sample_at_time(1)["gas_remaining_kg"])]:
        warm_pressure=gas_remaining*R*TI/VOLUME
        primitive=lambda pressure:VOLUME*(pressure*math.log(pressure/PA)-pressure)
        work=primitive(PI)-primitive(warm_pressure)
        result["ideal_isothermal_refill_after_warming"][label]={
            "warm_residual_pressure_Pa_abs":warm_pressure,"work_J":work,"work_kWh":work/3.6e6,
            "time_at_5kW_s":work/5000,
            "boundary":"Ideal reversible compression from ambient into a tank held at 300 K; excludes actual compressor losses, electrical efficiency and finite-rate heat transfer."}
    result["illustrative_added_mass_penalty"]={"assumed_hardware_mass_kg":100.,"initial_added_mass_kg":100+M0,
        "reference_acceleration_m_s2":26.8224,"initial_force_penalty_N":(100+M0)*26.8224,
        "boundary":"Chosen 100 kg hardware budget plus modeled 100 L initial gas. Not a Tesla mass estimate or a complete variable-mass vehicle trajectory."}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parents[1]/"content"/"assets"/"roadster-cold-gas")
    args=parser.parse_args();checks=check_model();args.output_dir.mkdir(parents=True,exist_ok=True)
    for lang in ["en","zh"]:
        for mobile in [False,True]:
            for render in [render_thrust,render_state,render_volume]:render(args.output_dir,lang,mobile)
    render_cover(args.output_dir)
    result=summary(checks)
    (args.output_dir/"model.json").write_text(json.dumps(result,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=="__main__":main()
