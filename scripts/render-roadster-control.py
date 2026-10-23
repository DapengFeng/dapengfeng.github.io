#!/usr/bin/env python3
"""Reproduce Roadster 05's illustrative force allocation and delay figures.

Run with Python, NumPy, Matplotlib and Noto Sans CJK JP:
    python3 scripts/render-roadster-control.py [--output-dir DIRECTORY]

Not Tesla geometry, tyre data or a production controller. This frozen planar
snapshot groups the left and right tyre longitudinal forces at y=+b and -b.
Their independently selectable signed bounds +/-C are a relaxed envelope,
including braking. Actual motor topology, axle coupling, power, rate limits,
slip and individual tyre constraints can only further restrict this set.
All comparisons retain the same fully installed mass and wheel-force bounds.
The smaller C is a chosen remaining longitudinal capacity, not a cornering
simulation. Fixed lateral forces and their yaw moment have already been
accounted for; Mz denotes the remaining yaw-moment request in this subproblem.

Forward jet force T at y=d contributes -d*T yaw moment. The force balance is
Fx=L+R+T and Mz=b*(R-L)-d*T. The delay example integrates an uncompensated
constant moment in isolation. It is not a closed-loop vehicle trajectory.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "roadster-control-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.patches import Polygon, FancyBboxPatch, Circle
import numpy as np

B, D, T, MASS, IZZ = .8, .6, 10000., 2000., 3000.
PAPER, INK, MUTED = "#101b24", "#f2eee3", "#b5c3cc"
GRID, MINT, BLUE, AMBER = "#506778", "#bcecc0", "#8cd1f4", "#f4b780"


def allocation(fx, mz, thrust=T, offset=D, half_track=B):
    if half_track <= 0:
        raise ValueError("Half-track must be positive")
    total = fx - thrust
    difference = (mz + offset * thrust) / half_track
    return (total - difference) / 2, (total + difference) / 2


def forward_limit(capacity, mz=0., thrust=T, offset=D, half_track=B):
    if capacity < 0 or half_track <= 0:
        raise ValueError("Nonnegative capacity and positive half-track required")
    correction = abs(mz + offset * thrust) / half_track
    if correction > 2 * capacity:
        return None
    return thrust + 2 * capacity - correction


def envelope(capacity, thrust=T, offset=D):
    # Directly transform the independent L/R box; order traces its perimeter.
    wheel_vertices = np.array([[-capacity, -capacity], [capacity, -capacity],
                               [capacity, capacity], [-capacity, capacity]])
    return np.array([[l + r + thrust, B * (r-l) - offset*thrust] for l, r in wheel_vertices])


def intersect_zero_moment(vertices):
    intersections = []
    for a, c in zip(vertices, np.roll(vertices, -1, axis=0)):
        if a[1] == 0:
            intersections.append(a[0])
        if a[1] * c[1] < 0:
            intersections.append(a[0] - a[1] * (c[0] - a[0]) / (c[1] - a[1]))
    return max(intersections) if intersections else None


def delay_state(seconds):
    seconds = np.asarray(seconds)
    if np.any(seconds < 0):
        raise ValueError("Delay cannot be negative")
    acceleration = -D * T / IZZ
    return acceleration * seconds, .5 * acceleration * seconds**2


def moment_limits_at_force(fx, capacity=12000., thrust=T, offset=D):
    total = fx-thrust
    left_min, left_max = max(-capacity,total-capacity), min(capacity,total+capacity)
    if left_min > left_max:
        return None
    return (B*(total-2*left_max)-offset*thrust,
            B*(total-2*left_min)-offset*thrust)


def check_model():
    cases = []
    max_balance_error = 0.
    # The closed form is independently checked by intersecting the transformed
    # wheel-box edges with the zero-moment line, including braking solutions.
    for capacity in [1000., 3750., 5000., 8000., 12000., 20000.]:
        for thrust in [0., 3000., 10000.]:
            for offset in [-.6, 0., .6]:
                expected = forward_limit(capacity, thrust=thrust, offset=offset)
                actual = intersect_zero_moment(envelope(capacity, thrust, offset))
                assert (expected is None) == (actual is None)
                if expected is None:
                    continue
                assert math.isclose(actual, expected, abs_tol=1e-9)
                left, right = allocation(expected, 0., thrust, offset)
                assert abs(left) <= capacity + 1e-9 and abs(right) <= capacity + 1e-9
                error = max(abs(left+right+thrust-expected), abs(B*(right-left)-offset*thrust))
                max_balance_error = max(max_balance_error, error)
                assert error < 1e-9
                # A strictly larger forward request cannot remain in the box.
                too_left, too_right = allocation(expected + 1., 0., thrust, offset)
                assert max(abs(too_left), abs(too_right)) > capacity
                cases.append([capacity, thrust, offset])
    # Validate the entire diamond criterion against inverse wheel allocation.
    rng = np.random.default_rng(5)
    agreement = 0
    for _ in range(3000):
        fx, mz = rng.uniform(-40000, 45000), rng.uniform(-35000, 25000)
        left, right = allocation(fx, mz)
        inverse_feasible = max(abs(left), abs(right)) <= 12000.
        diamond_feasible = abs(fx-T) + abs(mz+D*T)/B <= 24000.
        assert inverse_feasible == diamond_feasible
        agreement += 1
    # Independently integrate r'=Mz/I, psi'=r with exact midpoint increments.
    yaw_rate, heading = 0., 0.
    steps, dt = 20000, .2/20000
    alpha = -D*T/IZZ
    for _ in range(steps):
        heading += (yaw_rate + .5*alpha*dt)*dt
        yaw_rate += alpha*dt
    target_rate, target_heading = delay_state(.2)
    error = max(abs(yaw_rate-float(target_rate)), abs(heading-float(target_heading)))
    assert error < 1e-11
    assert forward_limit(12000.) == 26500.
    assert forward_limit(8000.) == 18500.
    # Uncoordinated right-side saturation differs from a jointly reallocated
    # command. The latter sacrifices left-side force to recover moment balance.
    initial_left, initial_right = allocation(forward_limit(12000.), 0.)
    clipped_right = min(initial_right, 8000.)
    assert initial_left + clipped_right + T == 22500.
    assert B*(clipped_right-initial_left)-D*T == -3200.
    corrected_left = clipped_right-D*T/B
    assert corrected_left + clipped_right + T == 18500.
    assert B*(clipped_right-corrected_left)-D*T == 0.
    left_26,right_26 = allocation(26000.,0.)
    assert left_26 == 4250. and right_26 == 11750.
    # Independent vertical slice of the transformed actuator box validates
    # the asymmetric moment headroom at the same forward demand.
    slice_moments=[]
    vertices=envelope(12000.)
    for start,end in zip(vertices,np.roll(vertices,-1,axis=0)):
        if (start[0]-26000.)*(end[0]-26000.) < 0:
            slice_moments.append(start[1]+(26000.-start[0])*(end[1]-start[1])/(end[0]-start[0]))
    lower,upper=moment_limits_at_force(26000.)
    assert lower == -12400. and upper == 400.
    assert math.isclose(min(slice_moments),lower,abs_tol=1e-9)
    assert math.isclose(max(slice_moments),upper,abs_tol=1e-9)
    assert moment_limits_at_force(26500.)[1] == 0.
    return {"edge_intersection_cases": len(cases), "inverse_feasibility_cases": agreement,
            "fixed_forward_force_moment_slice_checked_N":26000.,
            "maximum_balance_residual_N_or_Nm": max_balance_error,
            "delay_integration_absolute_error_rad_or_rad_s": error}


def configure(lang):
    family = "DejaVu Sans" if lang == "en" else "Noto Sans CJK JP"
    findfont(FontProperties(family=family), fallback_to_default=False)
    plt.rcParams.update({"font.family": family, "font.size": 20, "svg.fonttype": "path",
                        "svg.hashsalt": f"roadster-control-v1-{lang}", "axes.unicode_minus": False,
                        "mathtext.fontset": "dejavusans", "lines.solid_capstyle": "round"})


def save(fig, output, stem, lang, mobile, description):
    name = f"{stem}-{'mobile-' if mobile else ''}{lang}"
    path = output / f"{name}.svg"
    fig.savefig(path, metadata={"Date": None, "Creator": "FENG · Matplotlib", "Description": description})
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    ET.parse(path)
    fig.savefig(Path(tempfile.gettempdir()) / f"roadster-control-{name}.png", dpi=100)
    plt.close(fig)


def line_legend(fig, entries, y, mobile):
    for i, (color, label) in enumerate(entries):
        x = .08 if mobile else .08 + i*.30
        yy = y-i*.043 if mobile else y
        fig.add_artist(plt.Line2D([x,x+.055],[yy,yy], transform=fig.transFigure, color=color,lw=4))
        fig.text(x+.075,yy,label,color=color,fontsize=23 if mobile else 20,va="center")


def style_axis(ax, mobile, grid="y"):
    ax.set_facecolor(PAPER)
    ax.grid(axis=grid,color=GRID,alpha=.32,lw=.7)
    ax.set_axisbelow(True)
    for key, spine in ax.spines.items():
        spine.set_visible(key in ("left", "bottom")); spine.set_color(GRID)
    ax.tick_params(colors=MUTED,labelsize=22 if mobile else 20,length=0,pad=8)


def render_budget(output,lang,mobile):
    configure(lang); zh=lang=="zh"
    fig=plt.figure(figsize=(5.4,13.0) if mobile else (12,8.5),facecolor=PAPER)
    fig.text(.08,.955,"偏航平衡约束下的纵向力上限" if zh else ("Forward-force limit\nwith yaw balance" if mobile else "Forward-force limit with yaw balance"),
             color=INK,fontsize=25 if mobile else 29,va="top",linespacing=1.25)
    sub=("零偏航力矩 · 同一已装车质量" if zh else "Zero yaw moment · Same installed mass")
    if mobile and not zh: sub="Zero yaw moment\nSame installed mass"
    fig.text(.08,.85 if mobile else .855,sub,color=MUTED,fontsize=21 if mobile else 20,va="top",linespacing=1.4)
    line_legend(fig,[(BLUE,"轮胎" if zh else "Tyres"),(AMBER,"喷气" if zh else "Jet")],.765 if mobile else .78,mobile)
    labels=["无喷气","居中喷气","偏置喷气"] if zh else ["No jet","Centred jet","Offset jet"]
    for i,c in enumerate([12000.,8000.]):
        rect=[.10,.45-i*.28,.82,.18] if mobile else [.22,.49-i*.32,.72,.20]
        ax=fig.add_axes(rect); style_axis(ax,mobile,"x")
        title=(f"每侧上限 {c/1000:g} kN" if zh else f"{c/1000:g} kN capacity per side")
        ax.text(0,1.22,title,transform=ax.transAxes,color=INK,fontsize=23 if mobile else 22)
        for j,(jet,offset) in enumerate([(0.,0.),(T,0.),(T,D)]):
            total=forward_limit(c,thrust=jet,offset=offset)/1000
            tyre=total-jet/1000
            y=2-j
            ax.barh(y,tyre,height=.48,color=BLUE,zorder=3)
            if jet:
                ax.barh(y,jet/1000,left=tyre,height=.48,color=AMBER,zorder=3)
            ax.text(total+.65,y,f"{total:g}",color=INK,fontsize=23 if mobile else 23,va="center")
            if mobile:
                ax.text(.0,y+.31,labels[j],color=MUTED,fontsize=21,va="bottom")
            else:
                ax.text(-1.,y,labels[j],color=MUTED,fontsize=21,va="center",ha="right")
        ax.set(xlim=(0,38.5),ylim=(-.5,2.8),xticks=[0,10,20,30],yticks=[])
        ax.spines["left"].set_visible(False)
        if i == 1:
            ax.set_xlabel("可用前向合力 / kN" if zh else "Available forward force / kN",color=MUTED,fontsize=22 if mobile else 20,labelpad=10)
    footer=("偏置 0.6 m，产生 6 kN·m 力矩\n零偏航约束减少 7.5 kN 纵向能力" if zh else "0.6 m offset: 6 kN·m moment.\nYaw balance reduces forward\ntyre-force capacity by 7.5 kN.") if mobile else (
        "0.6 m 偏置产生 6 kN·m 力矩；零偏航约束使纵向轮胎力上限减少 7.5 kN。" if zh else "0.6 m offset; 6 kN·m moment; 7.5 kN reduction in forward tyre-force capacity.")
    fig.text(.08,.025 if mobile else .035,footer,color=MUTED,fontsize=20 if mobile else 19,va="bottom",linespacing=1.45)
    save(fig,output,"force-budget",lang,mobile,
         "Illustrative fixed-state relaxed force allocation, not Tesla data. At zero residual yaw moment and per-side capacity 12 kN, total forward force is 24 without jet, 34 with centred 10 kN jet, 26.5 with jet at lateral offset 0.6 m. With per-side capacity 8 kN, totals are 16, 26, 18.5 kN. Half-track is 0.8 m; yaw correction gives up 7.5 kN of tyre force in both comparisons. Same installed mass and bounds.")


def render_envelope(output,lang,mobile):
    configure(lang); zh=lang=="zh"
    fig=plt.figure(figsize=(5.4,11.3) if mobile else (12,8.7),facecolor=PAPER)
    fig.text(.08,.955,"纵向力与偏航力矩可行域" if zh else ("Forward-force / yaw-moment\nfeasible set" if mobile else "Forward-force and yaw-moment feasible set"),
             color=INK,fontsize=24 if mobile else 29,va="top",linespacing=1.25)
    fig.text(.08,.85 if mobile else .855,("每侧 ±12 kN · 喷气 10 kN" if zh else "±12 kN per side · 10 kN jet"),color=MUTED,fontsize=21 if mobile else 20,va="top")
    line_legend(fig,[(BLUE,"居中喷气" if zh else "Centred jet"),(AMBER,"偏置 0.6 m" if zh else "Offset 0.6 m")],.77 if mobile else .775,mobile)
    ax=fig.add_axes([.19,.275,.74,.34] if mobile else [.11,.21,.62,.46]);style_axis(ax,mobile)
    for offset,color in [(0.,BLUE),(D,AMBER)]:
        points=envelope(12000.,offset=offset)/1000
        ax.add_patch(Polygon(points,closed=True,facecolor=color,edgecolor=color,lw=2.5,alpha=.11))
        closed=np.vstack([points,points[0]])
        ax.plot(closed[:,0],closed[:,1],color=color,lw=2.5)
    ax.axhline(0,color=INK,lw=1.2,ls=(0,(4,5)),alpha=.8)
    ax.scatter([26.5,34],[0,0],s=70,c=[AMBER,BLUE],edgecolors=PAPER,zorder=6)
    ax.set(xlim=(-17,39),ylim=(-29,23),xticks=[-10,0,10,20,30],yticks=[-20,0,20])
    ax.set_xlabel("前向合力 / kN" if zh else "Forward force / kN",color=MUTED,fontsize=22,labelpad=10)
    ax.text(0,1.10,"偏航力矩 / kN·m" if zh else "Yaw moment / kN·m",transform=ax.transAxes,color=INK,fontsize=22)
    if mobile:
        fig.text(.08,.17,"力矩为零时" if zh else "At zero yaw moment",color=INK,fontsize=23)
        fig.text(.08,.125,"26.5 kN",color=AMBER,fontsize=26)
        fig.text(.55,.125,"34 kN",color=BLUE,fontsize=26)
        fig.text(.08,.025,"阴影内可分配\n边界上已有一侧饱和\n负前向力允许制动" if zh else "Inside: feasible.\nBoundary: a side saturates.\nSigned forces allow braking.",color=MUTED,fontsize=19,va="bottom",linespacing=1.45)
    else:
        fig.text(.78,.59,"零力矩截面" if zh else "$M_z=0$",color=INK,fontsize=22)
        fig.text(.78,.49,"26.5 kN",color=AMBER,fontsize=29)
        fig.text(.78,.39,"34 kN",color=BLUE,fontsize=29)
        fig.text(.08,.035,"阴影内可分配，边界上至少一侧饱和。\n负前向力允许制动；实际执行器约束可能进一步缩小可行域。" if zh else "Inside is feasible; a boundary saturates a side.\nBraking is allowed; actual actuator limits may shrink the set.",
                 color=MUTED,fontsize=18,va="bottom")
    save(fig,output,"force-envelope",lang,mobile,
         "Feasible total forward-force and residual yaw-moment diamonds for independent signed left/right bounds ±12 kN, half-track 0.8 m and a fixed 10 kN forward jet. A centred jet leaves the moment centre at zero; a jet at y=0.6 m translates it to -6 kNm. Their zero-yaw maximum forward forces are 34 and 26.5 kN. This relaxed envelope includes braking but omits real motor topology, slip, power and rate constraints.")


def render_delay(output,lang,mobile):
    configure(lang);zh=lang=="zh"
    fig=plt.figure(figsize=(5.4,11.3) if mobile else (12,8.0),facecolor=PAPER)
    fig.text(.08,.955,"补偿延迟引起的偏航响应" if zh else ("Yaw response to\ncompensation delay" if mobile else "Yaw response to compensation delay"),
             color=INK,fontsize=24 if mobile else 29,va="top",linespacing=1.25)
    subtitle=("未抵消力矩 −6 kN·m\n$I_z$ = 3,000 kg·m²" if zh else "Unmatched moment −6 kN·m\n$I_z$ = 3,000 kg·m²") if mobile else (
        "未抵消力矩 −6 kN·m · $I_z$ = 3,000 kg·m²" if zh else "Unmatched moment −6 kN·m · $I_z$ = 3,000 kg·m²")
    fig.text(.08,.85 if mobile else .84,subtitle,color=MUTED,fontsize=21 if mobile else 20,va="top",linespacing=1.45)
    times=np.linspace(0,.2,201);samples=np.array([.05,.1,.2])
    rates,angles=delay_state(times);sample_rates,sample_angles=delay_state(samples)
    for i,(curve,points,color,label,limits,ticks) in enumerate([
        (np.abs(rates)*180/np.pi,np.abs(sample_rates)*180/np.pi,BLUE,"偏航角速度增量 / °/s" if zh else "Yaw-rate change / °/s",(0,27),[0,10,20]),
        (np.abs(angles)*180/np.pi,np.abs(sample_angles)*180/np.pi,AMBER,"航向角增量 / °" if zh else "Heading increment / °",(0,2.7),[0,1,2])]):
        rect=[.18,.51-i*.285,.73,.19] if mobile else [.09+i*.48,.27,.37,.36]
        ax=fig.add_axes(rect);style_axis(ax,mobile)
        ax.plot(times*1000,curve,color=color,lw=3)
        ax.fill_between(times*1000,curve,0,color=color,alpha=.08)
        ax.scatter(samples*1000,points,color=color,s=45,zorder=5)
        for x,y in zip(samples*1000,points):
            ax.annotate(f"{y:.2f}",(x,y),xytext=(-9,10) if x==200 else (0,11),textcoords="offset points",
                        ha="right" if x==200 else "center",color=color,fontsize=21 if mobile else 20)
        ax.set(xlim=(0,215),ylim=limits,xticks=[0,100,200],yticks=ticks)
        ax.text(0,1.18,label,transform=ax.transAxes,color=color,fontsize=21 if mobile else 21)
        if not mobile or i == 1:
            ax.set_xlabel("未补偿时间 / ms" if zh else "Uncompensated time / ms",color=MUTED,fontsize=22 if mobile else 20,labelpad=9)
    footer=("显示增量的绝对值；初始增量为零\n只积分力矩，不预测车辆轨迹" if zh else "Magnitudes; zero initial increments.\nMoment integration only,\nnot a vehicle trajectory.") if mobile else (
        "显示增量的绝对值，初始增量为零。\n仅积分恒定未抵消力矩；未模拟轮胎反馈或闭环车辆轨迹。" if zh else "Magnitudes, from zero initial increments. Constant-moment integration only;\ntyre feedback and closed-loop motion omitted.")
    fig.text(.08,.025 if mobile else .055,footer,color=MUTED,fontsize=19 if mobile else 18,va="bottom",linespacing=1.5)
    save(fig,output,"delay-budget",lang,mobile,
         "Illustrative open-loop moment accumulation: uncompensated -6 kNm, constant yaw inertia 3000 kg m2, zero initial yaw-rate and heading increments. After 50,100,200 ms, yaw-rate magnitudes are 5.73,11.46,22.92 degrees/s and heading magnitudes .143,.573,2.292 degrees. No tyre response or closed-loop trajectory is modeled; these values are not safety thresholds.")


def render_cover(output):
    configure("en")
    fig=plt.figure(figsize=(12,6.3),facecolor=PAPER)
    ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,120),ylim=(0,63));ax.axis("off")
    ax.text(8,54,"ROADSTER 05",color=MUTED,fontsize=18,weight="bold")
    ax.text(8,45,"FORCE ALLOCATION AND YAW CONTROL",color=INK,fontsize=29,weight="bold")
    # Geometric force/moment sets are the subject; no invented production car.
    for x,y,color in [(76,23,BLUE),(76,16,AMBER)]:
        ax.add_patch(Polygon([[x-24,y],[x,y+15],[x+24,y],[x,y-15]],closed=True,
                            edgecolor=color,facecolor=color,alpha=.13,lw=2))
        ax.plot([x-24,x,x+24,x,x-24],[y,y+15,y,y-15,y],color=color,lw=2)
    ax.plot([47,106],[23,23],color=MUTED,ls=(0,(4,5)),lw=1)
    ax.scatter([88.8,100],[23,23],color=[AMBER,BLUE],s=70,zorder=5)
    ax.text(8,27,"10 kN",color=AMBER,fontsize=34)
    ax.text(8,20,"jet thrust",color=MUTED,fontsize=19)
    ax.text(8,10,"2.5 kN",color=INK,fontsize=34)
    ax.text(8,3.5,"forward-force gain at zero yaw moment",color=MUTED,fontsize=17)
    fig.savefig(output/"cover.png",dpi=100,metadata={"Description":"Roadster 05 illustrative force allocation: a 10 kN offset jet adds only 2.5 kN of forward capacity when tyres cancel its yaw moment. Not Tesla data."})
    plt.close(fig)


def render_geometry(output, mobile):
    """Shared-language geometric schematic; vector lengths show direction only."""
    configure("en")
    fig=plt.figure(figsize=(5.4,8.2) if mobile else (12,6.5),facecolor=PAPER)
    ax=fig.add_axes([.03,.14,.94,.81] if mobile else [0,.13,.67,.80])
    ax.set(xlim=(-1.75,1.75),ylim=(-2.15,2.05),aspect="equal");ax.axis("off")
    ax.add_patch(FancyBboxPatch((-.65,-1.45),1.3,2.9,boxstyle="round,pad=0,rounding_size=.45",
                              facecolor="#1c303e",edgecolor=GRID,lw=1.5))
    ax.add_patch(FancyBboxPatch((-.43,.45),.86,.55,boxstyle="round,pad=0,rounding_size=.18",
                              facecolor=PAPER,edgecolor=GRID,lw=1))
    for x in [-B,B]:
        for y in [-.9,.9]:
            ax.add_patch(FancyBboxPatch((x-.10,y-.30),.20,.60,boxstyle="round,pad=0,rounding_size=.055",
                                      facecolor="#233d4c",edgecolor=MUTED,lw=1.5))
    ax.plot([0,0],[-1.55,1.55],color=MUTED,lw=.9,ls=(0,(3,6)),alpha=.8)
    # Screen left is body +y; upward is body +x. Actual line-of-action offsets
    # are preserved while the stylized vehicle outline is not a dimensional drawing.
    for x,label in [(-B,"L"),(B,"R")]:
        ax.annotate("",(x,.49),(x,-.38),arrowprops={"arrowstyle":"-|>","color":BLUE,"lw":3,"mutation_scale":22})
        ax.text(x+(-.23 if x<0 else .23),.28,f"${label}$",color=BLUE,fontsize=30,
                ha="right" if x<0 else "left",va="center")
    ax.annotate("",(-D,-.15),(-D,-1.30),arrowprops={"arrowstyle":"-|>","color":AMBER,"lw":3,"mutation_scale":22})
    ax.text(-D+.18,-.75,"$T$",color=AMBER,fontsize=30,va="center")
    ax.add_patch(Circle((0,0),.055,facecolor=INK,edgecolor=PAPER,lw=1,zorder=9))
    ax.text(.12,.05,"$O$",color=INK,fontsize=23,va="bottom")
    # Dimension guides convey the two independent lateral offsets.
    for x in [-B,0,B]:ax.plot([x,x],[-1.64,-1.90],color=GRID,lw=1)
    for lo,hi in [(-B,0),(0,B)]:
        ax.plot([lo,hi],[-1.78,-1.78],color=MUTED,lw=1.1)
        ax.text((lo+hi)/2,-1.96,"$b$",color=INK,fontsize=24,ha="center",va="center")
    ax.plot([-D,-D],[1.5,1.85],color=AMBER,lw=1,ls=(0,(3,4)),alpha=.8)
    ax.plot([0,0],[1.5,1.85],color=MUTED,lw=1)
    ax.plot([-D,0],[1.70,1.70],color=AMBER,lw=1.1)
    ax.text(-D/2,1.88,"$d$",color=AMBER,fontsize=24,ha="center",va="center")
    ax.annotate("",(1.38,1.88),(1.38,1.40),arrowprops={"arrowstyle":"-|>","color":MUTED,"lw":1.3,"mutation_scale":12})
    ax.annotate("",(.98,1.40),(1.38,1.40),arrowprops={"arrowstyle":"-|>","color":MUTED,"lw":1.3,"mutation_scale":12})
    ax.text(1.40,1.97,"$+x$",color=MUTED,fontsize=21,ha="center")
    ax.text(1.02,1.59,"$+y$",color=MUTED,fontsize=21,ha="center",va="center")
    if mobile:
        fig.text(.5,.095,r"$b=0.8\,\mathrm{m}\qquad d=0.6\,\mathrm{m}$",ha="center",color=MUTED,fontsize=22)
        fig.text(.5,.035,"$M_z=b(R-L)-dT$",ha="center",color=INK,fontsize=24)
    else:
        fig.text(.65,.70,r"$b=0.8\,\mathrm{m}$",color=MUTED,fontsize=24)
        fig.text(.65,.61,r"$d=0.6\,\mathrm{m}$",color=AMBER,fontsize=24)
        fig.text(.65,.44,"$F_x=L+R+T$",color=INK,fontsize=24)
        fig.text(.65,.33,"$M_z=b(R-L)-dT$",color=INK,fontsize=24)
    stem="force-geometry-mobile" if mobile else "force-geometry"
    description="Shared top-view schematic; x forward/up the drawing, y left. Left and right grouped forward forces act at y=+b and -b, b=0.8 m. The jet acts at y=d=0.6 m. Point O is the centre of mass. Force arrows show direction, not magnitude; the stylized body is not Tesla geometry."
    path=output/f"{stem}.svg"
    fig.savefig(path,metadata={"Date":None,"Creator":"FENG · Matplotlib","Description":description})
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    ET.parse(path)
    fig.savefig(Path(tempfile.gettempdir())/f"roadster-control-{stem}.png",dpi=100)
    plt.close(fig)


def summary(checks):
    result={"boundary":"Illustrative frozen planar force-allocation problem; not Tesla data or a production controller.",
            "parameters":{"half_track_m":B,"jet_lateral_offset_m":D,"forward_jet_N":T,
                          "same_fully_installed_mass_kg":MASS,"assumed_yaw_inertia_kg_m2":IZZ},
            "sign_convention":"x forward, y left, z up; left wheel group y=+b, right=-b; forward jet y=d gives yaw moment -dT.",
            "equations":{"forward_force":"Fx=L+R+T","yaw_moment":"Mz=b(R-L)-dT",
                         "feasible_set":"abs(Fx-T)+abs(Mz+dT)/b <= 2C",
                         "maximum_forward_force":"Fx_max=T+2C-abs(Mz+dT)/b, provided abs(Mz+dT)<=2bC"},
            "assumptions":["Independent signed side forces in [-C,C] permit braking; optimistic relaxed envelope, not a claim about Tesla motor topology",
                           "Chosen C=12 or 8 kN is residual longitudinal capacity after any prescribed lateral forces; not a tyre measurement or corner simulation",
                           "Prescribed lateral-force yaw moments are already accounted for; Mz is the remaining moment request",
                           "Loads, geometry, capacity and jet thrust frozen; all comparisons retain the same installed mass",
                           "No axle coupling, power/speed limits, force slew, slip or wheel dynamics in the static envelope",
                           "Delay integrates constant unmatched moment with constant inertia and zero initial increments; not a closed-loop trajectory",
                           "Restoring zero net moment stops the incremental angular acceleration, not the already accumulated yaw rate"],
            "checks":checks,"snapshots":[],"delay_samples":[]}
    for c in [12000.,8000.]:
        cases=[]
        for name,thrust,offset in [("no_jet",0.,0.),("centred_jet",T,0.),("offset_jet",T,D)]:
            force=forward_limit(c,thrust=thrust,offset=offset)
            left,right=allocation(force,0.,thrust,offset)
            cases.append({"name":name,"left_force_N":left,"right_force_N":right,"jet_force_N":thrust,
                          "total_forward_N":force,"net_yaw_Nm":B*(right-left)-offset*thrust,
                          "instantaneous_acceleration_m_s2_if_no_other_longitudinal_force":force/MASS})
        result["snapshots"].append({"residual_capacity_per_side_N":c,"cases":cases})
    for t in [.05,.1,.2]:
        rate,heading=delay_state(t)
        result["delay_samples"].append({"delay_s":t,"delta_yaw_rate_rad_s":float(rate),"delta_heading_rad":float(heading),
                                       "delta_yaw_rate_deg_s":float(rate*180/np.pi),"delta_heading_deg":float(heading*180/np.pi)})
    result["asymmetric_right_capacity_reduction"]={
        "new_right_cap_N":8000., "unchanged_left_cap_N":12000.,
        "clipped_only":{"left_N":4500.,"right_N":8000.,"jet_N":T,"forward_N":22500.,"yaw_Nm":-3200.},
        "jointly_reallocated":{"left_N":500.,"right_N":8000.,"jet_N":T,"forward_N":18500.,"yaw_Nm":0.}}
    left,right=allocation(26000.,0.)
    lower,upper=moment_limits_at_force(26000.)
    result["fixed_26kN_operating_point"]={"forward_N":26000.,"yaw_Nm":0.,"left_N":left,"right_N":right,
        "right_unused_forward_N":12000.-right,"at_fixed_forward_force_yaw_min_Nm":lower,
        "at_fixed_forward_force_yaw_max_Nm":upper,
        "boundary":"Moment increments measured from zero moment while total forward force remains 26 kN; signed independent side forces bounded ±12 kN."}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parents[1]/"content"/"assets"/"roadster-control")
    args=parser.parse_args();checks=check_model();args.output_dir.mkdir(parents=True,exist_ok=True)
    for lang in ["en","zh"]:
        for mobile in [False,True]:
            for render in [render_budget,render_envelope,render_delay]:render(args.output_dir,lang,mobile)
    for mobile in [False,True]:render_geometry(args.output_dir,mobile)
    render_cover(args.output_dir)
    result=summary(checks)
    (args.output_dir/"model.json").write_text(json.dumps(result,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
    print(json.dumps({"checks":checks,"snapshots":result["snapshots"],"delay_samples":result["delay_samples"]},indent=2))


if __name__=="__main__":main()
