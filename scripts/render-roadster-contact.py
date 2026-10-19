#!/usr/bin/env python3
"""Reproduce the Roadster 03 contact-load diagrams and illustrative bounds.

Requires Python 3, NumPy, Matplotlib and Noto Sans CJK for Chinese labels.
Run: python3 scripts/render-roadster-contact.py [--output-dir DIRECTORY]

This is an illustrative quasi-static force budget, NOT Tesla specifications,
measured tyre data, a handling prediction or a suspension simulation. The
concave power law is an explicit chosen sensitivity, not a fitted tyre model.
All forces act on a flat road. Four normal loads follow vertical, pitch and
roll balances; lambda prescribes the front share of lateral load transfer.
D is imposed normal force, with beta at the front axle. No device mass, drag,
power draw, suspension compliance, tyre temperature, slip-angle coupling,
transient response or actuator delays are included.

The corner bound assumes an equal front/rear lateral demand for a centered CG
and independently usable tyre peak lateral forces, with zero longitudinal
force. Each axle is evaluated at the candidate acceleration, not at a frozen
load split. Roots stop before wheel unloading. Unused axle curves are drawn
only where all four normal loads remain nonnegative. They do not imply that
steering, tyre slip and vehicle dynamics can realize independent tyre peaks.

Assets include six en/zh SVG pairs with narrow-screen editions, a PNG cover,
and model.json. PNG inspection previews are written only to /tmp.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "roadster-contact-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.patches import FancyBboxPatch, Circle, PathPatch
from matplotlib.path import Path as PlotPath
import numpy as np

MASS, G, WHEELBASE, TRACK, HEIGHT = 2000., 9.81, 2.8, 1.7, .5
D, BETA, LAMBDA = 4905., .5, .6
FREF, MU_REF, EXPONENT = 4905., 1.5, .88
PAPER, INK, MUTED = "#111a20", "#f2eee3", "#b6c1c0"
GRID, MINT, BLUE, AMBER = "#60737a", "#bcecc0", "#97d5f5", "#f4bb83"


def parameters(*, downforce=D, beta=BETA, share=LAMBDA, p=EXPONENT, height=HEIGHT):
    if not all(math.isfinite(x) for x in (downforce, beta, share, p, height)):
        raise ValueError("Finite model parameters required")
    if height <= 0 or downforce < 0 or not 0 <= beta <= 1 or not 0 <= share <= 1 or not 0 < p <= 1:
        raise ValueError("h > 0, D >= 0, beta/lambda in [0,1], p in (0,1] required")


def normal_loads(ax=0., ay=0., *, downforce=D, beta=BETA, share=LAMBDA, height=HEIGHT):
    """Order: front outside, front inside, rear outside, rear inside; ay >= 0."""
    parameters(downforce=downforce, beta=beta, share=share, height=height)
    if not math.isfinite(ax) or not math.isfinite(ay) or ay < 0:
        raise ValueError("Finite accelerations and ay >= 0 required")
    front = MASS * G / 2 + beta * downforce - MASS * ax * height / WHEELBASE
    rear = MASS * G + downforce - front
    transfer = MASS * ay * height / TRACK  # amount shifted, NOT outer-minus-inner
    result = np.array([front / 2 + share * transfer, front / 2 - share * transfer,
                       rear / 2 + (1 - share) * transfer, rear / 2 - (1 - share) * transfer])
    if result.min() < -1e-8:
        raise ValueError("Wheel unloading: four-contact model no longer applies")
    return np.maximum(result, 0.)


def tyre_capacity(load, p=EXPONENT):
    parameters(p=p)
    value = np.asarray(load, dtype=float)
    if not np.all(np.isfinite(value)) or np.any(value < 0):
        raise ValueError("Finite nonnegative normal loads required")
    return MU_REF * FREF * (value / FREF) ** p


def axle_loads(ay, axle, *, downforce=D, beta=BETA, share=LAMBDA, height=HEIGHT):
    """Individual-axle inputs, including before the other axle's lift cutoff."""
    normal = MASS * G / 2 + (beta if axle == 0 else 1 - beta) * downforce
    fraction = share if axle == 0 else 1 - share
    shifted = fraction * MASS * ay * height / TRACK
    return normal / 2 + shifted, normal / 2 - shifted


def lift_limit(axle, *, downforce=D, beta=BETA, share=LAMBDA, height=HEIGHT):
    fraction = share if axle == 0 else 1 - share
    if fraction == 0:
        return math.inf
    normal = MASS * G / 2 + (beta if axle == 0 else 1 - beta) * downforce
    return normal * TRACK / (2 * fraction * MASS * height)


def capacity_residual(ay, axle, *, downforce=D, beta=BETA, share=LAMBDA, p=EXPONENT, height=HEIGHT):
    loads = np.maximum(axle_loads(ay, axle, downforce=downforce, beta=beta, share=share, height=height), 0.)
    return float(tyre_capacity(loads, p).sum() - MASS * ay / 2)


def axle_bound(axle, *, downforce=D, beta=BETA, share=LAMBDA, p=EXPONENT, height=HEIGHT):
    parameters(downforce=downforce, beta=beta, share=share, p=p, height=height)
    if axle not in (0, 1):
        raise ValueError("Axle must be 0 or 1")
    opts = dict(downforce=downforce, beta=beta, share=share, p=p, height=height)
    contact = lift_limit(axle, downforce=downforce, beta=beta, share=share, height=height)
    hi = min(contact, 100 * G)
    if capacity_residual(hi, axle, **opts) > 0:
        return {"ay_m_s2": contact, "constraint": "wheel-unloading"}
    lo = 0.
    for _ in range(80):
        mid = (lo + hi) / 2
        if capacity_residual(mid, axle, **opts) > 0:
            lo = mid
        else:
            hi = mid
    return {"ay_m_s2": (lo + hi) / 2, "constraint": "tyre-force"}


def corner_bound(beta=BETA, *, downforce=D, share=LAMBDA, p=EXPONENT, height=HEIGHT):
    opts = dict(downforce=downforce, beta=beta, share=share, p=p, height=height)
    bounds = [axle_bound(a, **opts) for a in (0, 1)]
    first = min(range(2), key=lambda a: bounds[a]["ay_m_s2"])
    ay = bounds[first]["ay_m_s2"]
    loads = normal_loads(ay=ay, downforce=downforce, beta=beta, share=share, height=height)
    return {"beta": beta, "ay_m_s2": ay, "ay_g": ay / G,
            "limiting_axle": "front" if first == 0 else "rear",
            "constraint": bounds[first]["constraint"], "axle_bounds": bounds,
            "normal_loads_N": loads.tolist(), "axle_force_demand_N": MASS * ay / 2,
            "axle_capacity_N": [float(tyre_capacity(loads[:2], p).sum()), float(tyre_capacity(loads[2:], p).sum())]}


def optimal_balance(*, downforce=D, share=LAMBDA, p=EXPONENT, height=HEIGHT):
    lo, hi = .1, .9
    for _ in range(65):
        mid = (lo + hi) / 2
        b = [axle_bound(a, beta=mid, downforce=downforce, share=share, p=p, height=height)["ay_m_s2"] for a in (0, 1)]
        if b[0] < b[1]:
            lo = mid
        else:
            hi = mid
    return corner_bound((lo + hi) / 2, downforce=downforce, share=share, p=p, height=height)


MIN_LOAD_EXAMPLE = 500.0  # Illustrative screening floor, NOT a safety standard.
BRAKE_FRONT_SHARE = .6
CONTACT_NAMES = ["front-outer", "front-inner", "rear-outer", "rear-inner"]


def root_between(fn, lo, hi):
    """Bisection with an explicit sign-changing bracket."""
    flo, fhi = fn(lo), fn(hi)
    if abs(flo) < 1e-10: return lo
    if abs(fhi) < 1e-10: return hi
    if flo * fhi > 0: raise ValueError("Root is not bracketed")
    for _ in range(70):
        mid = (lo + hi) / 2
        fm = fn(mid)
        if fm * flo > 0:
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def screened_beta_range(*, target_ay=1.65*G, p=EXPONENT, floor=MIN_LOAD_EXAMPLE):
    """Pure-lateral axle force bounds plus four illustrative normal-load floors."""
    def margins(beta):
        f = normal_loads(ay=target_ay, beta=beta)
        c = tyre_capacity(f, p)
        return [float(c[:2].sum()-MASS*target_ay/2), float(c[2:].sum()-MASS*target_ay/2),
                float(f[1]-floor), float(f[3]-floor)]
    # First use the analytic load-floor bounds so no negative loads are sampled.
    transfer = MASS * target_ay * HEIGHT / TRACK
    low_load = (2*(floor + LAMBDA*transfer) - MASS*G/2) / D
    high_load = 1-(2*(floor + (1-LAMBDA)*transfer) - MASS*G/2) / D
    lo, hi = max(0.,low_load), min(1.,high_load)
    if lo > hi: raise ValueError("No beta satisfies the normal-load floors")
    if margins(lo)[0] < 0:
        lo = root_between(lambda b:margins(b)[0],lo,hi)
    if margins(hi)[1] < 0:
        hi = root_between(lambda b:margins(b)[1],lo,hi)
    labels=["front-axle-force","rear-axle-force","front-inner-normal-floor","rear-inner-normal-floor"]
    return {"p":p,"target_ay_g":target_ay/G,"illustrative_minimum_load_N":floor,
        "beta_min":lo,"beta_max":hi,
        "lower_constraint":labels[int(np.argmin(np.abs(margins(lo))))],
        "upper_constraint":labels[int(np.argmin(np.abs(margins(hi))))],
        "lower_margins_N":margins(lo),"upper_margins_N":margins(hi)}


def design_comparison():
    optimum=optimal_balance()
    cases=[("baseline",.5,.6,.5),("aero-balance",optimum["beta"],.6,.5),
           ("lower-cg",.5,.6,.45),("equal-roll-share",.5,.5,.5)]
    result=[]
    for name,beta,share,height in cases:
        b=corner_bound(beta,share=share,height=height)
        floor_limits=[]
        for axle in [0,1]:
            n=MASS*G/2+(beta if axle==0 else 1-beta)*D
            fraction=share if axle==0 else 1-share
            floor_limits.append((n/2-MIN_LOAD_EXAMPLE)*TRACK/(fraction*MASS*height))
        result.append({"id":name,"beta":beta,"lambda":share,"height_m":height,
            "peak_ay_g":b["ay_g"],"minimum_load_at_own_peak_N":min(b["normal_loads_N"]),
            "peak_normal_loads_N":b["normal_loads_N"],
            "screened_ay_g":min(b["ay_m_s2"],*floor_limits)/G,
            "screening_floor_N":MIN_LOAD_EXAMPLE})
    return result


def combined_budget(braking_g=.6, lateral_g=1., *, downforce=D, beta=BETA, gamma=BRAKE_FRONT_SHARE,
                    share=LAMBDA, height=HEIGHT, p=EXPONENT):
    """Explicit prescribed allocation, NOT an optimizer or real control law.

    Front/rear brake totals are gamma/(1-gamma) times m|ax|; each axle splits
    Fx equally left/right. Front/rear Fy totals are both m*ay/2; within each
    axle, Fy is proportional to the chosen peak capacity C_i. This enforces
    zero net yaw in vehicle axes but does NOT prove steering/slip compatibility.
    """
    if braking_g < 0 or lateral_g < 0 or not 0 <= gamma <= 1:
        raise ValueError("Nonnegative acceleration magnitudes, gamma in [0,1] required")
    fz=normal_loads(-braking_g*G,lateral_g*G,downforce=downforce,beta=beta,share=share,height=height)
    c=tyre_capacity(fz,p)
    total_braking=MASS*braking_g*G
    fx=np.array([-gamma*total_braking/2]*2+[-(1-gamma)*total_braking/2]*2)
    fy=np.r_[c[:2]/c[:2].sum(),c[2:]/c[2:].sum()]*(MASS*lateral_g*G/2)
    resultant=np.hypot(fx,fy)
    use=np.divide(resultant,c,out=np.full(4,np.inf),where=c>0)
    # x: front positive; y: inside positive for this chosen positive-Fy turn.
    xs=np.array([WHEELBASE/2]*2+[-WHEELBASE/2]*2)
    ys=np.array([-TRACK/2,TRACK/2]*2)
    mz=float(np.sum(xs*fy-ys*fx))
    first=int(np.argmax(use))
    return {"braking_g":braking_g,"lateral_g":lateral_g,"downforce_N":downforce,"beta":beta,"gamma":gamma,
        "normal_loads_N":fz.tolist(),"capacity_N":c.tolist(),"Fx_N":fx.tolist(),"Fy_N":fy.tolist(),
        "resultant_N":resultant.tolist(),"utilisation":use.tolist(),"maximum_utilisation":float(use[first]),
        "most_used_contact":CONTACT_NAMES[first],"sum_Fx_N":float(fx.sum()),"sum_Fy_N":float(fy.sum()),
        "net_yaw_moment_Nm":mz,"minimum_normal_load_N":float(min(fz)),
        "within_rule_budget":bool(max(use)<=1+1e-12 and min(fz)>=MIN_LOAD_EXAMPLE)}


def disturbance_boundaries():
    bmax=root_between(lambda b:combined_budget(b)["maximum_utilisation"]-1,.6,1.2)
    dmin=root_between(lambda d:combined_budget(downforce=d)["maximum_utilisation"]-1,0,D)
    blo=root_between(lambda b:combined_budget(beta=b)["maximum_utilisation"]-1,0,.5)
    bhi=root_between(lambda b:combined_budget(beta=b)["maximum_utilisation"]-1,.5,1)
    # At this fixed ax/ay and allocation, axle requirements depend separately
    # on front and rear downforce. These two constants reproduce the 2-D boundary.
    front_needed=blo*D;rear_needed=(1-bhi)*D
    return {"operating_point":combined_budget(),"braking_boundary":combined_budget(bmax),
        "downforce_boundary":combined_budget(downforce=dmin),
        "beta_lower_boundary":combined_budget(beta=blo),"beta_upper_boundary":combined_budget(beta=bhi),
        "required_front_downforce_N":front_needed,"required_rear_downforce_N":rear_needed,
        "domain":{"downforce_N":[0,6000],"beta":[0,1]},
        "force_constraint_is_first_in_plotted_domain":True,
        "downforce_loss_cases":[{"lost_fraction":f,**combined_budget(downforce=D*(1-f))} for f in [0,.2,.5,.6]],
        "braking_sweep":[combined_budget(float(b)) for b in np.linspace(0,1.,101)]}


def check_extended_model():
    for design in design_comparison():
        assert design["screened_ay_g"]<=design["peak_ay_g"]+1e-10
        assert design["minimum_load_at_own_peak_N"]>0
    for p in [.84,.88,.92]:
        interval=screened_beta_range(p=p)
        assert 0<interval["beta_min"]<interval["beta_max"]<1
        assert min(interval["lower_margins_N"])>-1e-7
        assert min(interval["upper_margins_N"])>-1e-7
    for brake in [0,.3,.6,.9]:
        for d in [0,D,6000]:
            for beta in [0,.5,1]:
                b=combined_budget(brake,downforce=d,beta=beta)
                assert math.isclose(b["sum_Fx_N"],-MASS*brake*G,abs_tol=1e-8)
                assert math.isclose(b["sum_Fy_N"],MASS*G,abs_tol=1e-8)
                assert abs(b["net_yaw_moment_Nm"])<1e-8
                assert math.isclose(sum(b["normal_loads_N"]),MASS*G+d,abs_tol=1e-8)
    bounds=disturbance_boundaries()
    for k in ["braking_boundary","downforce_boundary","beta_lower_boundary","beta_upper_boundary"]:
        assert math.isclose(bounds[k]["maximum_utilisation"],1.,abs_tol=1e-10)
        assert bounds[k]["minimum_normal_load_N"]>MIN_LOAD_EXAMPLE
    dmin=bounds["downforce_boundary"]["downforce_N"]
    assert combined_budget(downforce=dmin+1)["within_rule_budget"]
    assert not combined_budget(downforce=dmin-1)["within_rule_budget"]
    for beta in np.linspace(.05,.95,19):
        needed=max(bounds["required_front_downforce_N"]/beta,bounds["required_rear_downforce_N"]/(1-beta))
        if needed<6000:
            state=combined_budget(downforce=needed,beta=float(beta))
            assert abs(state["maximum_utilisation"]-1)<1e-9
            assert state["minimum_normal_load_N"]>MIN_LOAD_EXAMPLE
    # Independent reference numbers computed by a separate four-wheel audit.
    assert np.allclose(combined_budget()["utilisation"],[.5404128103,.7797291972,.6908709818,.8408927658],atol=1e-9)
    assert math.isclose(dmin,2276.284754,abs_tol=1e-5)


def check_model():
    # Force and moment balances; same sign conventions as the article.
    for ax in [-.6 * G, 0., .6 * G]:
        for ay in [0., .4 * G, G]:
            for beta in [.2, .5, .8]:
                f = normal_loads(ax, ay, beta=beta)
                assert math.isclose(sum(f), MASS * G + D, abs_tol=1e-8)
                pitch = (sum(f[2:]) - sum(f[:2])) * WHEELBASE / 2 - (1 - 2 * beta) * D * WHEELBASE / 2
                assert math.isclose(pitch, MASS * ax * HEIGHT, abs_tol=1e-8)
                roll = (f[0] + f[2] - f[1] - f[3]) * TRACK / 2
                assert math.isclose(roll, MASS * ay * HEIGHT, abs_tol=1e-8)
                assert math.isclose((f[0] - f[1]) * TRACK / 2, LAMBDA * MASS * ay * HEIGHT, abs_tol=1e-8)
    assert np.allclose(normal_loads(downforce=0), [FREF] * 4)
    assert float(tyre_capacity(0)) == 0.
    assert np.allclose(tyre_capacity([0, FREF, 2 * FREF], p=1), np.array([0, FREF, 2 * FREF]) * MU_REF)
    assert float(tyre_capacity([1962, 7848]).sum()) < float(tyre_capacity([FREF, FREF]).sum())
    assert math.isclose(float(tyre_capacity([1962, 7848], 1).sum()), 2 * MU_REF * FREF)
    # A symmetric load-transfer distribution has an optimum at beta = .5.
    symmetric = corner_bound(share=.5)
    assert math.isclose(symmetric["axle_bounds"][0]["ay_m_s2"], symmetric["axle_bounds"][1]["ay_m_s2"], abs_tol=1e-9)
    for beta in np.linspace(.1, .9, 41):
        b = corner_bound(float(beta))
        assert b["constraint"] == "tyre-force"
        assert min(b["normal_loads_N"]) > 0
        # Independent root check: explicit per-contact arithmetic and secant solve.
        a = 0 if b["limiting_axle"] == "front" else 1
        n = MASS * G / 2 + (beta if a == 0 else 1 - beta) * D
        fraction = LAMBDA if a == 0 else 1 - LAMBDA
        def independent(y):
            z = fraction * MASS * y * HEIGHT / TRACK
            return MU_REF * FREF ** (1 - EXPONENT) * ((n / 2 + z) ** EXPONENT + (n / 2 - z) ** EXPONENT) - MASS * y / 2
        x0, x1 = b["ay_m_s2"] * .99, b["ay_m_s2"] * 1.001
        for _ in range(8):
            q0, q1 = independent(x0), independent(x1)
            if abs(q1) < 1e-8:
                break
            x0, x1 = x1, x1 - q1 * (x1 - x0) / (q1 - q0)
        assert math.isclose(x1, b["ay_m_s2"], abs_tol=1e-8)
        assert independent(b["ay_m_s2"] * .999) > 0
        assert independent(b["ay_m_s2"] * 1.001) < 0
        # More prescribed downforce at the same balance cannot reduce this bound.
        assert corner_bound(float(beta), downforce=D / 2)["ay_m_s2"] < b["ay_m_s2"]
    best = optimal_balance()
    assert abs(best["axle_capacity_N"][0] - best["axle_capacity_N"][1]) < 1e-7
    for beta in [.1, .3, .5, .7, .9]:
        assert corner_bound(beta)["ay_m_s2"] <= best["ay_m_s2"] + 1e-9
    for fn in [lambda: tyre_capacity(-1), lambda: tyre_capacity(100, 1.2), lambda: normal_loads(ay=-1),
               lambda: normal_loads(ay=3 * G), lambda: corner_bound(beta=1.1), lambda: corner_bound(downforce=-1)]:
        try:
            fn()
        except ValueError:
            pass
        else:
            raise AssertionError("Expected invalid model input to fail")
    assert math.isclose(math.sqrt(8000**2 - 6000**2), 5291.502622129181, rel_tol=1e-12)


def configure(language):
    family = "DejaVu Sans" if language == "en" else "Noto Sans CJK JP"
    findfont(FontProperties(family=family), fallback_to_default=False)
    plt.rcParams.update({"font.family": family, "font.size": 20, "svg.fonttype": "path",
        "svg.hashsalt": f"roadster-contact-v1-{language}", "axes.unicode_minus": False,
        "mathtext.fontset": "dejavusans", "lines.solid_capstyle": "round"})


def canvas(language, mobile, title, subtitle, *, mobile_height=10, mobile_subtitle_size=19):
    configure(language)
    fig = plt.figure(figsize=(6, mobile_height) if mobile else (12, 7.5), facecolor=PAPER)
    fig.text(.07 if mobile else .055, .96, title, color=INK, fontsize=25 if mobile else 29, va="top", linespacing=1.35)
    fig.text(.07 if mobile else .055, .86 if mobile else .865, subtitle, color=MUTED, fontsize=mobile_subtitle_size if mobile else 20, va="top", linespacing=1.4)
    return fig


def axis_style(ax, mobile=False):
    ax.set_facecolor(PAPER)
    ax.grid(color=GRID, alpha=.25, lw=.7)
    ax.set_axisbelow(True)
    for where, spine in ax.spines.items():
        spine.set_visible(where in ("left", "bottom"))
        spine.set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=19 if mobile else 20, length=0, pad=8)


def save(fig, output, stem, lang, mobile, description):
    name = f"{stem}-{'mobile-' if mobile else ''}{lang}"
    path = output / f"{name}.svg"
    fig.savefig(path, metadata={"Date": None, "Creator": "FENG · Matplotlib", "Description": description})
    path.write_text("\n".join(x.rstrip() for x in path.read_text().splitlines()) + "\n")
    ET.parse(path)
    fig.savefig(Path(tempfile.gettempdir()) / f"roadster-contact-{name}.png", dpi=100)
    plt.close(fig)


def render_loads(output, lang, mobile):
    """One spatial map plus a shared-scale load sequence; tyre shapes are schematic."""
    zh = lang == "zh"
    configure(lang)
    fig = plt.figure(figsize=(6, 13.5) if mobile else (12, 8), facecolor=PAPER)
    left = .07 if mobile else .055
    fig.text(left, .966, "总量相同，四轮各不相同" if zh else
             ("One total.\nFour different loads." if mobile else "One total. Four different loads."),
             color=INK, fontsize=26 if mobile else 29, va="top", linespacing=1.25)
    eq_y = (.877 if not zh else .911) if mobile else .864
    fig.text(left, eq_y, "19.620 + 4.905 = 24.525 kN", color=MINT,
             fontsize=23 if mobile else 26, va="top")
    fig.text(left, eq_y - (.038 if mobile else .063),
             "重力 + 下压力 = 四轮法向载荷" if zh else "Weight + downforce = total normal load",
             color=MUTED, fontsize=18 if mobile else 19, va="top")

    # A neutral vehicle outline, not a claimed Roadster body or contact-area model.
    car = fig.add_axes([.07, .465, .86, .35] if mobile else [.035, .17, .40, .59])
    car.set(xlim=(0, 100), ylim=(0, 112)); car.set_aspect("equal"); car.axis("off")
    def outline(vertices, codes, **kw):
        car.add_patch(PathPatch(PlotPath(vertices, codes), **kw))
    v = [(50,100),(38,100),(32,94),(31,82),(30,69),(34,56),(32,37),
         (30,23),(31,10),(39,7),(50,7),(61,7),(69,10),(70,23),
         (68,37),(66,56),(70,69),(69,82),(68,94),(62,100),(50,100)]
    # Quadratic spline keeps the silhouette restrained and smooth.
    outline(v, [PlotPath.MOVETO]+[PlotPath.CURVE3]*20,
            facecolor="#1a272f", edgecolor="#73868d", lw=1.25)
    car.plot([50,50],[10,99], color=GRID, alpha=.30, lw=.8, dashes=(3,5))
    for y in [78,29]:
        car.plot([29,71],[y,y],color=GRID,alpha=.55,lw=1)
    # Windscreen, cabin and rear glass make orientation readable without ornament.
    outline([(36,69),(50,76),(64,69),(60,56),(50,60),(40,56),(36,69)],
            [PlotPath.MOVETO,PlotPath.CURVE3,PlotPath.CURVE3,PlotPath.LINETO,PlotPath.CURVE3,PlotPath.CURVE3,PlotPath.CLOSEPOLY],
            facecolor="#253943",edgecolor="#6a8592",lw=.9)
    outline([(40,51),(50,54),(60,51),(61,36),(50,32),(39,36),(40,51)],
            [PlotPath.MOVETO,PlotPath.CURVE3,PlotPath.CURVE3,PlotPath.LINETO,PlotPath.CURVE3,PlotPath.CURVE3,PlotPath.CLOSEPOLY],
            facecolor="#142027",edgecolor="#506770",lw=.85)
    car.plot([37,63],[88,88],color=BLUE,alpha=.6,lw=1.1)
    car.plot([38,62],[17,17],color=AMBER,alpha=.6,lw=1.1)
    car.scatter([50],[54],s=18,facecolor=PAPER,edgecolor=MUTED,linewidth=.8,zorder=5)
    car.text(50,108,"车头" if zh else "FRONT",ha="center",va="center",
             fontsize=19 if mobile else 16,color=MUTED)
    values = normal_loads(-.6 * G, G) / 1000
    names = ["前外", "前内", "后外", "后内"] if zh else ["Front outer", "Front inner", "Rear outer", "Rear inner"]
    for j, value in enumerate(values):
        front, outer = j < 2, j % 2 == 0
        c = BLUE if front else AMBER
        x, y = (25 if outer else 69), (70 if front else 21)
        car.add_patch(FancyBboxPatch((x,y),6,16,boxstyle="round,pad=0,rounding_size=1.6",
                      facecolor="#10171d",edgecolor=c,lw=1.5))
        for dx in [1.9,4.1]:
            car.plot([x+dx,x+dx],[y+3,y+13],color=c,alpha=.35,lw=.65)
        tx, ha = (2,"left") if outer else (98,"right")
        ty = 96 if front else 28
        car.text(tx,ty,names[j] if zh else names[j].replace(" ", "\n"),
                 ha=ha,va="bottom",color=c,fontsize=21 if mobile else 16,linespacing=1.1)
        car.text(tx,ty-10,f"{value:.2f}",ha=ha,va="center",color=INK,fontsize=26 if mobile else 24)
        # Leaders point to the tyre, while numbers stay outside the silhouette.
        xs = [3,17,25] if outer else [97,83,75]
        ys = [ty-17,ty-17,y+8]
        car.plot(xs,ys,color=c,alpha=.65,lw=.9,ls="-" if outer else (0,(3,2)))
    car.text(50,-1,"制动 0.6g · 转弯 1g" if zh else "Brake 0.6g · corner 1g",
             ha="center",color=MUTED,fontsize=19 if mobile else 17)

    # Every stage has exactly the same total load; only its distribution changes.
    stages = np.array([normal_loads(), normal_loads(-.6*G), normal_loads(-.6*G,G)]) / 1000
    assert np.allclose(stages.sum(axis=1), (MASS*G+D)/1000)
    ax = fig.add_axes([.145, .125, .70, .265] if mobile else [.51, .25, .37, .44])
    ax.set_facecolor(PAPER)
    ax.set(xlim=(-.06,2.30),ylim=(0,12))
    for side in ax.spines.values(): side.set_visible(False)
    ax.set_yticks([0,4,8,12]); ax.set_xticks([0,1,2])
    ax.set_xticklabels((["静态","制动","制动+转弯"] if zh else ["Static","Brake","Brake +\ncorner"]),
                       fontsize=20 if mobile else 18)
    ax.tick_params(axis="both",colors=MUTED,length=0,pad=9,labelsize=20 if mobile else 18)
    ax.grid(axis="y",color=GRID,alpha=.25,lw=.7)
    for x in [0,1,2]: ax.axvline(x,color=GRID,alpha=.20,lw=.7)
    ax.text(0,1.12,"每轮载荷 / kN" if zh else "Normal load / kN",transform=ax.transAxes,
            color=INK,fontsize=22 if mobile else 22)
    ax.text(0,1.035,"下压力固定为 4.905 kN" if zh else "Downforce held at 4.905 kN",
            transform=ax.transAxes,color=MUTED,fontsize=18)
    # Draw inner traces last so dashed segments remain visible while paths coincide.
    for j in [0,2,1,3]:
        c = BLUE if j < 2 else AMBER
        if j % 2:
            ax.plot([0,1,2],stages[:,j],color=PAPER,lw=4.5,zorder=3)
        ax.plot([0,1,2],stages[:,j],color=c,lw=2.4,ls="-" if j%2==0 else (0,(4,3)),zorder=4)
        ax.scatter([2],[stages[2,j]],s=28,facecolor=c if j%2==0 else PAPER,edgecolor=c,lw=1.2,zorder=5)
        ax.text(2.10,stages[2,j],f"{stages[2,j]:.2f}",va="center",color=c,fontsize=22 if mobile else 20)
    ax.scatter([0,1,1],[stages[0,0],stages[1,0],stages[1,2]],s=22,facecolor=PAPER,edgecolor=MUTED,zorder=6,lw=1)
    ax.annotate("6.13",(0,stages[0,0]),xytext=(0,-24),textcoords="offset points",ha="center",color=MUTED,fontsize=19)
    # Line style and position carry the same meaning as the spatial map.
    legend_y = .062 if mobile else .135
    for k,(style,label) in enumerate([("-","外侧" if zh else "Outer"),((0,(4,3)),"内侧" if zh else "Inner")]):
        lx = (.15 if mobile else .51) + k*(.39 if mobile else .21)
        line = plt.Line2D([lx,lx+.065],[legend_y,legend_y],transform=fig.transFigure,color=MUTED,lw=2,ls=style)
        fig.add_artist(line)
        fig.text(lx+.085,legend_y,label,color=MUTED,fontsize=20 if mobile else 18,va="center")
    fig.text(left,.022 if mobile else .041,
             "β = 0.5 · λ = 0.6 · 准静态算例" if zh else "β = 0.5 · λ = 0.6 · Quasi-static example",
             color=MUTED,fontsize=18 if mobile else 17,va="bottom")
    save(fig, output, "four-contact-loads", lang, mobile,
         "Normal loads in kN, not contact-patch areas. All curve stages hold downforce at 4905 N: static, braking at 0.6g, then braking while cornering at 1g. Outer traces are solid; inner traces dashed. The top view shows the final state. These are prescribed accelerations, not measured Tesla performance.")


def render_sensitivity(output, lang, mobile):
    zh = lang == "zh"
    title = "载荷加倍，抓地力未必加倍" if zh else "Double the load,\nless than double the grip" if mobile else "Double the load, less than double the grip"
    fig = canvas(lang, mobile, title, "同一车轴：总载荷始终为 9.81 kN" if zh else "Axle normal load stays 9.81 kN")
    ax = fig.add_axes([.16, .40, .76, .36] if mobile else [.105, .39, .84, .33])
    axis_style(ax, mobile)
    loads = np.linspace(0, 11000, 300)
    ax.plot(loads / 1000, MU_REF * loads / 1000, color=GRID, ls=(0, (3, 4)), lw=1.8)
    ax.plot(loads / 1000, tyre_capacity(loads) / 1000, color=MINT, lw=3)
    for load in [1962, 4905, 7848]:
        value = float(tyre_capacity(load)) / 1000
        ax.plot([load / 1000] * 2, [0, value], color=BLUE if load == FREF else AMBER, lw=1, alpha=.7)
        ax.scatter(load / 1000, value, s=62, color=BLUE if load == FREF else AMBER, zorder=4)
    ax.set_xlim(0, 11); ax.set_ylim(0, 17)
    ax.set_xticks([0, 5, 10]); ax.set_yticks([0, 5, 10, 15])
    ax.set_xlabel("单胎载荷 / kN" if zh else "Normal load per tyre / kN", color=MUTED, fontsize=19, labelpad=10)
    ax.text(0, 1.075, "单胎力上限 / kN" if zh else "Tyre force capacity / kN", color=MUTED, fontsize=19, transform=ax.transAxes)
    ax.text(.045, .90, "p = 1" if mobile else "p = 1  ·  固定 μ" if zh else "p = 1  ·  constant μ", transform=ax.transAxes, color=MUTED, fontsize=18)
    ax.text(.58, .22, "p = 0.88", transform=ax.transAxes, color=MINT, fontsize=20)
    entries = [("50 : 50", "4.905 + 4.905", "14.715", BLUE), ("80 : 20", "7.848 + 1.962", "14.411", AMBER)]
    if mobile:
        for i, (split, normal, total, color) in enumerate(entries):
            y = .255 - .095 * i
            fig.text(.075, y, split, fontsize=21, color=color)
            fig.text(.40, y, f"{total} kN", fontsize=24, color=color)
            fig.text(.075, y - .035, ("法向载荷 " if zh else "Normal loads ") + normal + " kN", fontsize=18, color=MUTED)
        fig.text(.075, .055, "同样总载荷，力上限少 2.06%" if zh else "Same total load; 2.06% less force", color=INK, fontsize=19)
    else:
        for i, (split, normal, total, color) in enumerate(entries):
            x = .105 + i * .46
            fig.text(x, .23, f"{split}     {total} kN", color=color, fontsize=25)
            fig.text(x, .16, ("法向载荷 " if zh else "Normal loads ") + normal + " kN", fontsize=18, color=MUTED)
        fig.text(.105, .055, "相同总载荷：80 : 20 分配损失 2.06% 的车轴力上限" if zh else "Same total load: the 80 : 20 split loses 2.06% of axle force capacity", color=INK, fontsize=20)
    save(fig, output, "load-sensitivity", lang, mobile, "Chosen illustrative capacity law C=1.5*4905*(Fz/4905)^0.88. It is not fitted Tesla tyre data. Equal and 80:20 loads have the same 9810 N axle total.")


def render_balance(output, lang, mobile):
    zh=lang=="zh"
    fig=canvas(lang,mobile,"先筛出一个有余量的区间" if zh else "Screen a range,\nthen test the assumptions" if mobile else "Screen a range, then test the assumptions",
        "$a_y$ = 1.65g · D = 4.905 kN\n$F_z$ ≥ 500 N" if mobile else "$a_y$ = 1.65g · D = 4.905 kN · $F_z$ ≥ 500 N",
        mobile_height=11.5, mobile_subtitle_size=22)
    ax=fig.add_axes([.26,.43,.65,.30] if mobile else [.11,.38,.82,.36])
    axis_style(ax,mobile)
    intervals=[screened_beta_range(p=p) for p in [.84,.88,.92]]
    shared_lo=max(x["beta_min"] for x in intervals);shared_hi=min(x["beta_max"] for x in intervals)
    ax.axvspan(shared_lo,shared_hi,color=MINT,alpha=.09)
    ax.axvline(.5,color=MUTED,ls=(0,(2,4)),lw=1.4)
    for i,row in enumerate(intervals):
        y=2-i;lo,hi=row["beta_min"],row["beta_max"]
        ax.plot([lo,hi],[y,y],color=[BLUE,MINT,AMBER][i],lw=7)
        ax.scatter([lo,hi],[y,y],color=[BLUE,MINT,AMBER][i],s=50,zorder=4)
        ax.text(hi,y+.22,f"{100*hi:.1f}%",color=[BLUE,MINT,AMBER][i],ha="center",fontsize=22 if mobile else 20)
    ax.set_xlim(.45,.8);ax.set_ylim(-.5,2.6)
    ax.set_xticks([.5,.6,.7,.8]);ax.set_yticks([2,1,0],["p = .84","p = .88","p = .92"])
    ax.tick_params(axis="y",labelsize=22 if mobile else 18,pad=7)
    if mobile:
        ax.tick_params(axis="x",labelsize=22)
    ax.set_xlabel("前轴下压力占比 β" if zh else "Front downforce fraction β",color=MUTED,fontsize=22 if mobile else 19,labelpad=10)
    if mobile:
        fig.text(.075,.31,"53.3–66.8%",fontsize=29,color=MINT)
        fig.text(.075,.265,"三种假设都满足的范围" if zh else "Common to all three p values",fontsize=22,color=MUTED)
        fig.text(.075,.205,"下界：内前轮载荷\n上界：后轴力" if zh else "Lower: inner-front load\nUpper: rear-axle force",fontsize=22,color=MUTED,linespacing=1.5,va="top")
        fig.text(.075,.085,"500 N 为演示门槛\n并非安全标准" if zh else "500 N: illustrative floor,\nnot a safety standard",fontsize=22,color=INK,linespacing=1.4,va="top")
    else:
        fig.text(.11,.19,"53.3–66.8%",fontsize=28,color=MINT)
        fig.text(.41,.19,"三种假设都满足的范围" if zh else "Passes all three assumed p values",fontsize=20,color=MUTED)
        fig.text(.11,.115,"下界：内前轮载荷；上界：后轴力" if zh else "Lower boundary: inner-front load · Upper boundary: rear-axle force",fontsize=19,color=MUTED)
        fig.text(.11,.045,"500 N 仅为演示门槛，非安全标准" if zh else "The 500 N load floor is illustrative, not a safety standard",fontsize=19,color=INK)
    save(fig,output,"aero-balance",lang,mobile,"Pure-lateral quasi-static force budget at 1.65g, prescribed D=4905 N and minimum normal load 500 N. The three chosen p values are scenarios, not statistical confidence intervals. Lower beta bound comes from inside-front 500 N load; upper from rear-axle capacity. No vehicle safety standard or handling prediction.")


def render_designs(output,lang,mobile):
    zh=lang=="zh"
    fig=canvas(lang,mobile,"相同下压力，比较四种方案" if zh else "Same downforce.\nFour design choices." if mobile else "Same downforce. Four design choices.",
        "D = 4.905 kN · p = 0.88")
    names=(["基线","调整下压力分布","降低质心","平分侧向转移"] if zh else ["Baseline","Rebalance downforce","Lower the CG","Equal roll-transfer share"])
    data=design_comparison()
    if mobile:
        for i,(name,row) in enumerate(zip(names,data)):
            y=.76-i*.145
            fig.text(.075,y,name,fontsize=19,color=INK)
            fig.text(.925,y,f"{row['peak_ay_g']:.3f}g",fontsize=25,color=BLUE,ha="right")
            params=f"β {row['beta']:.3f} · λ {row['lambda']:.2f} · h {row['height_m']:.2f} m"
            fig.text(.075,y-.04,params,fontsize=18,color=MUTED)
            fig.text(.075,y-.08,("峰值处最小轮载 " if zh else "Min load at its peak ")+f"{row['minimum_load_at_own_peak_N']:.0f} N",fontsize=18,color=AMBER if i<2 else MINT)
            line=fig.add_axes([.075,y-.105,.85,.008]);line.set_axis_off()
            line.set_xlim(0,1200);line.plot([0,1200],[0,0],color=GRID,lw=3,alpha=.3)
            line.plot([0,row['minimum_load_at_own_peak_N']],[0,0],color=AMBER if i<2 else MINT,lw=3)
            line.axvline(500,color=INK,lw=1)
        fig.text(.075,.105,"每行在各自的横向加速度上界取值" if zh else "Each row is evaluated at its own peak",fontsize=18,color=MUTED)
        fig.text(.075,.055,"白线：500 N 演示筛选值" if zh else "White tick: illustrative 500 N floor",fontsize=18,color=INK)
    else:
        fig.text(.48,.76,"加速度上界" if zh else "Acceleration bound",fontsize=19,color=BLUE)
        fig.text(.73,.76,"峰值处最小轮载" if zh else "Min load at own peak",fontsize=19,color=MINT)
        for i,(name,row) in enumerate(zip(names,data)):
            y=.66-i*.14
            fig.text(.07,y,name,fontsize=22,color=INK)
            fig.text(.07,y-.052,f"β {row['beta']:.3f} · λ {row['lambda']:.2f} · h {row['height_m']:.2f} m",fontsize=18,color=MUTED)
            fig.text(.48,y,f"{row['peak_ay_g']:.3f}g",fontsize=28,color=BLUE)
            fig.text(.73,y,f"{row['minimum_load_at_own_peak_N']:.0f} N",fontsize=27,color=AMBER if i<2 else MINT)
            bar=fig.add_axes([.73,y-.05,.21,.012]);bar.set_axis_off();bar.set_xlim(0,1200)
            bar.plot([0,1200],[0,0],color=GRID,lw=4,alpha=.35)
            bar.plot([0,row['minimum_load_at_own_peak_N']],[0,0],color=AMBER if i<2 else MINT,lw=4)
            bar.axvline(500,color=INK,lw=1.2)
        fig.text(.07,.065,"每行在各自峰值取值；白线是 500 N 演示筛选值" if zh else "Each row uses its own peak; the white tick is an illustrative 500 N floor",fontsize=19,color=MUTED)
    save(fig,output,"design-tradeoffs",lang,mobile,"Four changes with the same prescribed downforce and tyre law. Peak lateral acceleration upper bound and minimum normal load at each design's OWN peak are shown; they are not compared at one common acceleration. 500 N is an illustrative filter, not a safety standard. No implementation cost, slip compatibility or vehicle measurement is claimed.")


def render_combined(output,lang,mobile):
    zh=lang=="zh";state=combined_budget();limits=disturbance_boundaries();limit=limits["braking_boundary"]
    fig=canvas(lang,mobile,"每个接地点，都要闭合合力预算" if zh else "Close the force budget\nat every contact" if mobile else "Close the force budget at every contact",
        "$a_x$ = −0.6g · $a_y$ = 1.0g\nγ = 0.60" if mobile else "$a_x$ = −0.6g · $a_y$ = 1.0g · γ = 0.60",
        mobile_height=11.5, mobile_subtitle_size=22)
    names=(["前外","前内","后外","后内"] if zh else ["Front outer","Front inner","Rear outer","Rear inner"])
    for i,name in enumerate(names):
        row,col=divmod(i,2)
        ax=fig.add_axes([.07+.46*col,.51-.265*row,.40,.255] if mobile else [.07+.46*col,.48-.27*row,.40,.24])
        ax.set_axis_off();ax.set_xlim(0,1);ax.set_ylim(0,1)
        color=BLUE if i<2 else AMBER
        ax.plot([0,1],[.99,.99],color=color,lw=2)
        ax.text(0,.80,name,color=color,fontsize=22 if mobile else 21)
        if mobile:
            for j,(label,key) in enumerate([(r"$F_z$","normal_loads_N"),(r"$F_x$","Fx_N"),(r"$F_y$","Fy_N")]):
                y=.64-j*.18
                ax.text(0,y,label,color=MUTED,fontsize=22)
                ax.text(.98,y,f"{state[key][i]/1000:.2f} kN",ha="right",color=INK,fontsize=22)
            ax.text(0,.10,"u",color=MUTED,fontsize=22)
            ax.text(.98,.10,f"{state['utilisation'][i]:.2f}",ha="right",color=color,fontsize=25)
            ax.plot([0,1],[.015,.015],color=GRID,lw=3,alpha=.35)
            ax.plot([0,state['utilisation'][i]],[.015,.015],color=color,lw=3)
            ax.plot([1,1],[.005,.035],color=INK,lw=1)
        else:
            for j,(label,key) in enumerate([(r"$F_z$","normal_loads_N"),(r"$F_x$","Fx_N"),(r"$F_y$","Fy_N")]):
                x=j*.34
                ax.text(x,.57,label+" / kN",color=MUTED,fontsize=18)
                ax.text(x,.32,f"{state[key][i]/1000:.2f}",color=INK,fontsize=24)
            ax.text(0,.02,f"u = {state['utilisation'][i]:.2f}",color=color,fontsize=23)
            ax.plot([.38,1],[.08,.08],color=GRID,lw=4,alpha=.35)
            ax.plot([.38,.38+.62*state['utilisation'][i]],[.08,.08],color=color,lw=4)
    if mobile:
        fig.text(.075,.22,"制动增至 0.754g\n后内轮先满" if zh else "At 0.754g braking,\nrear inner fills first",fontsize=22,color=AMBER,linespacing=1.15,va="top")
        fig.text(.075,.15,"ΣFx = −11.772 kN\nΣFy = 19.620 kN · Mz = 0",fontsize=22,color=MUTED,linespacing=1.15,va="top")
        fig.text(.075,.08,"固定分配规则\n并非控制器最优解" if zh else "Fixed allocation,\nnot a control optimum",fontsize=22,color=INK,linespacing=1.15,va="top")
    else:
        fig.text(.07,.14,"制动增至 0.754g：后内轮首先 u = 1" if zh else "Raise braking to 0.754g: the rear inner tyre reaches u = 1 first",fontsize=21,color=AMBER)
        fig.text(.07,.075,"ΣFx = −11.772 kN · ΣFy = 19.620 kN · Mz = 0",fontsize=20,color=MUTED)
        fig.text(.07,.025,"固定分配规则；不是控制器最优解" if zh else "Prescribed force allocation, not a control optimum",fontsize=18,color=INK)
    save(fig,output,"combined-contact-budget",lang,mobile,"Fixed front braking share gamma=.6; equal Fx left/right per axle. Each axle gets half the total lateral force, split between its tyres in proportion to each C_i. This closes total Fx, Fy and net yaw moment. u=hypot(Fx,Fy)/C. It does not prove tyre-slip compatibility, real control behavior or an optimal allocation.")


def render_pressure(output,lang,mobile):
    zh=lang=="zh";bounds=disturbance_boundaries()
    fig=canvas(lang,mobile,("这套分配能容忍\n多少下压力损失？" if mobile else "这套分配能容忍多少下压力损失？") if zh else "How much downforce\ncan this allocation lose?" if mobile else "How much downforce can this allocation lose?",
        "$a_x$ = −0.6g · $a_y$ = 1.0g\nγ = 0.60" if mobile else "$a_x$ = −0.6g · $a_y$ = 1.0g · γ = 0.60",
        mobile_height=11.5, mobile_subtitle_size=22)
    ax=fig.add_axes([.16,.43,.76,.29] if mobile else [.105,.38,.84,.33]);axis_style(ax,mobile)
    if mobile:
        ax.tick_params(labelsize=22)
    beta=np.linspace(.001,.999,800)
    f=bounds['required_front_downforce_N']/beta/1000
    r=bounds['required_rear_downforce_N']/(1-beta)/1000
    need=np.maximum(f,r)
    ax.fill_between(beta,need,6,where=need<=6,color=MINT,alpha=.10)
    ax.plot(beta,np.where(f>=r,f,np.nan),color=BLUE,lw=3)
    ax.plot(beta,np.where(r>=f,r,np.nan),color=AMBER,lw=3)
    blo=bounds['beta_lower_boundary']['beta'];bhi=bounds['beta_upper_boundary']['beta'];dmin=bounds['downforce_boundary']['downforce_N']/1000
    ax.plot([blo,bhi],[D/1000]*2,color=INK,lw=1.4,ls=(0,(3,4)))
    ax.plot([.5,.5],[dmin,D/1000],color=INK,lw=1.3,ls=(0,(3,4)))
    ax.scatter([blo,.5,bhi],[D/1000]*3,color=INK,s=[35,55,35],zorder=5)
    ax.scatter(.5,dmin,color=AMBER,s=55,zorder=5)
    ax.text(.48 if mobile else .51,D/1000+.18,"基线" if zh else "Baseline",fontsize=22 if mobile else 18,color=INK,ha="right" if mobile else "left")
    ax.text(.54,dmin-.30,"2.276",fontsize=22 if mobile else 19,color=AMBER)
    ax.text(.28,3.5,"允许区" if zh else "Feasible",color=MINT,fontsize=22 if mobile else 20)
    ax.set_xlim(0,1);ax.set_ylim(0,6)
    ax.set_xticks([0,.25,.5,.75,1]);ax.set_yticks([0,2,4,6])
    ax.set_xlabel("前轴下压力占比 β" if zh else "Front downforce fraction β",color=MUTED,fontsize=22 if mobile else 19,labelpad=10)
    ax.text(0,1.09,"总下压力 D / kN" if zh else "Total downforce D / kN",transform=ax.transAxes,color=MUTED,fontsize=22 if mobile else 19)
    if mobile:
        fig.text(.075,.325,"内前轮先满" if zh else "Front inner first",color=BLUE,fontsize=22)
        fig.text(.54,.325,"内后轮先满" if zh else "Rear inner first",color=AMBER,fontsize=22)
        fig.text(.075,.265,"β = 0.5 : D ≥ 2.276 kN",color=INK,fontsize=22)
        fig.text(.075,.22,"D = 4.905 kN:\nβ = 0.059–0.768",color=MUTED,fontsize=22,linespacing=1.3,va="top")
        fig.text(.075,.115,"本图首限为轮胎合力\n准静态规则边界\n并非故障响应" if zh else "Tyre-force limits bind first.\nStatic allocation budget,\nnot a fault response.",color=INK,fontsize=22,linespacing=1.3,va="top")
    else:
        fig.text(.105,.205,"内前轮先满" if zh else "Front inner first",color=BLUE,fontsize=19)
        fig.text(.43,.205,"内后轮先满" if zh else "Rear inner first",color=AMBER,fontsize=19)
        fig.text(.105,.125,"β = 0.5 : D ≥ 2.276 kN",color=INK,fontsize=24)
        fig.text(.57,.125,"D = 4.905 kN : β 0.059–0.768",color=MUTED,fontsize=19)
        fig.text(.105,.045,"本图首限为轮胎合力；准静态规则边界不代表故障响应" if zh else "Tyre force binds first; this static allocation boundary is not a fault response",color=MUTED,fontsize=19)
    save(fig,output,"pressure-reserve",lang,mobile,"Fixed braking .6g and cornering 1g under the specified gamma=.6 allocation. Force feasibility requires beta*D >= 290.480416 N and (1-beta)*D >= 1138.142377 N. All normal loads in this domain exceed 500 N; the load-floor constraint is inactive. The domain is not control-independent and has no time-response meaning.")


def render_cover(output):
    configure("en")
    fig = plt.figure(figsize=(12, 6.3), facecolor=PAPER)
    fig.text(.065, .88, "ROADSTER / 03", fontsize=21, color=MUTED)
    fig.text(.065, .72, "FOUR CONTACTS.", fontsize=39, color=INK)
    fig.text(.065, .60, "ONE FORCE BUDGET.", fontsize=39, color=INK)
    ax = fig.add_axes([.065, .13, .51, .31]); ax.set_axis_off()
    vals = normal_loads(-.6 * G, G)
    for j, val in enumerate(vals):
        x, y = (j % 2) * .55, .58 - (j // 2) * .55
        color = BLUE if j < 2 else AMBER
        ax.add_patch(FancyBboxPatch((x, y), .43, .38, boxstyle="round,pad=.01,rounding_size=.04", facecolor=color, alpha=.25, edgecolor=color))
        ax.text(x + .215, y + .19, f"{val/1000:.2f} kN", color=INK, ha="center", va="center", fontsize=22)
    ax = fig.add_axes([.66, .18, .28, .53]); ax.set_aspect("equal"); ax.set_axis_off()
    theta = np.linspace(0, 2*np.pi, 360)
    ax.plot(np.cos(theta), np.sin(theta), color=MINT, lw=2)
    ax.axhline(0, color=GRID, lw=.8); ax.axvline(0, color=GRID, lw=.8)
    ax.plot([0, .75], [0, math.sqrt(1-.75**2)], color=INK, lw=2.5)
    ax.plot([.75, .75], [0, math.sqrt(1-.75**2)], color=AMBER, lw=4)
    ax.scatter(.75, math.sqrt(1-.75**2), color=INK, s=65)
    fig.savefig(output / "cover.png", dpi=100, metadata={"Software": "FENG · Matplotlib", "Description": "Illustrative four-contact load and force budget; not Tesla specifications."})
    plt.close(fig)


def summary():
    examples = []
    for name, d, ax, ay in [("static-no-downforce",0.,0.,0.),("static-with-downforce",D,0.,0.),("braking-and-cornering",D,-.6*G,G)]:
        loads = normal_loads(ax, ay, downforce=d)
        examples.append({"id":name,"downforce_N":d,"ax_m_s2":ax,"ay_m_s2":ay,
            "normal_loads_N":loads.tolist(),"total_normal_load_N":float(loads.sum()),
            "sum_of_independent_peak_capacities_N":float(tyre_capacity(loads).sum())})
    equal, unequal = float(tyre_capacity([FREF,FREF]).sum()), float(tyre_capacity([1962,7848]).sum())
    return {"kind":"illustrative-four-contact-force-budget-not-tesla-data-or-handling-prediction",
        "design_comparison":design_comparison(),
        "screened_balance_ranges":[screened_beta_range(p=p) for p in [.84,.88,.92]],
        "combined_allocation_rule":{"braking_front_fraction_gamma":BRAKE_FRONT_SHARE,"Fx":"equal left/right within each axle; fixed front/rear fraction gamma","Fy":"equal front/rear axle totals; within axle proportional to chosen tyre capacity C_i","scope":"constructed force allocation, not a control optimum, slip-compatible handling solution or Tesla specification","minimum_load_screen_N":MIN_LOAD_EXAMPLE,"minimum_load_screen_is_safety_standard":False},
        "disturbance_boundaries":disturbance_boundaries(),
        "assumptions":{"mass_kg":MASS,"g_m_s2":G,"wheelbase_m":WHEELBASE,"equal_front_rear_track_m":TRACK,
            "cg_height_m":HEIGHT,"front_static_weight_fraction":.5,"prescribed_downforce_N":D,
            "front_downforce_fraction_beta":BETA,"front_lateral_transfer_share_lambda":LAMBDA,
            "tyre_reference_load_N":FREF,"tyre_reference_mu":MU_REF,"chosen_capacity_exponent_p":EXPONENT,
            "load_order":["front-outer","front-inner","rear-outer","rear-inner"],
            "ax_positive":"accelerating-forward","ay_positive":"cornering; outer side receives transfer",
            "lateral_transfer_definition":"delta=m*ay*h/track is load shifted, not outside-minus-inside",
            "scope":["flat-road quasi-static balance","prescribed aero distribution; no additional aero moments",
                "pure-lateral bound: ax=Fx=0 and equal front/rear Fy demand for centered CG",
                "independent individual-tyre peak allocation assumed; steering/slip coupling omitted",
                "no device mass, drag, auxiliary power, camber, compliance, temperature or transient model",
                "all usable bounds end before any wheel unloads; near-zero loads expose model limits",
                "chosen load-sensitivity power law is not fitted data","circle assumes isotropic combined-force capacity"]},
        "snapshots":examples,
        "load_sensitivity_example":{"total_axle_normal_load_N":9810.,"equal_loads_N":[4905.,4905.],
            "unequal_loads_N":[7848.,1962.],"equal_capacity_N":equal,"unequal_capacity_N":unequal,
            "relative_capacity_loss":1-unequal/equal},
        "cornering_force_bounds":{"no_downforce":corner_bound(downforce=0),"best_beta":optimal_balance(),
            "samples":[corner_bound(b) for b in [.1,.3,.5,.56,.7,.9]],
            "beta_sweep":[corner_bound(float(b)) for b in np.linspace(.1,.9,81)]},
        "friction_circle_example":{"fixed_capacity_N":8000.,"lateral_force_N":6000.,
            "remaining_longitudinal_force_N":math.sqrt(8000**2-6000**2)},
        "verification":["vertical/pitch/roll force and moment balances","zero-load and linear-p=1 cases",
            "concavity penalty at fixed total load","symmetric-roll-share case","independent secant root checks",
            "residual sign on both sides of each usable root","monotonicity with fixed-distribution downforce",
            "wheel-load domain and invalid-input rejection","XML parse and finite JSON"]}


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parents[1]/"content"/"assets"/"roadster-contact")
    args=parser.parse_args()
    check_model()
    check_extended_model()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    for lang in ["en","zh"]:
        for mobile in [False,True]:
            for renderer in [render_loads,render_sensitivity,render_balance,render_designs,render_combined,render_pressure]:
                renderer(args.output_dir,lang,mobile)
    render_cover(args.output_dir)
    result=summary()
    (args.output_dir/"model.json").write_text(json.dumps(result,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
    print(json.dumps({"checks":"passed","output":str(args.output_dir),"best_beta":result["cornering_force_bounds"]["best_beta"],"no_downforce":result["cornering_force_bounds"]["no_downforce"]},indent=2))

if __name__ == "__main__":
    main()
