#!/usr/bin/env python3
"""Render the Roadster article's two illustrative force-budget figures.

Requires Python 3, NumPy and Matplotlib. From any working directory:
    python3 scripts/render-roadster-models.py [--output-dir DIRECTORY]

The default output directory is content/assets/roadster relative to this file's
repository. SVGs have a 360 by 360 viewBox with 24-point labels. PNG previews are
written to the system temporary directory, never to the asset directory.

Both models use g = 9.80665 m/s² and a constant friction coefficient mu = 1.
The launch model assumes fixed mass, zero resistance, zero lateral/vertical
thrust, sufficient wheel drive to use all grip, and constant suction and forward
thrust. The target is 0 to 60 mph (26.8224 m/s) in 1.000 s without rollout.
Writing d = D/(mg) and s = Tx/(mg), the target requires d + s >= target/g - 1.
Point A has (d, s) = (1, 0.5); B adds thrust and C adds suction to reach the line.

The fan model is one fixed operating point, not a simulation of a launch:
    mass = 1,000 kg; speed = 20 m/s; shared DC power pool = 500 kW;
    effective suction area = 4 m²; reference pressure difference = 2,500 Pa;
    reference volume flow = 12 m³/s; fan aerodynamic efficiency = 0.6;
    fan motor efficiency = 0.9; wheel-drive efficiency = 0.9;
    wheel torque-force ceiling = 25 kN; wheel mechanical-power ceiling = 450 kW;
    resistance = 0; direct thrust = 0; downforce range = 0 to 20 kN.
For dp = D/area, Q = 12 sqrt(dp/2500) and Pfan,e = dp Q/(0.6*0.9).
Fdrive = min(25 kN, 450 kW/v, 0.9*(500 kW - Pfan,e)/v).
The plotted acceleration is min(mu*(mg+D), Fdrive)/m. The chosen range leaves
positive DC power for the drivetrain; the model is not extrapolated beyond it.
The peak is found by bisection at the grip/drive intersection, without SciPy.

These are assumed values for an explanatory model, not Tesla specifications or
measurements. The article supplies the bilingual captions and legends.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import tempfile


# Keep Matplotlib's disposable font/configuration cache outside the repository.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "roadster-mpl"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


G = 9.80665
TARGET_SPEED = 26.8224
TARGET_TIME = 1.0
TARGET_G = TARGET_SPEED / (G * TARGET_TIME)
LAUNCH_THRESHOLD = TARGET_G - 1.0

MASS = 1_000.0
SPEED = 20.0
MU = 1.0
DC_POWER = 500_000.0
EFFECTIVE_AREA = 4.0
REFERENCE_PRESSURE = 2_500.0
REFERENCE_FLOW = 12.0
FAN_EFFICIENCY = 0.6
FAN_MOTOR_EFFICIENCY = 0.9
DRIVE_EFFICIENCY = 0.9
TORQUE_FORCE_LIMIT = 25_000.0
WHEEL_POWER_LIMIT = 450_000.0
MAX_DOWNFORCE = 20_000.0

PAPER = "#111a20"
INK = "#f2eee3"
MUTED = "#b6c1c0"
GRID = "#60737a"
MINT = "#bcecc0"
BLUE = "#97d5f5"
AMBER = "#f4bb83"


def fan_operating_point(downforce: float) -> tuple[float, float, float, float]:
    """Return fan DC power, grip force, drive force and acceleration in SI units."""
    if not 0.0 <= downforce <= MAX_DOWNFORCE:
        raise ValueError("Downforce must be within the illustrated 0–20 kN range")
    pressure = downforce / EFFECTIVE_AREA
    flow = REFERENCE_FLOW * math.sqrt(pressure / REFERENCE_PRESSURE)
    fan_power = pressure * flow / (FAN_EFFICIENCY * FAN_MOTOR_EFFICIENCY)
    available_dc = DC_POWER - fan_power
    assert available_dc >= 0.0
    grip_force = MU * (MASS * G + downforce)
    drive_force = min(
        TORQUE_FORCE_LIMIT,
        WHEEL_POWER_LIMIT / SPEED,
        DRIVE_EFFICIENCY * available_dc / SPEED,
    )
    return fan_power, grip_force, drive_force, min(grip_force, drive_force) / MASS


def grip_drive_intersection() -> float:
    """Bisect the unique crossing of rising grip and falling drive capacity."""
    def difference(downforce: float) -> float:
        _, grip, drive, _ = fan_operating_point(downforce)
        return grip - drive

    low, high = 0.0, MAX_DOWNFORCE
    assert difference(low) < 0.0 < difference(high)
    for _ in range(80):
        middle = (low + high) / 2.0
        if difference(middle) < 0.0:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def verify_models() -> float:
    """Check the independently calculated values used in the article and figures."""
    assert math.isclose(TARGET_G, 2.735123615097918, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(TARGET_SPEED / (2.5 * G), 1.0940494460391674, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(1.0 + (LAUNCH_THRESHOLD - 1.0), LAUNCH_THRESHOLD, abs_tol=1e-12)
    assert math.isclose((LAUNCH_THRESHOLD - 0.5) + 0.5, LAUNCH_THRESHOLD, abs_tol=1e-12)

    crossing = grip_drive_intersection()
    power_a, grip_a, drive_a, acceleration_a = fan_operating_point(0.0)
    power_b, grip_b, drive_b, acceleration_b = fan_operating_point(crossing)
    power_c, grip_c, drive_c, acceleration_c = fan_operating_point(MAX_DOWNFORCE)
    assert power_a == 0.0
    assert math.isclose(drive_a, 22_500.0, abs_tol=1e-9)
    assert math.isclose(acceleration_a / G, 1.0, abs_tol=1e-12)
    assert math.isclose(crossing, 10_140.483933, rel_tol=0.0, abs_tol=1e-6)
    assert math.isclose(grip_b, drive_b, rel_tol=0.0, abs_tol=1e-8)
    assert math.isclose(acceleration_b / G, 2.034041587, rel_tol=0.0, abs_tol=1e-8)
    assert math.isclose(acceleration_c / G, 1.573313230, rel_tol=0.0, abs_tol=1e-8)
    assert power_a < power_b < power_c < DC_POWER
    assert grip_a < grip_b < grip_c
    assert drive_a > drive_b > drive_c
    assert acceleration_b > max(acceleration_a, acceleration_c)
    # Values on either side must be lower: B is a maximum of the lower envelope.
    assert fan_operating_point(crossing - 100.0)[3] < acceleration_b
    assert fan_operating_point(crossing + 100.0)[3] < acceleration_b
    return crossing


def make_axes():
    fig = plt.figure(figsize=(5, 5), facecolor=PAPER)
    ax = fig.add_axes([0.17, 0.18, 0.76, 0.71], facecolor=PAPER)
    ax.set_axisbelow(True)
    ax.grid(color=GRID, alpha=0.25, lw=0.6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=24, length=0, pad=7)
    return fig, ax


def save_figure(fig, output_dir: Path, filename: str, preview_name: str, model: str) -> None:
    svg_path = output_dir / filename
    fig.savefig(
        svg_path,
        metadata={"Date": None, "Creator": f"FENG · Matplotlib · {model}"},
    )
    # Matplotlib adds trailing spaces inside multiline SVG paths.
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8")
    preview = Path(tempfile.gettempdir()) / preview_name
    fig.savefig(preview, dpi=150)
    plt.close(fig)
    print(f"Saved {output_dir / filename}; preview {preview}")


def render_launch(output_dir: Path) -> None:
    # Retain the initial launch figure's geometry, marker shapes and stable IDs.
    plt.rcParams["svg.hashsalt"] = "roadster-launch-v1"
    fig, ax = make_axes()
    x = np.linspace(0, 2, 401)
    lower_boundary = np.maximum(0, LAUNCH_THRESHOLD - x)
    ax.fill_between(x, lower_boundary, 2, color=MINT, alpha=0.12, linewidth=0)
    ax.plot([0, LAUNCH_THRESHOLD], [LAUNCH_THRESHOLD, 0], color=MINT, lw=2.4)
    ax.set(xlim=(0, 2), ylim=(0, 2), xticks=[0, 1, 2], yticks=[0, 1, 2])
    ax.text(0, 1.08, r"$T_x/(mg)$", transform=ax.transAxes, color=INK, fontsize=24)
    ax.set_xlabel(r"$D/(mg)$", color=INK, labelpad=10, fontsize=24)
    ax.plot([1, 1], [0.5, LAUNCH_THRESHOLD - 1], color=AMBER, lw=1.2, ls=(0, (3, 3)))
    ax.plot([1, LAUNCH_THRESHOLD - 0.5], [0.5, 0.5], color=BLUE, lw=1.2, ls=(0, (3, 3)))
    for label, px, py, color, marker, offset in [
        ("A", 1, 0.5, INK, "o", (-27, -25)),
        ("B", 1, LAUNCH_THRESHOLD - 1, AMBER, "D", (-34, 0)),
        ("C", LAUNCH_THRESHOLD - 0.5, 0.5, BLUE, "s", (14, 0)),
    ]:
        ax.plot(px, py, marker=marker, ms=8, mec=PAPER, mew=1.2, color=color, zorder=5)
        ax.annotate(label, (px, py), xytext=offset, textcoords="offset points", color=color, fontsize=24)
    save_figure(fig, output_dir, "launch-design-space.svg", "roadster-launch-map.png", "illustrative constant-force model")


def render_fan(output_dir: Path, crossing: float) -> None:
    plt.rcParams["svg.hashsalt"] = "roadster-fan-v1"
    fig, ax = make_axes()
    # Include the exact crossing so the solid envelope preserves its sharp peak.
    downforce = np.sort(np.append(np.linspace(0, MAX_DOWNFORCE, 401), crossing))
    points = np.array([fan_operating_point(float(d)) for d in downforce])
    grip_g = points[:, 1] / (MASS * G)
    drive_g = points[:, 2] / (MASS * G)
    acceleration_g = points[:, 3] / G
    ax.plot(downforce / 1_000, grip_g, color=BLUE, lw=1.7, ls=(0, (5, 3)), zorder=2)
    ax.plot(downforce / 1_000, drive_g, color=AMBER, lw=1.7, ls=(0, (5, 3)), zorder=2)
    ax.plot(downforce / 1_000, acceleration_g, color=MINT, lw=3.2, zorder=3)
    ax.set(xlim=(0, 20), ylim=(0, 3.3), xticks=[0, 10, 20], yticks=[0, 1, 2, 3])
    ax.text(0, 1.08, r"$a/g$", transform=ax.transAxes, color=INK, fontsize=24)
    ax.set_xlabel(r"$D/\mathrm{kN}$", color=INK, labelpad=10, fontsize=24)
    for label, d, color, marker, offset in [
        ("A", 0.0, INK, "o", (12, -26)),
        ("B", crossing, MINT, "D", (-8, 15)),
        ("C", MAX_DOWNFORCE, AMBER, "s", (-31, -29)),
    ]:
        value = fan_operating_point(d)[3] / G
        ax.plot(d / 1_000, value, marker=marker, ms=8, mec=PAPER, mew=1.2, color=color, zorder=5, clip_on=False)
        ax.annotate(label, (d / 1_000, value), xytext=offset, textcoords="offset points", color=color, fontsize=24)
    save_figure(fig, output_dir, "fan-power-tradeoff.svg", "roadster-fan-power-tradeoff.png", "illustrative shared-power model")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "content" / "assets" / "roadster",
        help="SVG destination; defaults to the repository's Roadster asset directory",
    )
    args = parser.parse_args()
    crossing = verify_models()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 24,
        "svg.fonttype": "path",
        "mathtext.fontset": "dejavusans",
    })
    render_launch(args.output_dir)
    render_fan(args.output_dir, crossing)
    print(f"Launch: target={TARGET_G:.9f}g; d+s≥{LAUNCH_THRESHOLD:.9f}; A minimum time={TARGET_SPEED / (2.5 * G):.9f}s")
    for name, downforce in (("A", 0.0), ("B", crossing), ("C", MAX_DOWNFORCE)):
        fan_power, grip, drive, acceleration = fan_operating_point(downforce)
        print(f"Fan {name}: D={downforce / 1000:.9f}kN; Pfan,e={fan_power / 1000:.9f}kW; grip={grip / 1000:.9f}kN; drive={drive / 1000:.9f}kN; a={acceleration / G:.9f}g")


if __name__ == "__main__":
    main()
