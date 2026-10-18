#!/usr/bin/env python3
"""Render the Roadster series' illustrative torque, power and gearing figures.

Requires Python 3, NumPy and Matplotlib; Chinese labels use Noto Sans CJK.
Run from any working directory:
    python3 scripts/render-roadster-power.py [--output-dir DIRECTORY]

These curves describe an assumed equivalent drivetrain, NOT a Tesla Roadster
specification or measurement. The equivalent motor aggregates all drive motors.
Its torque limit is 1,200 N m, shaft power limit 931 kW, speed limit 18,000 rpm.
The 931 kW comes from a 1 MW pack-terminal DC budget, 98% inverter efficiency
and 95% motor efficiency. A 97% efficient fixed reduction leaves 903.07 kW at
the wheels. G = motor speed / wheel speed is either 9 or 12; wheel radius is
0.34 m. Mass is 2,000 kg, g = 9.81 m/s², and constant tyre coefficient mu = 1.5.
The baseline has no aerodynamic or rolling resistance, rotating inertia,
additional downforce, battery sag or temperature dependence. The power budget is a ceiling,
not a claim that the pack draws 1 MW at zero speed or at every operating point.

T(omega) = min(Tmax, Pshaft/omega); Fdrive(v) = eta_gear G T(Gv/r)/r.
The delivered force is min(Fdrive, mu m g). Curves end at the motor-speed limit;
their endpoints are kinematic cutoffs, NOT predicted road top speeds.

Two independent calculations accompany the plots in model.json. The first
back-solves a 2,000 kg car accelerating uniformly from rest to 26.8224 m/s in
1 s, without resistance or rotating inertia. It does not apply the plotted
drivetrain's torque, power or grip limits: it calculates what would be needed.
The second rewires 4,000 identical 3.6 V, 5 Ah, 2 milliohm cells as 100s40p or
200s20p. Currents are referenced to nominal pack voltage at 360 kW; voltage sag
is ignored for that comparison. Cell heat is counted across all 4,000 cells;
external conductor loss uses an unchanged 2 milliohm external resistance.

Assets: motor-envelope-{en,zh}.svg, gearing-envelope-{en,zh}.svg, cover.png,
and model.json. target-envelope, gear-sensitivity, downforce-budget and
one-second-feasibility each have {en,zh}.svg (864×540 desktop) and
mobile-{en,zh}.svg (432×720 mobile) editions. PNG previews go to the system
temporary directory. The decision maps test gear sensitivity, the mass/power
cost of prescribed downforce and required capacities for a one-second launch.
The model includes controlled power, gear, mass and constant-downforce
comparisons; its piecewise analytic launch times are checked against a separate
numerical integral derived from motor speed and torque. Other checks cover
power conservation, transition continuity, limits and monotonicity.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "roadster-power-mpl"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont
import numpy as np


MASS = 2_000.0
GRAVITY = 9.81
MU = 1.5
RADIUS = 0.34
TORQUE_MAX = 1_200.0
RPM_MAX = 18_000.0
DC_POWER = 1_000_000.0
ETA_INVERTER = 0.98
ETA_MOTOR = 0.95
ETA_GEAR = 0.97
SHAFT_POWER = DC_POWER * ETA_INVERTER * ETA_MOTOR
WHEEL_POWER = SHAFT_POWER * ETA_GEAR
OMEGA_BASE = SHAFT_POWER / TORQUE_MAX
RPM_BASE = OMEGA_BASE * 60 / (2 * math.pi)
OMEGA_MAX = RPM_MAX * 2 * math.pi / 60
GRIP = MU * MASS * GRAVITY
GEARS = (9, 12)

PAPER = "#111a20"
INK = "#f2eee3"
MUTED = "#b6c1c0"
GRID = "#60737a"
MINT = "#bcecc0"
BLUE = "#97d5f5"
AMBER = "#f4bb83"


def motor_torque(omega):
    """Equivalent motor shaft torque; no extrapolation past rated maximum speed."""
    omega = np.asarray(omega, dtype=float)
    if np.any(omega < 0) or np.any(omega > OMEGA_MAX * (1 + 1e-12)):
        raise ValueError("Motor speed outside the illustrated operating range")
    return np.minimum(TORQUE_MAX, SHAFT_POWER / np.maximum(omega, 1e-12))


def speed_limit(gear):
    return OMEGA_MAX * RADIUS / gear


def drive_force(speed, gear):
    """Available longitudinal wheel force before applying the tyre grip limit."""
    speed = np.asarray(speed, dtype=float)
    return ETA_GEAR * gear * motor_torque(gear * speed / RADIUS) / RADIUS


def road_force(speed, gear):
    return np.minimum(drive_force(speed, gear), GRIP)


def launch_requirement():
    """Back-solve the target; independent of the plotted envelope's capacity."""
    speed, duration, gear = 26.8224, 1.0, 9
    acceleration = speed / duration
    force = MASS * acceleration
    energy = 0.5 * MASS * speed**2
    end_wheel_power = force * speed
    return {
        "kind": "required-performance-not-a-prediction-of-the-plotted-drivetrain",
        "assumptions": {"mass_kg": MASS, "start_speed_m_s": 0, "end_speed_m_s": speed,
            "duration_s": duration, "constant_acceleration": True, "g_m_s2": GRAVITY,
            "road_resistance_N": 0, "rotational_inertia_kg_m2": 0,
            "downforce_N": 0, "direct_thrust_N": 0, "radius_m": RADIUS, "gear_ratio": gear,
            "inverter_efficiency": ETA_INVERTER, "motor_efficiency": ETA_MOTOR, "gear_efficiency": ETA_GEAR},
        "acceleration_m_s2": acceleration, "distance_m": speed * duration / 2,
        "kinetic_energy_J": energy, "kinetic_energy_kWh": energy / 3_600_000,
        "mean_wheel_power_W": energy / duration, "end_wheel_power_W": end_wheel_power,
        "required_force_N": force, "required_mu_without_downforce": force / (MASS * GRAVITY),
        "required_aggregate_motor_torque_Nm": force * RADIUS / (gear * ETA_GEAR),
        "end_dc_power_W": end_wheel_power / (ETA_INVERTER * ETA_MOTOR * ETA_GEAR),
    }


def nominal_pack_comparison():
    """Compare wiring using nominal voltage, not a voltage-sag operating solution."""
    count, voltage, capacity, resistance = 4_000, 3.6, 5.0, 0.002
    power, external_resistance = 360_000.0, 0.002
    packs = []
    for series, parallel in [(100, 40), (200, 20)]:
        pack_voltage = series * voltage
        pack_current = power / pack_voltage
        cell_current = pack_current / parallel
        pack_resistance = series * resistance / parallel
        packs.append({
            "series_cells": series, "parallel_cells": parallel,
            "nominal_voltage_V": pack_voltage, "capacity_Ah": parallel * capacity,
            "nominal_energy_kWh": pack_voltage * parallel * capacity / 1000,
            "pack_current_A": pack_current, "cell_current_A": cell_current,
            "cell_C_rate": cell_current / capacity,
            "equivalent_cell_resistance_ohm": pack_resistance,
            "all_cells_loss_W": count * cell_current**2 * resistance,
            "fixed_external_conductor_loss_W": pack_current**2 * external_resistance,
        })
    return {
        "kind": "nominal-voltage-wiring-comparison-not-a-loaded-pack-solution",
        "assumptions": {"cell_count": count, "cell_nominal_voltage_V": voltage,
            "cell_capacity_Ah": capacity, "cell_resistance_ohm": resistance,
            "nominal_voltage_referenced_power_W": power,
            "voltage_sag_ignored": True, "external_conductor_resistance_ohm": external_resistance,
            "identical_cells": True},
        "packs": packs,
    }



def scenario(key, *, mass=MASS, mu=MU, gear=9, power_scale=1.0, downforce=0.0,
             group="single-variable", changed="none"):
    """Quasi-static full-envelope launch, with an explicit valid speed domain.

    D is a constant normal force available even at rest. It is NOT a passive
    aerodynamic-downforce model. C = min(mu (mg+D), eta_g G Tmax/r) is the
    low-speed force cap. Full-force time follows from dt = m dv/F(v).
    """
    values = [mass, mu, gear, power_scale, TORQUE_MAX, RADIUS, ETA_GEAR]
    if any(value <= 0 for value in values) or downforce < 0:
        raise ValueError("Scenario parameters must be positive; downforce cannot be negative")
    target_speed = launch_requirement()["assumptions"]["end_speed_m_s"]
    max_speed = speed_limit(gear)
    if target_speed > max_speed:
        raise ValueError("The requested launch exceeds the motor-speed cutoff")
    grip = mu * (mass * GRAVITY + downforce)
    torque_force = ETA_GEAR * gear * TORQUE_MAX / RADIUS
    force_cap = min(grip, torque_force)
    power = WHEEL_POWER * power_scale
    transition = power / force_cap
    first = min(target_speed, transition)
    constant_force_time = mass * first / force_cap
    constant_power_time = mass * max(0.0, target_speed**2 - transition**2) / (2 * power)
    initial_limit = "tyre-grip" if grip < torque_force else "motor-torque"
    final_limit = initial_limit if target_speed <= transition else "wheel-power"
    end_force = min(force_cap, power / target_speed)
    return {
        "id": key, "group": group, "changed": changed,
        "assumptions": {"mass_kg": mass, "mu": mu, "gear_ratio": gear,
            "pack_dc_W": DC_POWER * power_scale, "wheel_power_W": power,
            "aggregate_motor_torque_Nm": TORQUE_MAX, "constant_downforce_N": downforce,
            "motor_rpm_limit": RPM_MAX, "radius_m": RADIUS},
        "grip_limit_N": grip, "torque_force_limit_N": torque_force,
        "low_speed_force_cap_N": force_cap,
        "full_envelope_0_to_60_mph_time_s": constant_force_time + constant_power_time,
        "constant_force_phase_s": constant_force_time,
        "constant_power_phase_s": constant_power_time,
        "acceleration_at_0_m_s2": force_cap / mass,
        "acceleration_at_60_mph_m_s2": end_force / mass,
        "force_at_0_N": force_cap, "force_at_60_mph_N": end_force,
        "initial_limiting_constraint": initial_limit, "final_limiting_constraint": final_limit,
        "force_to_power_transition_kmh": transition * 3.6,
        "motor_torque_to_power_transition_kmh": power / torque_force * 3.6,
        "transition_within_motor_speed_range": transition <= max_speed,
        "motor_speed_cutoff_kmh": max_speed * 3.6,
    }


def design_comparison():
    constant_downforce = 0.25 * MASS * GRAVITY
    return {
        "kind": "controlled-illustrative-force-envelope-comparison-not-vehicle-performance",
        "target_speed_m_s": 26.8224, "target_speed_kmh": 96.56064,
        "time_definition": "minimum theoretical rest-to-60-mph time from dt=m*dv/F(v), no rollout",
        "time_formula": "C=min(mu*(m*g+D),eta_g*G*Tmax/r); vb=Pwheel/C; t=m*v/C for v<=vb, otherwise m*vb/C+m*(v^2-vb^2)/(2*Pwheel)",
        "force_domain": "At v=0 use C; at 0<v<=omega_max*r/G use min(C,Pwheel/v). No extrapolation above maximum motor speed.",
        "constraints": [
            "Each single-variable case changes only the named parameter; baseline mass is 2000 kg.",
            "The power increase assumes the battery, inverter and motor all support 25 percent more power, with unchanged torque and rpm limits and efficiencies.",
            "The mass reduction retains all powertrain limits and tyre mu; it is not a claim that a lighter physical battery has identical power.",
            "Downforce is a constant 4905 N referenced to the baseline mass, available at rest; not speed-squared passive aero downforce.",
            "Downforce has zero added mass, zero power draw and zero drag here; its benefit is an optimistic isolated-force bound, not a realizable system estimate.",
            "Combined modifications have their own group and are not a single-variable comparison.",
            "The separate mu=2 cases test when the original gear-ratio conclusion reverses; they do not assert a real tyre coefficient.",
        ],
        "omissions": ["aerodynamic drag", "rolling resistance", "rotational inertia", "tyre load sensitivity",
            "axle load transfer", "traction-control transients", "battery sag", "thermal limits",
            "power needed to create downforce", "mass or drag added by downforce hardware"],
        "scenarios": [
            scenario("baseline", group="baseline"),
            scenario("power-plus-25-percent", power_scale=1.25, changed="DC and wheel power +25 percent; torque unchanged"),
            scenario("gear-12", gear=12, changed="gear ratio 9 to 12"),
            scenario("mass-minus-10-percent", mass=0.9 * MASS, changed="mass 2000 to 1800 kg"),
            scenario("constant-downforce", downforce=constant_downforce, changed="constant normal downforce +4905 N"),
            scenario("downforce-and-gear-12", gear=12, downforce=constant_downforce,
                group="combined", changed="constant normal downforce +4905 N and gear ratio 9 to 12"),
            scenario("mu-2-gear-9", mu=2.0, group="sensitivity", changed="mu 1.5 to 2.0"),
            scenario("mu-2-gear-12", mu=2.0, gear=12, group="sensitivity", changed="mu 1.5 to 2.0 and gear ratio 9 to 12"),
        ],
    }


def envelope_time(mass, force_cap, wheel_power, target_speed=26.8224):
    """Exact rest-to-speed integral for F(v)=min(C,P/v), including v=0.

    This capacity model has no gear-speed extrapolation. Callers that map it to
    a motor must separately verify target_speed <= omega_max*r/G.
    Arrays are accepted for sensitivity maps; every capacity must be finite and
    strictly positive, while target speed may be zero.
    """
    mass, force_cap, wheel_power = np.broadcast_arrays(mass, force_cap, wheel_power)
    if any(np.any(~np.isfinite(value)) or np.any(value <= 0)
           for value in [mass, force_cap, wheel_power]):
        raise ValueError("Mass, force and power must be finite and positive")
    if not math.isfinite(target_speed) or target_speed < 0:
        raise ValueError("Target speed must be finite and nonnegative")
    boundary = wheel_power / force_cap
    return np.where(target_speed <= boundary, mass * target_speed / force_cap,
        mass * (target_speed**2 / wheel_power + wheel_power / force_cap**2) / 2)


def required_wheel_power(mass, force_cap, duration, target_speed=26.8224):
    """Minimum constant wheel-power ceiling meeting the time target.

    None denotes an impossible force cap. The stable smaller quadratic root
    is the only root on the constant-power branch (P <= C*vf).
    """
    if any(not math.isfinite(v) or v <= 0 for v in [mass, force_cap, duration, target_speed]):
        raise ValueError("Inversion inputs must be finite and positive")
    force_only_time = mass * target_speed / force_cap
    if force_only_time > duration * (1 + 1e-12):
        return None
    discriminant = 0.0 if math.isclose(force_only_time, duration, rel_tol=1e-12) else max(0.0, duration**2 - force_only_time**2)
    return mass * target_speed**2 / (duration + math.sqrt(discriminant))


def gear_sensitivity_sample(mu, pack_dc):
    if not math.isfinite(mu) or mu <= 0 or not math.isfinite(pack_dc) or pack_dc <= 0:
        raise ValueError("Coefficient and DC power must be finite and positive")
    force_caps = [min(mu * MASS * GRAVITY, ETA_GEAR * gear * TORQUE_MAX / RADIUS) for gear in GEARS]
    times = [float(envelope_time(MASS, cap, pack_dc * ETA_INVERTER * ETA_MOTOR * ETA_GEAR)) for cap in force_caps]
    return {"mu": mu, "pack_dc_W": pack_dc, "gear_9_time_s": times[0],
        "gear_12_time_s": times[1], "time_reduction_percent": 100 * (1 - times[1] / times[0])}


def downforce_budget_sample(added_mass, auxiliary_dc):
    if not math.isfinite(added_mass) or added_mass < 0 or not math.isfinite(auxiliary_dc) or not 0 <= auxiliary_dc < DC_POWER:
        raise ValueError("Added mass must be nonnegative and auxiliary DC draw in [0, 1 MW)")
    mass, downforce = MASS + added_mass, 0.25 * MASS * GRAVITY
    cap = min(MU * (mass * GRAVITY + downforce), ETA_GEAR * 12 * TORQUE_MAX / RADIUS)
    wheel_power = (DC_POWER - auxiliary_dc) * ETA_INVERTER * ETA_MOTOR * ETA_GEAR
    duration = float(envelope_time(mass, cap, wheel_power))
    baseline = float(envelope_time(MASS, GRIP, WHEEL_POWER))
    return {"added_mass_kg": added_mass, "auxiliary_dc_W": auxiliary_dc,
        "total_mass_kg": mass, "low_speed_force_cap_N": cap,
        "drive_dc_W": DC_POWER - auxiliary_dc, "wheel_power_W": wheel_power,
        "time_s": duration, "time_reduction_percent": 100 * (1 - duration / baseline)}


def auxiliary_break_even(added_mass):
    """Return available DC allowance; a negative value means even Paux=0 loses."""
    mass = MASS + added_mass
    cap = min(MU * (mass * GRAVITY + 0.25 * MASS * GRAVITY), ETA_GEAR * 12 * TORQUE_MAX / RADIUS)
    baseline = float(envelope_time(MASS, GRIP, WHEEL_POWER))
    required = required_wheel_power(mass, cap, baseline)
    return None if required is None else DC_POWER - required / (ETA_INVERTER * ETA_MOTOR * ETA_GEAR)


def decision_boundaries():
    baseline = float(envelope_time(MASS, GRIP, WHEEL_POWER))
    # At zero auxiliary draw the actual power budget becomes restrictive before
    # the force-only maximum-mass bound. Solve the monotone time equation.
    lo, hi = 0.0, 900.0
    for _ in range(70):
        middle = (lo + hi) / 2
        if downforce_budget_sample(middle, 0)["time_s"] > baseline:
            hi = middle
        else:
            lo = middle
    constant = launch_requirement()
    varying_cap = 70_000.0
    varying_power = required_wheel_power(MASS, varying_cap, 1.0)
    return {
        "kind": "conditional-capacity-boundaries-not-real-vehicle-predictions",
        "shared_model": {"target_speed_m_s": 26.8224, "mass_kg": MASS,
            "time_formula": "t=m*vf/C if P>=C*vf; otherwise t=m/2*(vf^2/P+P/C^2)",
            "power_inversion": "Pmin=m*vf^2/(T+sqrt(T^2-(m*vf/C)^2)); feasible only C>=m*vf/T",
            "force_definition": "C=min(mu*(m*g+D),eta_g*G*Tmax/r)",
            "force_at_rest": "C; P/v is not evaluated at v=0",
            "omissions": ["drag", "rolling resistance", "rotational inertia", "tyre load sensitivity",
                "axle load transfer", "traction-control delay", "battery sag", "thermal derating"],
            "limits": "G=9 and G=12 examples obey the 18000 rpm cutoff through 60 mph. Capacity feasibility map is a requirements map; it does not assert those capacities are available from the assumed motor or tyres."},
        "gear_sensitivity": {
            "fixed": {"mass_kg": MASS, "aggregate_motor_torque_Nm": TORQUE_MAX, "radius_m": RADIUS,
                "downforce_N": 0, "efficiency": ETA_INVERTER * ETA_MOTOR * ETA_GEAR},
            "map_domain": {"mu": [1.0, 2.5], "available_drive_dc_MW": [0.3, 1.5]},
            "mu_where_gear_9_becomes_torque_limited": ETA_GEAR * 9 * TORQUE_MAX / RADIUS / (MASS * GRAVITY),
            "baseline_G9_torque_threshold_Nm": GRIP * RADIUS / (ETA_GEAR * 9),
            "baseline_G9_dc_power_threshold_W": GRIP * 26.8224 / (ETA_INVERTER * ETA_MOTOR * ETA_GEAR),
            "interpretation": "Within this positive-power shared-speed domain G12 has zero time benefit when mu<=the G9 torque/grip crossing. Above it G12 is faster, with the percentage benefit dependent on power. This says nothing about higher-speed gearing or tyre feasibility.",
            "samples": [gear_sensitivity_sample(mu, p) for mu, p in [(1.5, 1e6), (1.6, 1e6), (2.0, 1e6), (2.0, 0.4e6), (2.5, 1.5e6)]],
        },
        "downforce_budget": {
            "reference_time_s": baseline,
            "fixed": {"mu": MU, "gear": 12, "constant_downforce_N": 4905.0, "shared_pack_dc_W": DC_POWER,
                "wheel_power_formula": "(1000000-Paux)*0.98*0.95*0.97",
                "force_formula": "min(1.5*((2000+dm)*9.81+4905),0.97*12*1200/0.34)"},
            "map_domain": {"added_mass_kg": [0, 900], "auxiliary_dc_kW": [0, 600]},
            "interpretation": "A budget sensitivity only: D is prescribed independently of auxiliary draw and hardware mass. The map does not establish that any fan can generate D at the specified cost. Positive gain is relative to original G9 without downforce, not to G12 without hardware.",
            "constant_draw_at_rest": "Auxiliary power is drawn from rest; unused drive power is not automatically transferred to the auxiliary system. No time-dependent sharing or fan ramp-up is modelled.",
            "zero_auxiliary_break_even_added_mass_kg": (lo + hi) / 2,
            "break_even_samples": [{"added_mass_kg": dm, "auxiliary_dc_W": auxiliary_break_even(dm)} for dm in [0, 200, 500, 700]],
            "samples": [downforce_budget_sample(dm, aux) for dm, aux in [(0, 0), (100, 100_000), (200, 200_000), (500, 300_000)]],
        },
        "one_second_feasibility": {
            "mass_kg": MASS, "duration_s": 1, "target_speed_m_s": 26.8224,
            "map_domain": {"effective_force_cap_kN": [25, 100], "wheel_power_MW": [0.4, 1.8]},
            "minimum_effective_force_cap_N": MASS * 26.8224,
            "infinite_force_power_infimum_W": MASS * 26.8224**2 / 2,
            "existing_G12_torque_only_time_lower_bound_s": MASS * 26.8224 / (ETA_GEAR * 12 * TORQUE_MAX / RADIUS),
            "interpretation": "All points on/above the boundary and at/right of the minimum-force boundary meet the one-second ideal time. C requires both tyre force and torque-derived wheel force at least C; it is not a command or average force. Baseline 1200 Nm, mu=1.5 and 1 MW DC are not held fixed when designing capacities.",
            "uniform_acceleration_corner": {"force_cap_N": constant["required_force_N"],
                "wheel_power_W": constant["end_wheel_power_W"], "time_s": 1},
            "variable_acceleration_example": {"force_cap_N": varying_cap, "wheel_power_W": varying_power,
                "dc_power_at_assumed_efficiency_W": varying_power / (ETA_INVERTER * ETA_MOTOR * ETA_GEAR),
                "transition_speed_m_s": varying_power / varying_cap,
                "initial_acceleration_m_s2": varying_cap / MASS,
                "final_acceleration_m_s2": varying_power / (MASS * 26.8224),
                "time_s": float(envelope_time(MASS, varying_cap, varying_power)),
                "minimum_mu_without_downforce": varying_cap / (MASS * GRAVITY),
                "required_aggregate_torque_G9_Nm": varying_cap * RADIUS / (ETA_GEAR * 9)},
            "boundary_samples": [{"force_cap_N": c, "minimum_wheel_power_W": required_wheel_power(MASS, c, 1)}
                for c in [50_000, 53_644.8, 60_000, 70_000, 80_000, 100_000]],
        },
    }


def check_decision_boundaries():
    data = decision_boundaries()
    target = 26.8224
    # Midpoint quadrature independently samples F at positive speeds. Tests span
    # both branches, transition exactly at target, impossible target and zero v.
    for mass, cap, power in [(MASS, GRIP, WHEEL_POWER), (MASS, 70_000, 876_075.4765948907),
                            (MASS, 53_644.8, 1_438_882.28352), (2200, 38_259, 722_456),
                            (2500, 41_082.35294117647, 632_149)]:
        n = 400_000
        velocity = (np.arange(n) + 0.5) * target / n
        force = np.minimum(cap, power / velocity)
        numerical = float(np.mean(mass / force) * target)
        analytic = float(envelope_time(mass, cap, power))
        assert math.isclose(numerical, analytic, rel_tol=1e-10)
        recovered = required_wheel_power(mass, cap, analytic)
        assert recovered is not None and recovered <= power * (1 + 1e-10)
        assert math.isclose(float(envelope_time(mass, cap, recovered)), analytic, rel_tol=1e-10)
    assert float(envelope_time(MASS, GRIP, WHEEL_POWER, 0)) == 0
    assert required_wheel_power(MASS, 50_000, 1) is None
    assert math.isclose(required_wheel_power(MASS, 53_644.8, 1), 1_438_882.28352, rel_tol=1e-12)
    assert all(speed_limit(gear) > target for gear in GEARS)
    # Gear samples independently re-use motor-speed/torque scenarios rather than
    # deriving C again in the comparison being checked.
    for sample in data["gear_sensitivity"]["samples"]:
        for gear in GEARS:
            independent = scenario("check", gear=gear, mu=sample["mu"], power_scale=sample["pack_dc_W"] / DC_POWER)
            assert math.isclose(sample[f"gear_{gear}_time_s"], independent["full_envelope_0_to_60_mph_time_s"], rel_tol=1e-12)
        assert sample["time_reduction_percent"] >= -1e-12
    budget = data["downforce_budget"]
    n = 400_000
    velocity = (np.arange(n) + 0.5) * target / n
    for sample in budget["break_even_samples"]:
        case = downforce_budget_sample(sample["added_mass_kg"], sample["auxiliary_dc_W"])
        assert math.isclose(case["time_s"], budget["reference_time_s"], rel_tol=1e-12)
        # Independent midpoint quadrature reconstructs wheel force from motor
        # speed and torque, including the auxiliary load's DC subtraction.
        mass = MASS + sample["added_mass_kg"]
        omega = 12 * velocity / RADIUS
        shaft_power = (DC_POWER - sample["auxiliary_dc_W"]) * ETA_INVERTER * ETA_MOTOR
        torque = np.minimum(TORQUE_MAX, shaft_power / omega)
        wheel_force = ETA_GEAR * 12 * torque / RADIUS
        delivered = np.minimum(wheel_force, MU * (mass * GRAVITY + 4905))
        numerical = float(np.mean(mass / delivered) * target)
        assert math.isclose(numerical, budget["reference_time_s"], rel_tol=1e-10)
        for auxiliary_scale, should_win in [(.99, True), (1.01, False)]:
            shaft = (DC_POWER - sample["auxiliary_dc_W"] * auxiliary_scale) * ETA_INVERTER * ETA_MOTOR
            force = np.minimum(ETA_GEAR * 12 * np.minimum(TORQUE_MAX, shaft / omega) / RADIUS,
                MU * (mass * GRAVITY + 4905))
            duration = float(np.mean(mass / force) * target)
            assert (duration < budget["reference_time_s"]) == should_win
    assert math.isclose(downforce_budget_sample(budget["zero_auxiliary_break_even_added_mass_kg"], 0)["time_s"], budget["reference_time_s"], rel_tol=1e-12)
    assert budget["samples"][1]["time_s"] < budget["reference_time_s"] < budget["samples"][3]["time_s"]
    for sample in data["one_second_feasibility"]["boundary_samples"]:
        c, p = sample["force_cap_N"], sample["minimum_wheel_power_W"]
        if p is None:
            assert c < MASS * target
            continue
        assert math.isclose(float(envelope_time(MASS, c, p)), 1, rel_tol=1e-12)
        numerical = float(np.mean(MASS / np.minimum(c, p / velocity)) * target)
        assert math.isclose(numerical, 1, rel_tol=1e-10)
        assert float(envelope_time(MASS, c, p * .99)) > 1
        assert float(envelope_time(MASS, c, p * 1.01)) <= 1 + 1e-12
        assert float(np.mean(MASS / np.minimum(c, p * .99 / velocity)) * target) > 1
        assert float(np.mean(MASS / np.minimum(c, p * 1.01 / velocity)) * target) <= 1 + 1e-12
    for args in [(0, GRIP, WHEEL_POWER), (MASS, float("nan"), WHEEL_POWER), (MASS, GRIP, 0)]:
        try:
            envelope_time(*args)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid capacities must be rejected")


def check_design_comparison():
    data = design_comparison()
    cases = {item["id"]: item for item in data["scenarios"]}
    target = data["target_speed_m_s"]
    for case in cases.values():
        p = case["assumptions"]
        assert all(p[key] > 0 for key in ["mass_kg", "mu", "gear_ratio", "pack_dc_W", "wheel_power_W"])
        assert 0 < target * 3.6 < case["motor_speed_cutoff_kmh"]
        assert case["acceleration_at_0_m_s2"] > 0
        assert case["acceleration_at_60_mph_m_s2"] > 0
        # Independent numerical quadrature derives force via motor torque and gear
        # speed, rather than using C or the analytical time formula above.
        velocity = np.linspace(0, target, 200_001)
        omega = velocity * p["gear_ratio"] / RADIUS
        shaft_power = p["pack_dc_W"] * ETA_INVERTER * ETA_MOTOR
        torque = np.full_like(omega, TORQUE_MAX)
        torque[1:] = np.minimum(TORQUE_MAX, shaft_power / omega[1:])
        wheel_force = torque * ETA_GEAR * p["gear_ratio"] / RADIUS
        grip = p["mu"] * (p["mass_kg"] * GRAVITY + p["constant_downforce_N"])
        delivered = np.minimum(wheel_force, grip)
        dt_dv = p["mass_kg"] / delivered
        numerical = float(np.sum((dt_dv[1:] + dt_dv[:-1]) * 0.5 * np.diff(velocity)))
        assert math.isclose(numerical, case["full_envelope_0_to_60_mph_time_s"], rel_tol=1e-10)
        assert np.all(delivered * velocity <= p["wheel_power_W"] * (1 + 1e-12))
        assert np.all(delivered <= grip * (1 + 1e-12))
        assert math.isclose(delivered[-1] / p["mass_kg"], case["acceleration_at_60_mph_m_s2"], rel_tol=1e-12)
        # Constant-force/power branches are continuous at their shared boundary.
        boundary = case["force_to_power_transition_kmh"] / 3.6
        assert math.isclose(case["low_speed_force_cap_N"], p["wheel_power_W"] / boundary, rel_tol=1e-12)
        assert case["initial_limiting_constraint"] != "wheel-power"
    base_time = cases["baseline"]["full_envelope_0_to_60_mph_time_s"]
    assert math.isclose(base_time, target / (MU * GRAVITY), rel_tol=1e-12)
    for key in ["power-plus-25-percent", "gear-12", "mass-minus-10-percent"]:
        assert math.isclose(cases[key]["full_envelope_0_to_60_mph_time_s"], base_time, rel_tol=1e-12)
    assert cases["constant-downforce"]["full_envelope_0_to_60_mph_time_s"] < base_time
    assert cases["downforce-and-gear-12"]["full_envelope_0_to_60_mph_time_s"] < cases["constant-downforce"]["full_envelope_0_to_60_mph_time_s"]
    assert cases["mu-2-gear-12"]["full_envelope_0_to_60_mph_time_s"] < cases["mu-2-gear-9"]["full_envelope_0_to_60_mph_time_s"]
    assert cases["constant-downforce"]["initial_limiting_constraint"] == "motor-torque"
    assert cases["downforce-and-gear-12"]["final_limiting_constraint"] == "wheel-power"
    # Raising a capacity cannot reduce pointwise force at shared valid speeds.
    velocity = np.linspace(0, target, 10_001)
    def forces(case):
        return np.minimum(case["low_speed_force_cap_N"], case["assumptions"]["wheel_power_W"] / np.maximum(velocity, 1e-12))
    for key in ["power-plus-25-percent", "gear-12", "constant-downforce"]:
        assert np.all(forces(cases[key]) >= forces(cases["baseline"]) - 1e-8)
    for key in ["baseline", "constant-downforce", "downforce-and-gear-12"]:
        assert cases[key]["full_envelope_0_to_60_mph_time_s"] > 1
    try:
        scenario("invalid-gear", gear=100)
    except ValueError:
        pass
    else:
        raise AssertionError("A target above the motor-speed cutoff must be rejected")


def check_model():
    check_design_comparison()
    check_decision_boundaries()
    # Independent terminal accounting: three losses plus wheel output close the budget.
    losses = [20_000.0, 49_000.0, 27_930.0]
    assert math.isclose(WHEEL_POWER + sum(losses), DC_POWER, abs_tol=1e-8)
    assert math.isclose(SHAFT_POWER, 931_000.0, abs_tol=1e-8)
    assert math.isclose(WHEEL_POWER, 903_070.0, abs_tol=1e-8)
    assert math.isclose(GRIP, 29_430.0, abs_tol=1e-8)
    assert math.isclose(RPM_BASE, 7408.662600927728, abs_tol=1e-8)
    for gear in GEARS:
        base_speed = OMEGA_BASE * RADIUS / gear
        max_speed = speed_limit(gear)
        assert base_speed < WHEEL_POWER / GRIP < max_speed
        # The independent work relation is T dtheta -> eta T dtheta, with ds=r dtheta/G.
        for speed in [1, 10, 20, 40, 50]:
            omega = gear * speed / RADIUS
            torque = float(motor_torque(omega))
            wheel_force = float(drive_force(speed, gear))
            assert math.isclose(wheel_force * speed, ETA_GEAR * torque * omega, rel_tol=1e-12)
            assert wheel_force * speed <= WHEEL_POWER * (1 + 1e-12)
            assert float(road_force(speed, gear)) <= GRIP
        left = float(drive_force(base_speed * (1 - 1e-9), gear))
        right = float(drive_force(base_speed * (1 + 1e-9), gear))
        assert math.isclose(left, right, rel_tol=2e-9)
        assert math.isclose(float(road_force(0, gear)), GRIP, abs_tol=1e-8)
        assert math.isclose(float(road_force(50, gear)), WHEEL_POWER / 50, rel_tol=1e-12)
        assert math.isclose(gear * max_speed / RADIUS, OMEGA_MAX, rel_tol=1e-12)
        try:
            drive_force(max_speed + 1, gear)
        except ValueError:
            pass
        else:
            raise AssertionError("The envelope must stop at maximum motor speed")
    assert math.isclose(speed_limit(9) / speed_limit(12), 12 / 9, rel_tol=1e-12)
    assert math.isclose(float(drive_force(0, 12) / drive_force(0, 9)), 12 / 9, rel_tol=1e-12)
    # Gear ratios redistribute torque and speed; they do not multiply peak wheel power.
    shared_speed = np.linspace(0, speed_limit(12), 500)
    np.testing.assert_allclose(road_force(shared_speed, 9), road_force(shared_speed, 12), rtol=1e-12)

    target = launch_requirement()
    assert math.isclose(target["kinetic_energy_J"], 719_441.14176, abs_tol=1e-7)
    # Work over distance and the triangular power-time integral recover kinetic energy.
    assert math.isclose(target["required_force_N"] * target["distance_m"], target["kinetic_energy_J"], rel_tol=1e-12)
    assert math.isclose(0.5 * target["end_wheel_power_W"], target["mean_wheel_power_W"], rel_tol=1e-12)
    assert math.isclose(target["required_force_N"], 53_644.8, abs_tol=1e-8)
    assert math.isclose(target["kinetic_energy_kWh"] * 3_600_000, target["kinetic_energy_J"], rel_tol=1e-12)
    assert math.isclose(target["required_aggregate_motor_torque_Nm"] * 9 * ETA_GEAR / RADIUS, target["required_force_N"], rel_tol=1e-12)
    assert math.isclose(target["end_dc_power_W"] * ETA_INVERTER * ETA_MOTOR * ETA_GEAR, target["end_wheel_power_W"], rel_tol=1e-12)
    # This target exceeds the separate plotted example in all three relevant capacities.
    assert target["required_aggregate_motor_torque_Nm"] > TORQUE_MAX
    assert target["end_dc_power_W"] > DC_POWER
    assert target["required_mu_without_downforce"] > MU

    pack_data = nominal_pack_comparison()
    for pack, expected_current, expected_external_loss in zip(pack_data["packs"], [1000, 500], [2000, 500]):
        assert pack["series_cells"] * pack["parallel_cells"] == 4_000
        assert math.isclose(pack["nominal_energy_kWh"], 72.0, abs_tol=1e-12)
        assert math.isclose(pack["pack_current_A"], expected_current, abs_tol=1e-12)
        assert math.isclose(pack["cell_current_A"], 25.0, abs_tol=1e-12)
        assert math.isclose(pack["all_cells_loss_W"], 5000.0, abs_tol=1e-12)
        # Equivalent pack I²R and the sum of individual cell losses must agree.
        assert math.isclose(pack["pack_current_A"]**2 * pack["equivalent_cell_resistance_ohm"], pack["all_cells_loss_W"], rel_tol=1e-12)
        assert math.isclose(pack["fixed_external_conductor_loss_W"], expected_external_loss, abs_tol=1e-12)


def labels(language):
    return {
        "en": {
            "torque": "Motor torque / N·m", "power": "Shaft power / kW",
            "speed": "Motor speed / 1,000 rpm", "constant_torque": "Constant\ntorque",
            "constant_power": "Constant\npower", "base": "Base speed  7,409 rpm",
            "available": "Available drive force / kN", "delivered": "After tyre limit / kN",
            "vehicle": "Vehicle speed / km/h", "grip": "Grip  29.43 kN",
            "equal": "Same launch force", "cutoffs": "Motor-speed cutoffs",
            "limit9": "G = 9   256 km/h", "limit12": "G = 12   192 km/h",
        },
        "zh": {
            "torque": "电机转矩 / N·m", "power": "轴端功率 / kW",
            "speed": "电机转速 / 千转每分", "constant_torque": "恒转矩",
            "constant_power": "恒功率", "base": "基速  7,409 rpm",
            "available": "可用驱动力 / kN", "delivered": "轮胎限制后 / kN",
            "vehicle": "车速 / km/h", "grip": "抓地上限  29.43 kN",
            "equal": "起步驱动力相同", "cutoffs": "电机转速对应的车速上限",
            "limit9": "G = 9   256 km/h", "limit12": "G = 12   192 km/h",
        },
    }[language]


def configure(language):
    font = "DejaVu Sans" if language == "en" else "Noto Sans CJK JP"
    # Fail clearly instead of silently rendering missing-glyph boxes.
    findfont(FontProperties(family=font), fallback_to_default=False)
    plt.rcParams.update({
        "font.family": font, "font.size": 20, "svg.fonttype": "path",
        "svg.hashsalt": f"roadster-power-v1-{language}", "axes.unicode_minus": False,
        "mathtext.fontset": "dejavusans", "lines.solid_capstyle": "round",
    })


def style_axis(ax, title):
    ax.set_facecolor(PAPER)
    ax.set_axisbelow(True)
    ax.grid(color=GRID, alpha=0.23, lw=0.7)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
        ax.spines[spine].set_linewidth(0.7)
    ax.tick_params(colors=MUTED, labelsize=18, length=0, pad=7)
    ax.set_title(title, loc="left", pad=18, color=INK, fontsize=20)


def save_svg(fig, output_dir, name, description):
    path = output_dir / f"{name}.svg"
    fig.savefig(path, metadata={"Date": None, "Creator": "FENG · Matplotlib", "Description": description})
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    preview = Path(tempfile.gettempdir()) / f"roadster-{name}.png"
    fig.savefig(preview, dpi=120)
    plt.close(fig)
    print(f"{path}; preview {preview}")


def render_motor(output_dir, language):
    configure(language)
    label = labels(language)
    fig = plt.figure(figsize=(6, 8.3), facecolor=PAPER)
    axes = [fig.add_axes([0.18, 0.58, 0.76, 0.30]), fig.add_axes([0.18, 0.17, 0.76, 0.29])]
    rpm = np.sort(np.r_[np.linspace(0, RPM_MAX, 600), RPM_BASE])
    omega = rpm * 2 * np.pi / 60
    torque = motor_torque(omega)
    for ax, title in zip(axes, [label["torque"], label["power"]]):
        style_axis(ax, title)
        ax.set_xlim(0, 18)
        ax.set_xticks([0, 6, 12, 18])
        ax.axvspan(0, RPM_BASE / 1000, color=MINT, alpha=0.045, lw=0)
        ax.axvspan(RPM_BASE / 1000, 18, color=BLUE, alpha=0.045, lw=0)
        ax.axvline(RPM_BASE / 1000, color=MUTED, alpha=0.55, ls=(0, (2, 3)), lw=1)
    axes[0].plot(rpm / 1000, torque, color=MINT, lw=3)
    axes[0].set(ylim=(0, 1450), yticks=[0, 600, 1200])
    axes[0].tick_params(labelbottom=False)
    axes[0].text(3.3, 540, label["constant_torque"], ha="center", va="center", color=MINT, fontsize=20, linespacing=1.4)
    axes[0].text(13.0, 1080, label["constant_power"], ha="center", va="center", color=BLUE, fontsize=20, linespacing=1.4)
    axes[0].plot(RPM_BASE / 1000, TORQUE_MAX, "o", ms=6, color=INK, mec=PAPER, mew=1)
    axes[1].plot(rpm / 1000, torque * omega / 1000, color=BLUE, lw=3)
    axes[1].set(ylim=(0, 1120), yticks=[0, 500, 1000])
    axes[1].text(12.8, 680, "931 kW", ha="center", color=BLUE, fontsize=23)
    axes[1].set_xlabel(label["speed"], color=INK, fontsize=20, labelpad=14)
    fig.text(0.18, 0.027, label["base"], color=MUTED, fontsize=19)
    save_svg(fig, output_dir, f"motor-envelope-{language}", "Illustrative equivalent motor: 1200 Nm, 931 kW, 18000 rpm maximum. Not Tesla specifications.")


def gear_curve(gear):
    # Insert both breakpoints, so the sampled piecewise curves preserve their corners.
    speed = np.sort(np.r_[np.linspace(0, speed_limit(gear), 650), OMEGA_BASE * RADIUS / gear, WHEEL_POWER / GRIP])
    return speed * 3.6, drive_force(speed, gear) / 1000, road_force(speed, gear) / 1000


def render_gearing(output_dir, language):
    configure(language)
    label = labels(language)
    fig = plt.figure(figsize=(6, 9.8), facecolor=PAPER)
    axes = [fig.add_axes([0.16, 0.61, 0.78, 0.23]), fig.add_axes([0.16, 0.27, 0.78, 0.23])]
    fig.text(0.16, 0.945, "G = 9", color=MINT, fontsize=22)
    fig.text(0.57, 0.945, "G = 12", color=BLUE, fontsize=22)
    fig.lines.extend([
        plt.Line2D([0.16, 0.40], [0.927, 0.927], transform=fig.transFigure, color=MINT, lw=3.5),
        plt.Line2D([0.57, 0.84], [0.927, 0.927], transform=fig.transFigure, color=BLUE, lw=3.5, ls=(0, (4, 3))),
    ])
    for ax, title in zip(axes, [label["available"], label["delivered"]]):
        style_axis(ax, title)
        ax.set(xlim=(0, 270), ylim=(0, 46), xticks=[0, 100, 200], yticks=[0, 20, 40])
        ax.axhline(GRIP / 1000, color=AMBER, alpha=0.85, lw=1.2, ls=(0, (2, 3)))
    axes[0].tick_params(labelbottom=False)
    for gear, color, linestyle, width in [(9, MINT, "-", 3.8), (12, BLUE, (0, (4, 3)), 2.7)]:
        speed, drive, delivered = gear_curve(gear)
        for ax, force in zip(axes, [drive, delivered]):
            ax.plot(speed, force, color=color, ls=linestyle, lw=width)
            ax.plot(speed[-1], force[-1], "o" if gear == 9 else "s", color=color, ms=6, mec=PAPER, mew=1)
            ax.plot([speed[-1], speed[-1]], [0, force[-1]], color=color, alpha=0.55, ls=(0, (2, 3)), lw=1)
    axes[0].text(178, 35, label["grip"], color=AMBER, fontsize=18, ha="center")
    axes[1].text(0, 37, label["equal"], color=INK, fontsize=20)
    axes[1].set_xlabel(label["vehicle"], color=INK, fontsize=20, labelpad=14)
    fig.text(0.16, 0.128, label["cutoffs"], color=MUTED, fontsize=19)
    fig.text(0.16, 0.078, label["limit12"], color=BLUE, fontsize=20)
    fig.text(0.16, 0.031, label["limit9"], color=MINT, fontsize=20)
    save_svg(fig, output_dir, f"gearing-envelope-{language}", "Two gear ratios compared before and after tyre grip limits; model only, no drag or rolling resistance. Endpoints are motor-speed cutoffs, not road top speeds.")



CORAL = "#ffab98"
LILAC = "#c4b5ed"


def target_labels(language):
    return {
        "en": {
            "title": "The one-second force gap", "subtitle": "Constant acceleration · G = 9",
            "force": "Force / kN", "speed": "Speed / km/h", "gap": "24.2 kN short",
            "target": "1 s target · 53.64 kN", "tyre": "Tyres · 29.43 kN",
            "torque": "Torque-derived force · 30.81 kN", "power": "Power ceiling · 903 kW",
            "actual": "Delivered · lowest limit", "transition": "110.5 km/h: grip to power",
        },
        "zh": {
            "title": "一秒匀加速的驱动力缺口", "subtitle": "0–60 mph · 减速比 G = 9",
            "force": "驱动力 / kN", "speed": "车速 / km/h", "gap": "缺口 24.2 kN",
            "target": "匀加速目标 · 53.64 kN", "tyre": "轮胎抓地 · 29.43 kN",
            "torque": "转矩换算力 · 30.81 kN", "power": "轮端功率上限 · 903 kW",
            "actual": "实际可用 · 取最低上限", "transition": "110.5 km/h：抓地转为功率限制",
        },
    }[language]


def draw_target_curves(ax, *, annotate=True, language="en", font=22):
    """Plot true coordinates, clipping the singular power curve without flattening it."""
    label = target_labels(language)
    end = launch_requirement()["assumptions"]["end_speed_m_s"] * 3.6
    requirement = launch_requirement()["required_force_N"] / 1000
    transition = WHEEL_POWER / GRIP * 3.6
    speed_kmh = np.sort(np.r_[np.linspace(0, 120, 1000), end, transition,
        OMEGA_BASE * RADIUS / 9 * 3.6])
    available = road_force(speed_kmh / 3.6, 9) / 1000
    target_region = speed_kmh <= end
    ax.fill_between(speed_kmh, 0, available, color=MINT, alpha=0.075, lw=0)
    ax.fill_between(speed_kmh[target_region], available[target_region], requirement,
        color=CORAL, alpha=0.06, lw=0)
    # P/v is not evaluated at zero. Its values above the y-range are clipped by
    # Matplotlib, with no artificial plateau or finite zero-speed power limit.
    positive_kmh = np.linspace(1, 120, 1600)
    ax.plot(positive_kmh, WHEEL_POWER / (positive_kmh / 3.6) / 1000,
        color=BLUE, lw=2.1, ls=(0, (7, 3)), zorder=3)
    ax.axhline(float(drive_force(0, 9)) / 1000, color=LILAC, lw=2,
        ls=(0, (6, 2, 1, 2)), zorder=4)
    ax.plot(speed_kmh, available, color=MINT, lw=5, zorder=5)
    # The tyre line and delivered line coincide until the true transition. Draw
    # the narrow dashed grip line over the thick delivered line without offset.
    ax.axhline(GRIP / 1000, color=AMBER, lw=1.8, ls=(0, (2, 3)), zorder=6)
    ax.plot([0, end], [requirement, requirement], color=CORAL, lw=2.8, zorder=7)
    ax.plot([end, end], [0, requirement], color=CORAL, alpha=0.5,
        ls=(0, (1, 3)), lw=1.3, zorder=2)
    ax.plot(end, requirement, "s", color=CORAL, ms=6, mec=PAPER, mew=1, zorder=8)
    ax.plot([transition, transition], [0, GRIP / 1000], color=MINT,
        alpha=0.45, ls=(0, (1, 3)), lw=1.3, zorder=2)
    ax.plot(transition, GRIP / 1000, "o", ms=7, color=MINT, mec=PAPER, mew=1.5, zorder=9)
    if annotate:
        ax.text(4, 41.2, label["gap"], color=CORAL, fontsize=font,
            bbox={"facecolor": PAPER, "edgecolor": "none", "pad": 2, "alpha": 0.8}, zorder=8)
        ax.text(end, 58.5, "60 mph", color=CORAL, fontsize=font, ha="center", va="center")
        ax.text(114, 11.0, "110.5", color=MINT, fontsize=font, ha="right",
            bbox={"facecolor": PAPER, "edgecolor": "none", "pad": 2, "alpha": 0.8})
    return end, transition


def check_figure_text_bounds(fig):
    """Reject clipped labels in figure assets; visual inspection checks hierarchy."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    for text_artist in fig.findobj(match=matplotlib.text.Text):
        if not text_artist.get_visible() or not text_artist.get_text().strip():
            continue
        box = text_artist.get_window_extent(renderer)
        assert box.x0 >= bounds.x0 - 1 and box.y0 >= bounds.y0 - 1, text_artist.get_text()
        assert box.x1 <= bounds.x1 + 1 and box.y1 <= bounds.y1 + 1, text_artist.get_text()


def render_target(output_dir, language, mobile=False):
    configure(language)
    label = target_labels(language)
    fig = plt.figure(figsize=(6, 10) if mobile else (12, 7.5), facecolor=PAPER)
    fig.text(0.065 if mobile else 0.08, 0.944, label["title"], color=INK,
        fontsize=26 if mobile else 29, weight="medium")
    fig.text(0.065 if mobile else 0.08, 0.895 if mobile else 0.872,
        label["subtitle"], color=MUTED, fontsize=22 if mobile else 20)
    ax = fig.add_axes([0.16, 0.405, 0.77, 0.385] if mobile else [0.08, 0.31, 0.86, 0.47])
    style_axis(ax, "")
    ax.tick_params(labelsize=22 if mobile else 20)
    ax.set(xlim=(0, 120), ylim=(0, 64), xticks=[0, 60, 120] if mobile else [0, 30, 60, 90, 120],
        yticks=[0, 30, 60] if mobile else [0, 15, 30, 45, 60])
    ax.text(0, 1.055, label["force"], transform=ax.transAxes, color=MUTED, fontsize=22 if mobile else 20)
    ax.set_xlabel(label["speed"], color=MUTED, fontsize=22 if mobile else 20, labelpad=8)
    draw_target_curves(ax, language=language, font=22 if mobile else 20)
    items = [
        ("target", CORAL, "-", 2.8), ("tyre", AMBER, (0, (2, 3)), 1.8),
        ("torque", LILAC, (0, (6, 2, 1, 2)), 2), ("power", BLUE, (0, (7, 3)), 2.1),
        ("actual", MINT, "-", 5),
    ]
    positions = [(0.065, 0.294 - i * 0.047) for i in range(5)] if mobile else [
        (0.08, 0.192), (0.53, 0.192), (0.08, 0.122), (0.53, 0.122), (0.08, 0.052)]
    for (key, color, line_style, width), (x, y) in zip(items, positions):
        fig.lines.append(plt.Line2D([x, x + (0.07 if mobile else 0.048)], [y, y],
            transform=fig.transFigure, color=color, lw=width, ls=line_style))
        fig.text(x + (0.095 if mobile else 0.065), y, label[key], color=INK,
            fontsize=22 if mobile else 20, va="center")
    if mobile:
        transition_label = "110.5 km/h:\ngrip to power" if language == "en" else "110.5 km/h：\n抓地转为功率限制"
        fig.text(0.065, 0.043, transition_label, color=MUTED,
            fontsize=22, va="center", linespacing=1.2)
    else:
        fig.text(0.53, 0.052, label["transition"], color=MUTED, fontsize=20, va="center")
    check_figure_text_bounds(fig)
    suffix = f"mobile-{language}" if mobile else language
    save_svg(fig, output_dir, f"target-envelope-{suffix}",
        "Illustrative 2000 kg G=9 drivetrain. The uniform 0-60 mph in one second target requires 53.6448 kN up to 96.56064 km/h. Grip is 29.43 kN; torque force 30.8118 kN; wheel power 903.07 kW. The available force is the minimum; power curve is undefined at rest. Grip-to-power transition is 110.467 km/h. Not Tesla specifications.")


def boundary_layout(language, mobile, title, subtitle, ylabel, xlabel):
    configure(language)
    fig = plt.figure(figsize=(6, 10) if mobile else (12, 7.5), facecolor=PAPER)
    fontsize = 22 if mobile else 20
    fig.text(.07 if mobile else .085, .94, title, color=INK,
        fontsize=25 if mobile else 29, weight="medium", linespacing=1.2, va="top")
    fig.text(.07 if mobile else .085, .845 if mobile else .852, subtitle,
        color=MUTED, fontsize=fontsize, va="top", linespacing=1.25)
    ax = fig.add_axes([.17, .37, .76, .365] if mobile else [.10, .35, .82, .39])
    style_axis(ax, "")
    ax.tick_params(labelsize=fontsize)
    ax.text(0, 1.055, ylabel, transform=ax.transAxes, color=MUTED, fontsize=fontsize)
    ax.set_xlabel(xlabel, color=MUTED, fontsize=fontsize, labelpad=10)
    return fig, ax, fontsize


def point_label(ax, x, y, text, offset, color=INK, fontsize=20):
    ax.plot(x, y, "o", ms=8, color=color, mec=PAPER, mew=1.6, zorder=8, clip_on=False)
    ax.annotate(text, (x, y), xytext=offset, textcoords="offset points",
        ha="left" if offset[0] >= 0 else "right", va="center", color=color,
        fontsize=fontsize, zorder=9,
        bbox={"facecolor": PAPER, "edgecolor": "none", "alpha": .86, "pad": 2})


def render_gear_sensitivity(output_dir, language, mobile=False):
    zh = language == "zh"
    title = "何时值得增大减速比？" if zh else ("When does more\ngearing help?" if mobile else "When does more gearing help?")
    subtitle = "G = 12 相对 G = 9 的用时降幅" if zh else ("Time reduction: G = 12 vs G = 9")
    fig, ax, font = boundary_layout(language, mobile, title, subtitle,
        "可用驱动直流功率 / MW" if zh else "Drive DC power / MW",
        "轮胎附着系数 μ" if zh else "Tyre coefficient μ")
    mu, power = np.meshgrid(np.linspace(1, 2.5, 321), np.linspace(.3, 1.5, 251))
    times = []
    for gear in GEARS:
        cap = np.minimum(mu * MASS * GRAVITY, ETA_GEAR * gear * TORQUE_MAX / RADIUS)
        times.append(envelope_time(MASS, cap, power * 1e6 * ETA_INVERTER * ETA_MOTOR * ETA_GEAR))
    gain = 100 * (1 - times[1] / times[0])
    palette = matplotlib.colors.LinearSegmentedColormap.from_list("gear-gain", ["#23333b", "#336765", "#7eb89a", MINT])
    filled = ax.contourf(mu, power, gain, levels=np.linspace(0, 25, 26), cmap=palette, vmin=0, vmax=25)
    threshold = ETA_GEAR * 9 * TORQUE_MAX / RADIUS / (MASS * GRAVITY)
    ax.axvspan(1, threshold, color="#202a30", alpha=.92, lw=0, zorder=2)
    ax.axvline(threshold, color=AMBER, ls=(0, (3, 3)), lw=1.7, zorder=3)
    contours = ax.contour(mu, power, gain, levels=[5, 10, 15, 20], colors=INK, linewidths=1.1, alpha=.8)
    contour_labels = ax.clabel(contours, fmt=lambda v: f"{int(v)}%", fontsize=font - 1, inline=True, inline_spacing=5)
    for text_artist in contour_labels:
        text_artist.set_bbox({"facecolor": PAPER, "edgecolor": "none", "alpha": .65, "pad": 1})
    ax.set(xlim=(1, 2.5), ylim=(.3, 1.5), xticks=[1, 1.5, 2, 2.5], yticks=[.3, .6, .9, 1.2, 1.5])
    point_label(ax, 1.5, 1, "基线" if zh else "Baseline", (-12, 9), fontsize=font)
    ax.text(1.22, .48, "0%", color=MUTED, fontsize=font + 3, ha="center", zorder=4)
    cax = fig.add_axes([.17, .245, .76, .025] if mobile else [.10, .175, .82, .025])
    colorbar = fig.colorbar(filled, cax=cax, orientation="horizontal", ticks=[0, 10, 20, 25])
    colorbar.ax.tick_params(colors=MUTED, labelsize=font, length=0, pad=4)
    colorbar.outline.set_visible(False)
    colorbar.ax.set_xlabel("用时降幅 / %" if zh else "Time reduction / %", color=MUTED, fontsize=font, labelpad=6)
    note = ("μ ≤ 1.5704：增大减速比无收益\n边界右侧：收益随功率变化" if zh else "μ ≤ 1.5704: no time benefit\nAbove it, benefit depends on power")
    if not mobile:
        note = ("μ ≤ 1.5704：无用时收益；右侧收益大小取决于功率" if zh else "μ ≤ 1.5704: no time benefit; above it, the gain depends on power")
    fig.text(.07 if mobile else .10, .065 if mobile else .045, note,
        fontsize=font, color=AMBER, va="center", linespacing=1.45)
    check_figure_text_bounds(fig)
    suffix = f"mobile-{language}" if mobile else language
    save_svg(fig, output_dir, f"gear-sensitivity-{suffix}",
        "Conditional G12 versus G9 rest-to-60-mph time reduction at fixed 2000 kg, 1200 Nm, 0.34 m wheel radius and efficiencies. Both gears satisfy the rpm cutoff. Exactly zero benefit for mu at or below 1.570426; above that coefficient the gain is positive but depends on DC power. No real tyre capability is asserted.")


def render_downforce_budget(output_dir, language, mobile=False):
    zh = language == "zh"
    title = "下压力的收益能承受多少代价？" if zh else ("What can the\ndownforce gain pay for?" if mobile else "What can the downforce gain pay for?")
    if mobile and zh:
        title = "下压力的收益\n能承受多少代价？"
    subtitle = ("D = 4.905 kN · G = 12 · 总预算 1 MW" if zh else "D = 4.905 kN · G = 12 · 1 MW total")
    if mobile:
        subtitle = "D = 4.905 kN · G = 12\n" + ("总直流预算 1 MW" if zh else "Shared DC budget: 1 MW")
    fig, ax, font = boundary_layout(language, mobile, title, subtitle,
        "辅助系统直流功耗 / kW" if zh else "Auxiliary DC draw / kW",
        "装置新增质量 / kg" if zh else "Added hardware mass / kg")
    # Keep room for the mobile two-line subtitle.
    if mobile:
        ax.set_position([.17, .34, .76, .34])
    mass_grid = np.linspace(0, 900, 800)
    budget = np.array([np.nan if auxiliary_break_even(dm) is None else auxiliary_break_even(dm) / 1000 for dm in mass_grid])
    clipped = np.clip(np.nan_to_num(budget, nan=-1), 0, 600)
    ax.fill_between(mass_grid, 0, clipped, color=MINT, alpha=.14, lw=0)
    ax.fill_between(mass_grid, clipped, 600, color=CORAL, alpha=.09, lw=0)
    valid = budget >= 0
    intercept = decision_boundaries()["downforce_budget"]["zero_auxiliary_break_even_added_mass_kg"]
    boundary_x = np.r_[mass_grid[valid], intercept]
    boundary_y = np.r_[budget[valid], 0]
    ax.plot(boundary_x, boundary_y, color=INK, lw=2.8, zorder=5)
    ax.set(xlim=(0, 900), ylim=(0, 600), xticks=[0, 300, 600, 900], yticks=[0, 200, 400, 600])
    point_label(ax, 0, 0, "A", (10, 14), color=MINT, fontsize=font)
    point_label(ax, 100, 100, "B", (10, 3), color=MINT, fontsize=font)
    point_label(ax, 500, 300, "C", (10, 12), color=CORAL, fontsize=font)
    ax.text(410, 85, "更快" if zh else "Faster", ha="center", color=MINT, fontsize=font)
    ax.text(560, 515, "更慢" if zh else "Slower", ha="center", color=CORAL, fontsize=font)
    if mobile:
        notes = [("A  1.464 s   ·   B  1.512 s", MINT), ("C  1.891 s", CORAL),
            ("白线：与基线 1.823 s 持平" if zh else "White line: baseline 1.823 s", INK),
            ("D 固定；并非风扇性能预测" if zh else "Fixed D; not a fan prediction", MUTED)]
        for i, (text, color) in enumerate(notes):
            fig.text(.07, .213 - i * .051, text, fontsize=font, color=color, va="center")
    else:
        fig.text(.10, .18, "A  1.464 s     B  1.512 s", fontsize=font, color=MINT)
        fig.text(.64, .18, "C  1.891 s", fontsize=font, color=CORAL)
        fig.text(.10, .106, "白线：与原始基线 1.823 s 持平" if zh else "White line: equal to the original 1.823 s baseline", fontsize=font, color=INK)
        fig.text(.10, .042, "D 固定；图示是成本预算，不代表风扇能提供该下压力" if zh else "D is prescribed: this is a cost budget, not a fan-performance prediction", fontsize=font, color=MUTED)
    check_figure_text_bounds(fig)
    suffix = f"mobile-{language}" if mobile else language
    save_svg(fig, output_dir, f"downforce-budget-{suffix}",
        "Conditional added-mass and auxiliary-power budget for fixed 4905 N downforce from rest and G12. Shared DC power is 1 MW, auxiliary draw is subtracted before drivetrain efficiencies. Zero-gain curve is relative to original 2000 kg G9 no-downforce baseline 1.822793 s. A=(0 kg,0 kW), B=(100 kg,100 kW), C=(500 kg,300 kW). Downforce is not derived from the hardware cost; no fan feasibility is asserted.")


def render_one_second(output_dir, language, mobile=False):
    zh = language == "zh"
    title = "一秒需要哪些能力同时成立？" if zh else ("Which capacities\nmake one second possible?" if mobile else "Which capacities make one second possible?")
    if mobile and zh:
        title = "一秒需要哪些能力\n同时成立？"
    subtitle = "2,000 kg · 0–60 mph · 能力需求图" if zh else "2,000 kg · 0–60 mph · Capacity requirements"
    if mobile:
        subtitle = "2,000 kg · 0–60 mph\n" + ("能力需求，不是既有硬件" if zh else "Required, not existing hardware")
    fig, ax, font = boundary_layout(language, mobile, title, subtitle,
        "轮端功率上限 P / MW" if zh else "Wheel-power ceiling P / MW",
        "有效驱动力上限 C / kN" if zh else "Effective force cap C / kN")
    if mobile:
        ax.set_position([.17, .34, .76, .34])
    force_min = MASS * 26.8224
    capacity = np.r_[force_min, np.linspace(force_min + 1, 100_000, 1000)]
    power = np.array([required_wheel_power(MASS, c, 1) for c in capacity]) / 1e6
    ax.axvspan(25, force_min / 1000, color=CORAL, alpha=.08, lw=0)
    ax.fill_between(capacity / 1000, .4, power, color=CORAL, alpha=.08, lw=0)
    ax.fill_between(capacity / 1000, power, 1.8, color=MINT, alpha=.14, lw=0)
    ax.plot(capacity / 1000, power, color=MINT, lw=3)
    ax.plot([force_min / 1000] * 2, [power[0], 1.8], color=MINT, lw=3)
    ax.axvline(force_min / 1000, color=MUTED, ls=(0, (2, 4)), lw=1)
    ax.axhline(MASS * 26.8224**2 / 2 / 1e6, color=GRID, ls=(0, (2, 4)), lw=1)
    ax.set(xlim=(25, 100), ylim=(.4, 1.8), xticks=[25, 50, 75, 100], yticks=[.4, .8, 1.2, 1.6])
    point_label(ax, force_min / 1000, power[0], "A", (-12, 8), color=AMBER, fontsize=font)
    point_label(ax, 70, required_wheel_power(MASS, 70_000, 1) / 1e6, "B", (9, -8), color=BLUE, fontsize=font)
    # Original and zero-cost combined example both fail the minimum force cap.
    point_label(ax, GRIP / 1000, WHEEL_POWER / 1e6, "0", (-4, 19), color=MUTED, fontsize=font)
    point_label(ax, MU * (MASS * GRAVITY + 4905) / 1000, WHEEL_POWER / 1e6, "+", (7, -18), color=MUTED, fontsize=font)
    ax.text(79, 1.53, "t ≤ 1 s", color=MINT, fontsize=font + 2, ha="center")
    ax.text(81, .51, "t > 1 s", color=CORAL, fontsize=font + 2, ha="center")
    if mobile:
        # Numeric capacity pairs use a separate line to stay readable at 320 px.
        notes = [("A  53.64 kN / 1.439 MW", AMBER), ("B  70 kN / 0.876 MW", BLUE),
            ("0 基线 · + 下压力与 G = 12" if zh else "0 Baseline · + Downforce & G12", MUTED),
            ("C 同时受轮胎与转矩限制" if zh else "C needs both grip and torque", INK)]
        for i, (text, color) in enumerate(notes):
            fig.text(.07, .213 - i * .051, text, fontsize=font, color=color, va="center")
    else:
        fig.text(.10, .18, "A  53.64 kN / 1.439 MW", fontsize=font, color=AMBER)
        fig.text(.55, .18, "B  70 kN / 0.876 MW", fontsize=font, color=BLUE)
        fig.text(.10, .106, "0 原始基线 · + 下压力与 G = 12；两者均未达到一秒" if zh else "0 Original baseline · + Downforce & G12; neither reaches one second", fontsize=font, color=MUTED)
        fig.text(.10, .042, "C 必须同时得到轮胎抓地与转矩换算力的支持" if zh else "C must be supported by both tyre grip and torque-derived wheel force", fontsize=font, color=INK)
    check_figure_text_bounds(fig)
    suffix = f"mobile-{language}" if mobile else language
    save_svg(fig, output_dir, f"one-second-feasibility-{suffix}",
        "Required independent effective force and wheel-power capacities for an ideal 2000 kg rest-to-60-mph launch in one second. Force cap must be at least 53.6448 kN. At that cap, 1.438882 MW is needed; with 70 kN, 876.075 kW suffices using varying acceleration. The current assumed 1200 Nm motor is not held fixed: both grip and torque force must support C. Baseline and zero-cost downforce/G12 combination are below the force requirement.")


def render_cover(output_dir):
    configure("en")
    fig = plt.figure(figsize=(12, 6.3), facecolor=PAPER)
    fig.text(0.085, 0.91, "ROADSTER / 02", color=MUTED, fontsize=20)
    fig.text(0.915, 0.91, "ILLUSTRATIVE MODEL", color=MUTED, fontsize=16, ha="right")
    fig.text(0.085, 0.78, "53.6 kN", color=CORAL, fontsize=42)
    fig.text(0.375, 0.79, "REQUIRED", color=MUTED, fontsize=18)
    fig.text(0.58, 0.78, "29.4 kN", color=MINT, fontsize=42)
    fig.text(0.87, 0.79, "AT LAUNCH", color=MUTED, fontsize=16)
    ax = fig.add_axes([0.085, 0.155, 0.83, 0.50])
    style_axis(ax, "")
    ax.set(xlim=(0, 120), ylim=(0, 64), xticks=[0, 30, 60, 90, 120], yticks=[0, 30, 60])
    end = 26.8224 * 3.6
    speed = np.sort(np.r_[np.linspace(0, 120, 900), end, WHEEL_POWER / GRIP * 3.6])
    actual = road_force(speed / 3.6, 9) / 1000
    ax.fill_between(speed[speed <= end], actual[speed <= end], 53.6448, color=CORAL, alpha=0.07)
    ax.fill_between(speed, 0, actual, color=MINT, alpha=0.075)
    ax.plot([0, end], [53.6448, 53.6448], color=CORAL, lw=3)
    ax.plot([end, end], [0, 53.6448], color=CORAL, ls=(0, (2, 4)), alpha=0.6, lw=1.3)
    ax.plot(end, 53.6448, "s", color=CORAL, ms=8)
    ax.plot(speed, actual, color=MINT, lw=4)
    ax.text(3, 58.5, "0–60 mph IN 1 s · CONSTANT ACCELERATION", color=CORAL, fontsize=18)
    ax.text(3, 18, "AVAILABLE FORCE", color=MINT, fontsize=18)
    ax.text(end + 1, 43, "60 mph", color=CORAL, fontsize=18)
    ax.set_xlabel("SPEED / km/h", color=MUTED, fontsize=18, labelpad=8)
    ax.text(0, 1.04, "FORCE / kN", transform=ax.transAxes, color=MUTED, fontsize=18)
    check_figure_text_bounds(fig)
    fig.savefig(output_dir / "cover.png", dpi=100, metadata={"Software": "FENG · Matplotlib",
        "Description": "An assumed one-second constant-acceleration target needs 53.6 kN, while the illustrative G=9 drivetrain delivers 29.4 kN at launch. Not Tesla specifications."})
    plt.close(fig)


def summary():
    return {
        "kind": "illustrative-equivalent-drivetrain-not-tesla-specifications",
        "assumptions": {"mass_kg": MASS, "g_m_s2": GRAVITY, "mu": MU, "radius_m": RADIUS,
            "torque_max_Nm": TORQUE_MAX, "rpm_max": RPM_MAX, "pack_dc_W": DC_POWER,
            "inverter_efficiency": ETA_INVERTER, "motor_efficiency": ETA_MOTOR, "gear_efficiency": ETA_GEAR},
        "shaft_power_W": SHAFT_POWER, "wheel_power_W": WHEEL_POWER,
        "loss_W": {"inverter": 20_000, "motor": 49_000, "gear": 27_930},
        "motor_base_rpm": RPM_BASE, "tyre_grip_N": GRIP,
        "grip_to_power_kmh": WHEEL_POWER / GRIP * 3.6,
        "gears": [{"ratio": gear, "available_torque_force_N": float(drive_force(0, gear)),
            "motor_base_kmh": OMEGA_BASE * RADIUS / gear * 3.6,
            "motor_speed_cutoff_kmh": speed_limit(gear) * 3.6} for gear in GEARS],
        "examples": [{"speed_kmh": speed * 3.6, "road_force_N": float(road_force(speed, 9)),
            "acceleration_m_s2": float(road_force(speed, 9)) / MASS} for speed in [0, 20, 40, 50]],
        "launch_requirement": launch_requirement(),
        "nominal_pack_comparison": nominal_pack_comparison(),
        "design_comparison": design_comparison(),
        "decision_boundaries": decision_boundaries(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1] / "content" / "assets" / "roadster-power")
    args = parser.parse_args()
    check_model()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for language in ["en", "zh"]:
        render_motor(args.output_dir, language)
        render_gearing(args.output_dir, language)
        render_target(args.output_dir, language)
        render_target(args.output_dir, language, mobile=True)
        for mobile in [False, True]:
            render_gear_sensitivity(args.output_dir, language, mobile)
            render_downforce_budget(args.output_dir, language, mobile)
            render_one_second(args.output_dir, language, mobile)
    render_cover(args.output_dir)
    (args.output_dir / "model.json").write_text(json.dumps(summary(), ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary(), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
