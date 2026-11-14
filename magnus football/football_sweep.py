"""Steady Fluent football sweep. SI inputs; see football_sweep_README.md.

Run on a machine with licensed Fluent and ansys-fluent-core installed.
No Fluent installation is needed for --help or --plot-only.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime
from pathlib import Path


DRAG = "sweep-drag-force"
MAGNUS = "sweep-magnus-force"


def nondimensional(U, omega, drag, magnus, D, rho, mu):
    """Magnus is the signed force along fixed -z, for all spin signs."""
    area = math.pi * D**2 / 4
    force_scale = 0.5 * rho * U**2 * area
    return dict(
        U_m_s=U, omega_rad_s=omega, rpm=omega * 60 / (2 * math.pi),
        equatorial_speed_m_s=omega * D / 2,
        drag_N=drag, magnus_minus_z_N=magnus,
        Re=rho * U * D / mu, S=omega * D / (2 * U),
        C_D=drag / force_scale, C_M=magnus / force_scale,
    )


def report_scalar(result, name):
    """Extract by name; never mistake report metadata for the force."""
    if isinstance(result, dict) and name in result:
        value = result[name]
        while isinstance(value, (list, tuple)):
            if not value:
                raise ValueError(f"Empty report: {name}")
            value = value[0]
        value = float(value)
        if not math.isfinite(value):
            raise ValueError(f"Non-finite force for {name}: {value}")
        return value
    if isinstance(result, (list, tuple)):
        for item in result:
            try:
                return report_scalar(item, name)
            except KeyError:
                pass
    raise KeyError(f"Cannot find {name!r} in compute response: {result!r}")


def write_tables(rows, folder):
    with (folder / "results.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    groups = [
        ("DIMENSIONAL VALUES", ["U_m_s", "omega_rad_s", "rpm", "drag_N", "magnus_minus_z_N"]),
        ("DIMENSIONLESS VALUES", ["Re", "S", "C_D", "C_M"]),
    ]
    lines = []
    for title, columns in groups:
        lines += [title, " ".join(f"{c:>19}" for c in columns)]
        for row in rows:
            lines.append(" ".join(f"{row[c]:19.8g}" for c in columns))
        lines.append("")
    text = "\n".join(lines)
    (folder / "tables.txt").write_text(text, encoding="utf-8")
    print(text)


def plot_results(rows, folder):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Us = sorted({r["U_m_s"] for r in rows})
    spins = sorted({r["omega_rad_s"] for r in rows})
    lookup = {(r["U_m_s"], r["omega_rad_s"]): r for r in rows}
    for coefficient, title in [("C_D", "Drag coefficient"), ("C_M", "Signed Magnus coefficient (-z)")]:
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111, projection="3d")
        if len(Us) >= 2 and len(spins) >= 2 and len(lookup) == len(Us) * len(spins):
            # The Cartesian grid is in U and omega, NOT in Re and S.
            # Preserve its connectivity after transforming coordinates.
            Re = np.array([[lookup[u, w]["Re"] for w in spins] for u in Us])
            S = np.array([[lookup[u, w]["S"] for w in spins] for u in Us])
            C = np.array([[lookup[u, w][coefficient] for w in spins] for u in Us])
            ax.plot_surface(Re, S, C, cmap="viridis", alpha=0.8, edgecolor="0.5")
        ax.scatter([r["Re"] for r in rows], [r["S"] for r in rows],
                   [r[coefficient] for r in rows], color="black", s=15)
        ax.set(xlabel="Re = rho U D / mu", ylabel="S = omega D / (2 U)",
               zlabel=coefficient, title=title + " (steady iteration endpoints)")
        fig.tight_layout()
        fig.savefig(folder / f"{coefficient}_3D.png", dpi=200)
        fig.savefig(folder / f"{coefficient}_3D.pdf")
        plt.close(fig)


def set_force_report(collection, name, direction, output_type):
    collection.create(name=name)
    report = collection[name]
    report.zones = ["ball"]  # Exterior side only; never include ball-shadow.
    report.force_vector = direction
    report.average_over = 1
    # Current releases use report_output_type; older ones expose scaled.
    children = report.get_state()
    if "report_output_type" in children:
        report.report_output_type = output_type
    elif "scaled" in children:
        report.scaled = False
    else:
        raise RuntimeError(f"Cannot select dimensional force output: {children}")
    return report


def check_case(solver, settings):
    if solver.scheme_eval.scheme_eval("(rpgetvar 'rp-unsteady?)"):
        raise RuntimeError("This script uses steady iterations. Select steady mode first, or write a transient time-averaging workflow.")
    bc = settings.setup.boundary_conditions
    if "inlet" not in bc.velocity_inlet.get_object_names():
        raise RuntimeError("Expected velocity inlet named 'inlet'.")
    for zone in ("ball", "ball-shadow"):
        if zone not in bc.wall.get_object_names():
            raise RuntimeError(f"Missing wall zone {zone!r}.")
        m = bc.wall[zone].momentum
        if not m.rotating.get_state():
            raise RuntimeError(f"{zone}: enable rotational moving-wall motion first.")
        axis = m.rotation_axis_direction.get_state()
        origin = m.rotation_axis_origin.get_state()
        if any(not math.isclose(float(a), b, abs_tol=1e-10) for a, b in zip(axis, [0, 1, 0])):
            raise RuntimeError(f"{zone}: expected +y rotation axis, found {axis}.")
        if any(not math.isclose(float(a), 0, abs_tol=1e-10) for a in origin):
            raise RuntimeError(f"{zone}: expected rotation origin (0,0,0), found {origin}.")
    # Existing case must have stationary cell zones, as in OG_ball_DOM.cas.h5.


def run(args, folder):
    import ansys.fluent.core as pyfluent

    kwargs = dict(mode="solver", precision="double", dimension=3,
                  processor_count=args.processors, cwd=str(folder))
    if args.version:
        kwargs["product_version"] = args.version
    solver = pyfluent.launch_fluent(**kwargs)
    try:
        s = solver.settings
        s.file.read_case(file_name=str(args.case.resolve()))
        check_case(solver, s)
        reports = s.solution.report_definitions
        set_force_report(reports.drag, DRAG, [1, 0, 0], "Drag Force")
        set_force_report(reports.lift, MAGNUS, [0, 0, -1], "Lift Force")
        manifest = {
            "case": str(args.case.resolve()), "diameter_m": args.diameter,
            "rho_kg_m3": args.rho, "mu_Pa_s": args.mu,
            "velocities_m_s": args.velocities, "spins_rad_s": args.spins,
            "iterations_requested_per_point": args.iterations,
            "drag_direction": [1, 0, 0], "magnus_direction": [0, 0, -1],
            "initialization": "hybrid at every point", "pyfluent": pyfluent.__version__,
            "meaning": "Steady iteration endpoints; convergence must be checked separately.",
        }
        (folder / "configuration.json").write_text(json.dumps(manifest, indent=2))
        # Save settings for version-specific troubleshooting.
        (folder / "wall_settings.json").write_text(json.dumps(
            s.setup.boundary_conditions.wall["ball"].momentum.get_state(), indent=2))
        rows = []
        for i, U in enumerate(args.velocities):
            for j, omega in enumerate(args.spins):
                print(f"\nU = {U:g} m/s; omega = {omega:g} rad/s", flush=True)
                bc = s.setup.boundary_conditions
                # The supplied case already uses an inlet normal to its -x plane.
                bc.velocity_inlet["inlet"].momentum.velocity_magnitude.value = U
                for zone in ("ball", "ball-shadow"):
                    bc.wall[zone].momentum.rotation_speed.value = omega
                s.solution.initialization.hybrid_initialize()
                # Test force-return compatibility before a potentially long run.
                preflight = reports.compute(report_defs=[DRAG, MAGNUS])
                if preflight is None:
                    raise RuntimeError("Report compute returned None. Use a supported PyFluent/Fluent combination (prefer Fluent 2024 R2 or later).")
                report_scalar(preflight, DRAG)
                report_scalar(preflight, MAGNUS)
                s.solution.run_calculation.iterate(iter_count=args.iterations)
                result = reports.compute(report_defs=[DRAG, MAGNUS])
                row = nondimensional(U, omega, report_scalar(result, DRAG),
                                     report_scalar(result, MAGNUS),
                                     args.diameter, args.rho, args.mu)
                rows.append(row)
                # Persist after every point so a later failure loses no completed rows.
                write_tables(rows, folder)
                if args.save_solutions:
                    s.file.write_case_data(file_name=str(folder / f"point_{i:03d}_{j:03d}.cas.h5"))
        plot_results(rows, folder)
    finally:
        solver.exit()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", type=Path)
    p.add_argument("--velocities", nargs="+", type=float, default=[15, 20, 25, 30])
    p.add_argument("--spins", nargs="+", type=float, default=[-100, -50, 0, 50, 100])
    p.add_argument("--iterations", type=int, default=1000)
    p.add_argument("--diameter", type=float, default=0.229)
    p.add_argument("--rho", type=float, default=1.225)
    p.add_argument("--mu", type=float, default=1.7894e-5)
    p.add_argument("--processors", type=int, default=4)
    p.add_argument("--version", help="Installed Fluent version, e.g. 24.2 or 25.2; default: latest installed")
    p.add_argument("--output", type=Path, default=Path("football_results"))
    p.add_argument("--save-solutions", action="store_true")
    p.add_argument("--plot-only", type=Path, metavar="RESULTS_CSV")
    args = p.parse_args()
    if args.plot_only:
        with args.plot_only.open(newline="") as f:
            rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
        if not rows:
            p.error("CSV is empty")
        plot_results(rows, args.plot_only.resolve().parent)
        return
    if args.case is None or not args.case.is_file():
        p.error("Provide --case with the existing OG_ball_DOM.cas.h5 path")
    for name in ("diameter", "rho", "mu"):
        if not math.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
            p.error(f"--{name} must be finite and positive")
    if args.iterations < 1 or args.processors < 1:
        p.error("Iterations and processors must be positive")
    if any(not math.isfinite(u) or u <= 0 for u in args.velocities):
        p.error("All inlet speeds must be finite and positive (normalization is undefined at U=0)")
    if any(not math.isfinite(w) for w in args.spins):
        p.error("All spin speeds must be finite")
    if len(set(args.velocities)) != len(args.velocities) or len(set(args.spins)) != len(args.spins):
        p.error("Do not repeat sweep values")
    # Unique run folder prevents overwriting an earlier study.
    folder = args.output.resolve() / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    folder.mkdir(parents=True, exist_ok=False)
    run(args, folder)
    print(f"Results saved in {folder}")


if __name__ == "__main__":
    main()
