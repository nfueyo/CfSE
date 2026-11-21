"""Sweep elbow inlet speed in Fluent and calculate total-pressure-loss K.

Example (PowerShell):
  python elbow_k_Re_calculation.py --case "C:\\path\\elbow_no_h2.cas.h5" `
      --velocities 0.2 0.4 0.6 --iterations 500

The supplied case uses hot-air-inlet as the swept inlet, sets h2-inlet to zero,
and uses outlet as the outlet. Re uses volume-average fluid density and laminar
viscosity with the hot-air hydraulic diameter. K = (Pt,in - Pt,out) /
(0.5 * rho_avg * U_in**2), where Pt values are mass-weighted averages.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import uuid
from datetime import datetime
from pathlib import Path


def report_scalar(result, name):
    """Extract one finite scalar from a Fluent report-compute response."""
    if isinstance(result, dict):
        if name in result:
            value = result[name]
            if isinstance(value, (list, tuple)):
                if not value:
                    raise ValueError(f"Empty report result for {name!r}")
                value = value[0]
            value = float(value)
            if not math.isfinite(value):
                raise ValueError(f"Non-finite report result for {name!r}: {value}")
            return value
        for nested in result.values():
            try:
                return report_scalar(nested, name)
            except KeyError:
                continue
    elif isinstance(result, (list, tuple)):
        for nested in result:
            try:
                return report_scalar(nested, name)
            except KeyError:
                continue
    raise KeyError(f"Cannot find report {name!r} in Fluent response: {result!r}")


def mass_imbalance_percent(mass_flows):
    """Return absolute net boundary flow as a percent of total absolute flow."""
    total_absolute_flow = sum(abs(flow) for flow in mass_flows)
    if total_absolute_flow <= 0 or not math.isfinite(total_absolute_flow):
        raise ValueError("Cannot calculate mass imbalance: no finite boundary mass flow.")
    net_flow = sum(mass_flows)
    return 100.0 * abs(net_flow) / total_absolute_flow


def calculate_row(velocity, iterations, cumulative_iterations, inlet_total_pressure,
                  outlet_total_pressure, mass_flows, density, viscosity,
                  hydraulic_diameter):
    if density <= 0 or viscosity <= 0 or hydraulic_diameter <= 0:
        raise ValueError("Density, viscosity, and hydraulic diameter must be positive.")
    if velocity <= 0:
        raise ValueError("Swept inlet velocity must be positive.")
    pressure_loss = inlet_total_pressure - outlet_total_pressure
    dynamic_pressure = 0.5 * density * velocity**2
    return {
        "inlet_velocity_m_s": velocity,
        "inlet_total_pressure_Pa": inlet_total_pressure,
        "outlet_total_pressure_Pa": outlet_total_pressure,
        "total_pressure_loss_Pa": pressure_loss,
        "mass_imbalance_percent": mass_imbalance_percent(mass_flows),
        "iterations_this_velocity": iterations,
        "total_iterations_performed": cumulative_iterations,
        "volume_avg_density_kg_m3": density,
        "volume_avg_laminar_viscosity_Pa_s": viscosity,
        "hydraulic_diameter_m": hydraulic_diameter,
        "Re": density * velocity * hydraulic_diameter / viscosity,
        "k": pressure_loss / dynamic_pressure,
    }


def configure_reports(settings, inlet_name, disabled_inlet_name, outlet_name,
                      fluid_zone, prefix):
    definitions = settings.solution.report_definitions
    report_names = {}

    for key, boundary in (("inlet_total_pressure", inlet_name),
                          ("outlet_total_pressure", outlet_name)):
        name = f"{prefix}-{key}"
        definitions.surface.create(name=name)
        report = definitions.surface[name]
        report.report_type.set_state("surface-massavg")
        report.field.set_state("total-pressure")
        report.surface_names.set_state([boundary])
        report_names[key] = name

    for key, field in (("volume_avg_density", "density"),
                       ("volume_avg_viscosity", "viscosity-lam")):
        name = f"{prefix}-{key}"
        definitions.volume.create(name=name)
        report = definitions.volume[name]
        report.report_type.set_state("volume-average")
        report.field.set_state(field)
        report.cell_zones.set_state([fluid_zone])
        report_names[key] = name

    for key, boundary in (("inlet_mass_flow", inlet_name),
                          ("disabled_inlet_mass_flow", disabled_inlet_name),
                          ("outlet_mass_flow", outlet_name)):
        name = f"{prefix}-{key}"
        definitions.flux.create(name=name)
        definitions.flux[name].boundaries.set_state([boundary])
        report_names[key] = name
    return report_names


def compute_row(settings, report_names, velocity, iterations,
                cumulative_iterations, hydraulic_diameter):
    response = settings.solution.report_definitions.compute(
        report_defs=list(report_names.values()))
    if response is None:
        raise RuntimeError("Fluent returned no report results.")
    values = {key: report_scalar(response, name)
              for key, name in report_names.items()}
    mass_flows = [
        values["inlet_mass_flow"],
        values["disabled_inlet_mass_flow"],
        values["outlet_mass_flow"],
    ]
    return calculate_row(
        velocity=velocity,
        iterations=iterations,
        cumulative_iterations=cumulative_iterations,
        inlet_total_pressure=values["inlet_total_pressure"],
        outlet_total_pressure=values["outlet_total_pressure"],
        mass_flows=mass_flows,
        density=values["volume_avg_density"],
        viscosity=values["volume_avg_viscosity"],
        hydraulic_diameter=hydraulic_diameter,
    )


def write_results(rows, folder):
    if not rows:
        return
    with (folder / "results.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    columns = list(rows[0])
    lines = [" ".join(f"{column:>35}" for column in columns)]
    for row in rows:
        lines.append(" ".join(f"{row[column]:35.8g}" for column in columns))
    (folder / "tables.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines), flush=True)


def plot_results(rows, folder):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ordered = sorted(rows, key=lambda row: row["Re"])
    figure, axis = plt.subplots(figsize=(8, 6))
    axis.plot([row["Re"] for row in ordered],
              [row["k"] for row in ordered], marker="o")
    axis.set(xlabel="Reynolds number, Re", ylabel="Pressure-loss constant, k",
             title="Elbow pressure-loss constant vs Reynolds number")
    axis.grid(True, alpha=0.3)
    figure.tight_layout()
    figure.savefig(folder / "k_vs_Re.png", dpi=180)
    figure.savefig(folder / "k_vs_Re.pdf")
    plt.close(figure)


def run_sweep(args, folder):
    import ansys.fluent.core as pyfluent

    kwargs = {
        "mode": "solver",
        "precision": "double",
        "dimension": 3,
        "processor_count": args.processors,
        "cwd": str(folder),
    }
    if args.version:
        kwargs["product_version"] = args.version

    solver = pyfluent.launch_fluent(**kwargs)
    try:
        settings = solver.settings
        settings.file.read_case(file_name=str(args.case.resolve()))
        if solver.scheme.eval("(rpgetvar 'rp-unsteady?)"):
            raise RuntimeError("This script runs steady iterations; the supplied case is transient.")

        boundary_conditions = settings.setup.boundary_conditions
        inlets = boundary_conditions.velocity_inlet.get_object_names()
        outlets = boundary_conditions.pressure_outlet.get_object_names()
        for name in (args.inlet, args.disabled_inlet):
            if name not in inlets:
                raise ValueError(f"Velocity inlet {name!r} not found. Available: {inlets}")
        if args.outlet not in outlets:
            raise ValueError(f"Pressure outlet {args.outlet!r} not found. Available: {outlets}")
        fluid_zones = settings.setup.cell_zone_conditions.fluid.get_object_names()
        if args.fluid_zone not in fluid_zones:
            raise ValueError(f"Fluid cell zone {args.fluid_zone!r} not found. Available: {fluid_zones}")

        prefix = "elbow-" + uuid.uuid4().hex[:8]
        report_names = configure_reports(
            settings, args.inlet, args.disabled_inlet, args.outlet,
            args.fluid_zone, prefix)
        manifest = {
            "case": str(args.case.resolve()),
            "fluent_version": str(solver.get_fluent_version()),
            "pyfluent_version": pyfluent.__version__,
            "inlet_swept": args.inlet,
            "inlet_disabled": args.disabled_inlet,
            "outlet": args.outlet,
            "fluid_zone": args.fluid_zone,
            "hydraulic_diameter_m": args.hydraulic_diameter,
            "velocities_m_s": args.velocities,
            "iterations_per_velocity": args.iterations,
            "reynolds_definition": "volume_avg_density * inlet_velocity * hydraulic_diameter / volume_avg_laminar_viscosity",
            "k_definition": "(mass_avg_inlet_total_pressure - mass_avg_outlet_total_pressure) / (0.5 * volume_avg_density * inlet_velocity**2)",
            "mass_imbalance_definition": "100 * abs(sum(signed boundary mass flows)) / sum(abs(boundary mass flows))",
        }
        (folder / "configuration.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8")

        rows = []
        total_iterations = 0
        for velocity in args.velocities:
            print(f"\nSweeping {args.inlet} at {velocity:g} m/s", flush=True)
            boundary_conditions.velocity_inlet[
                args.inlet].momentum.velocity_magnitude.value = velocity
            boundary_conditions.velocity_inlet[
                args.disabled_inlet].momentum.velocity_magnitude.value = 0.0
            settings.solution.initialization.hybrid_initialize()
            settings.solution.run_calculation.iterate(iter_count=args.iterations)
            total_iterations += args.iterations
            rows.append(compute_row(
                settings, report_names, velocity, args.iterations,
                total_iterations, args.hydraulic_diameter))
            write_results(rows, folder)

        plot_results(rows, folder)
        manifest["status"] = "complete"
        manifest["total_iterations_performed"] = total_iterations
        (folder / "configuration.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8")
    except Exception as exc:
        (folder / "failure.txt").write_text(repr(exc), encoding="utf-8")
        raise
    finally:
        solver.exit()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True,
                        help="Fluent case file (.cas.h5 or .cas)")
    parser.add_argument("--velocities", nargs="+", type=float, required=True,
                        help="Hot-air inlet velocities in m/s")
    parser.add_argument("--iterations", type=int, default=500,
                        help="Steady iterations per velocity (default: 500)")
    parser.add_argument("--processors", type=int, default=4)
    parser.add_argument("--version", help="Fluent version, e.g. 26.1")
    parser.add_argument("--output", type=Path, default=Path("elbow_k_results"),
                        help="Output root; each run gets a unique subfolder")
    parser.add_argument("--inlet", default="hot-air-inlet")
    parser.add_argument("--disabled-inlet", default="h2-inlet")
    parser.add_argument("--outlet", default="outlet")
    parser.add_argument("--fluid-zone", default="fluid")
    parser.add_argument("--hydraulic-diameter", type=float, default=0.1016,
                        help="Reference hydraulic diameter in m (default: 0.1016)")
    args = parser.parse_args(argv)

    if not args.case.is_file():
        parser.error(f"Case file does not exist: {args.case}")
    if not args.velocities or any(not math.isfinite(v) or v <= 0 for v in args.velocities):
        parser.error("All sweep velocities must be finite and positive.")
    if len(set(args.velocities)) != len(args.velocities):
        parser.error("Do not repeat sweep velocities.")
    if (args.iterations < 1 or args.processors < 1
            or not math.isfinite(args.hydraulic_diameter)
            or args.hydraulic_diameter <= 0):
        parser.error("Iterations, processors, and hydraulic diameter must be positive.")

    run_folder = args.output.resolve() / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
    run_folder.mkdir(parents=True, exist_ok=False)
    run_sweep(args, run_folder)
    print(f"\nResults saved to {run_folder}")


if __name__ == "__main__":
    main()