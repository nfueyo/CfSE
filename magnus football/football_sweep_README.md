# Football velocity/spin sweep

`football_sweep.py` implements a steady parameter sweep for `OG_ball_DOM.cas.h5`.
Run it on Windows or Linux with a licensed Fluent installation. The script was
checked locally for Python syntax, force-response extraction, and normalization;
it has **not** been executed against Fluent. Boundary APIs can vary by release.
Prefer Fluent 2024 R2 or later with a compatible current PyFluent installation.
The original case was written by Fluent 2023 R1; that is the case-file version,
not a requirement to run the script in that version. Some older combinations
return `None` from report computation; the script detects this before iterating.

## Install and run

In a Python virtual environment on your Fluent computer:

```sh
python -m pip install ansys-fluent-core numpy matplotlib
```

First try ONE point with a small iteration budget to verify your installation
and API, then run the full sweep. These short runs are not converged studies.

```sh
python football_sweep.py --case "C:/cases/OG_ball_DOM.cas.h5" --velocities 25 --spins 50 --iterations 10
```

Example full sweep (edit the speeds, number of iterations, and case path):

```sh
python football_sweep.py --case "C:/cases/OG_ball_DOM.cas.h5" --velocities 15 20 25 30 --spins -100 -50 0 50 100 --iterations 1000 --processors 4 --save-solutions
```

This requests 20 calculations of up to 1000 iterations each. Fluent may stop
earlier when its convergence criteria are satisfied. A prescribed iteration
budget is not proof of convergence: inspect residuals and force stabilization,
and repeat with more iterations where needed. Optional case/data files preserve
each point for inspection. With `--save-solutions` omitted, only tables, metadata
and plots are retained, not restartable solution fields.

`--version 25.2`, for example, selects an installed Fluent version. Otherwise
PyFluent selects the latest available installation. Run `--help` for all options.

## Case assumptions and quantities

- Inlet `inlet`: normal inflow on its x = -2.5 m plane, with U > 0 along +x.
- Walls `ball` and `ball-shadow`: already configured as rotating moving walls
  about +y through the origin, with stationary adjacent fluid zones.
- The script changes the two rotation speeds to the SAME signed value. It does
  not reverse the shadow speed: angular velocity is a vector independent of
  which side of the wall is considered.
- Only `ball` contributes to the aerodynamic reports. The shadow side bounds
  the interior air zone and must not be included in the exterior force sum.
- Drag: total pressure plus viscous force along +x, in newtons.
- Magnus: signed total pressure plus viscous force along fixed -z, in newtons.
  This direction is +y cross +x. The direction remains fixed for negative spin,
  so reversing spin can reverse the reported coefficient. The script does not
  take an absolute value or impose the expected force sign. A nonzero lateral
  force at zero spin may reflect geometry asymmetry or numerical error; it is
  not automatically a pure spin-induced increment. No zero-spin subtraction is
  performed.
- The script creates two fresh dimensional force reports, rather than relying
  on possibly mis-scaled existing coefficient reports or Workbench parameters.
- The inherited SST model, numerical schemes, mesh, and turbulence inlet
  specification are retained. It does not silently switch a transient case to
  steady; a transient flag causes an explanatory stop.
- Every point is hybrid initialized independently, including the interior air.

Normalization uses diameter D = 0.229 m (rounded mesh extent), rho = 1.225 kg/m3
and mu = 1.7894e-5 Pa s. These are assumptions matching the inspected case;
if you change its material properties, update the corresponding script options.
The aerodynamic reference area is the nominal circular area pi D²/4, not a
numerically integrated projected area of the detailed ball geometry.

```
A   = pi D² / 4
Re  = rho U D / mu
S   = omega D / (2 U)       # signed equatorial speed / inlet speed
C_D = F_x / (0.5 rho U² A)
C_M = (-F_z) / (0.5 rho U² A)
```

Re and S are the two dimensionless independent coordinates; C_D or C_M is the
third coordinate. Spin ratio S depends on both U and omega, so the script
preserves the original sweep-grid connectivity rather than assuming uniform
spacing in Re and S.

## Outputs

Each run creates a unique folder under `football_results`:

- `results.csv`: all dimensional and dimensionless columns, updated each point.
- `tables.txt`: two printed tables, also displayed in the Python console.
- `configuration.json` and `wall_settings.json`: inputs and API settings.
- `C_D_3D.png`, `C_M_3D.png` and matching PDF figures: surface plots with sampled
  points overlaid. A one-point or one-dimensional sweep produces a scatter plot.
- Optional `point_*.cas.h5` and associated data files: per-point saved solutions.

All numeric output uses decimal dots. No example force values are fabricated.
The plotting stage runs only after a completed sweep; to replot completed rows
after an interrupted sweep:

```sh
python football_sweep.py --plot-only "football_results/run_.../results.csv"
```

These are steady iteration-endpoint forces. A transient study requires physical
time steps, a discarded startup interval, and time-averaged forces over a stated
sampling window; increasing steady iterations is not a substitute.

## Official references

- Wall momentum settings: https://fluent.docs.pyansys.com/version/stable/api/solver/_autosummary/settings/momentum_11.html
- Force-report example: https://systemcoupling.docs.pyansys.com/version/stable/examples/00-systemcoupling/turek_hron_fsi2.html
- Older report-return compatibility issue: https://github.com/ansys/pyfluent/issues/2919
