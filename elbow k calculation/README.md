# Elbow pressure-loss sweep

`elbow_k_Re_calculation.py` runs a steady Fluent calculation for each inlet
velocity supplied on the command line. It reads the case file, sets the
`hot-air-inlet` velocity to each sweep value, sets `h2-inlet` to zero, performs
hybrid initialization, and then runs the requested number of iterations.
The input case file is not overwritten.

## Results and calculations

For each velocity, the script records:

- Mass-weighted average total pressure at the swept inlet and outlet, and their
  difference.
- Mass imbalance as a percentage:
  `100 * abs(sum(signed mass flows)) / sum(abs(mass flows))`.
- Iterations for that velocity and cumulative iterations across the sweep.
- Reynolds number and pressure-loss constant `k`.

The supplied case defaults are:

- Swept inlet: `hot-air-inlet`
- Disabled inlet: `h2-inlet`
- Outlet: `outlet`
- Fluid cell zone: `fluid`
- Hydraulic diameter: `0.1016 m`

Density and laminar viscosity are volume-averaged over the fluid zone after
each velocity's iterations. The calculations are:

```text
Re = volume_average_density * inlet_velocity * hydraulic_diameter
     / volume_average_laminar_viscosity

k = (inlet_mass_average_total_pressure - outlet_mass_average_total_pressure)
    / (0.5 * volume_average_density * inlet_velocity^2)
```

The iteration counts include the requested solver iterations, not the internal
work performed by hybrid initialization. The script supports steady cases
only. Check convergence before using the reported pressure-loss results.

Each run creates a new timestamped folder under `elbow_k_results` in the
current terminal directory. It contains:

- `results.csv` and `tables.txt`
- `k_vs_Re.png` and `k_vs_Re.pdf`
- `configuration.json` with the case, settings, and calculation definitions

## Run from the VS Code PowerShell terminal

Change to this directory, set the Python executable that has both
`ansys-fluent-core` and `matplotlib` installed, then run the script. For
example:

```powershell
Set-Location 'G:\My Drive\_ nf wrk gdrive (nfgc)\cfse muiaic subject development on nfgc\pyton coding\git repo CfSE\CfSE\elbow k calculation'
$python = 'G:\My Drive\_ nf wrk gdrive (nfgc)\cfse muiaic subject development on nfgc\pyton coding\git repo CfSE\CfSE\magnus football\.venv\Scripts\python.exe'
$case = 'C:\Users\nfueyo\Desktop\cfseWorkDir\elbow\elbow_no_h2.msh.cas.h5'
& $python '.\elbow_k_Re_calculation.py' --case $case --velocities 0.2 0.4 0.6 --iterations 500
```

Velocities are in metres per second. `--iterations` is the number of steady
iterations per velocity; it defaults to 500. Use `--help` to see all options,
including processor count, Fluent version, output directory, boundary names,
fluid zone, and hydraulic diameter. The case and output paths can contain
spaces when quoted as shown.
