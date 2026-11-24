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

## Inspect cell adjacency

[`elbow_mesh_adjacency_graph.py`](../elbow%20cell%20adjacency%20graph/elbow_mesh_adjacency_graph.py)
loads a Fluent `.cas.h5` case through PyFluent, then reads its face
connectivity to build a cell adjacency graph. Cells are graph nodes; an edge
joins two cells that share an interior face. The script does not run iterations
or modify the case, but it does require a working Fluent installation and
license. See the [adjacency graph README](../elbow%20cell%20adjacency%20graph/README.md)
for dedicated usage instructions.

The script uses `ansys-fluent-core`, `h5py`, NumPy, and Matplotlib. Install
them into the Python environment used to run the script if needed:

```powershell
python -m pip install ansys-fluent-core h5py numpy matplotlib
```

Run it from the `elbow cell adjacency graph` directory. For example:

```powershell
python .\elbow_mesh_adjacency_graph.py `
  --case "C:\path\to\elbow_no_h2.msh.cas.h5" `
  --cell-id 5364 --rings 1 --show
```

`--cell-id` is a one-based Fluent cell ID (default: `1`). `--rings` may be `1`
or `2` (default: `1`). `--show` opens the saved graph in a Matplotlib window;
without it, the script only writes files. Use `--help` for the Fluent version,
processor-count, and output-directory options.

The default output folder is `adjacency_graph_results` in the current terminal
directory (the supplied graph preview and CSV are in that directory):

- `adjacency_graph.png` shows the full adjacency matrix and the selected cell's
  local node-link graph. The circular arrangement in the local view represents
  graph rings, not physical cell coordinates.
- `cell_adjacency.csv` lists interior faces with their face ID and the two
  adjacent one-based cell IDs. The plotted matrix collapses repeated pairs
  into one undirected edge.

The graph is read from the `.cas.h5` face-owner (`c0`) and face-neighbour
(`c1`) connectivity. Boundary faces have no second cell and are not graph
edges. The reader supports the single-section Fluent HDF5 layout used by the
supplied case.

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
