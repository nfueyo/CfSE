# Elbow mesh quality analysis

[`elbow_mesh_quality.py`](./elbow_mesh_quality.py) reads a Fluent `.cas.h5`
mesh and writes per-cell metrics, plots, and a PDF report for students. It
analyses the mesh geometry directly; it does not start Fluent or solve the CFD
case.

## Requirements

Use Python 3 with `h5py`, NumPy, and Matplotlib:

```powershell
python -m pip install h5py numpy matplotlib
```

No separate PDF package is needed; the report is written with Matplotlib.

## Run

From this directory, provide the path to the Fluent case:

```powershell
python .\elbow_mesh_quality.py "C:\cases\elbow.cas.h5"
```

By default, results are written to `.\elbow_mesh_results\`. To change the
output folder or the z-slices used in the centroid maps:

```powershell
python .\elbow_mesh_quality.py "C:\cases\elbow.cas.h5" `
  --output ".\results_run_01" `
  --slab-z-mm 0 15 `
  --slab-half-width-mm 2
```

Run `python .\elbow_mesh_quality.py --help` to see the available options.
The z-slab locations and half-width are in millimetres; the half-width must
be positive.

## Output files

Each run creates or updates these files in the output folder:

| File | Contents |
| --- | --- |
| `cell_metrics.csv` | One row per cell, with its ID, Fluent cell type, centroid coordinates, volume, and three quality indicators. Coordinates are in metres, volume in cubic metres, and characteristic size in millimetres. |
| `histograms.png` | Distributions of cell size, largest neighbour size ratio, and largest face angle, shown separately for hexahedra and polyhedra. Each cell type is normalized to its own percentage. |
| `centroid_slabs.png` | Metric values at cell centroids in the requested z-slices. |
| `worst_cells_3d.png` | The cells in the largest 1% of the computed face-angle metric, highlighted in the 3-D mesh view. |
| `mesh_quality_report.pdf` | An eight-page student report with definitions, formulas, data tables, the generated plots, interpretation guidance, and limitations. |

Running again with the same output folder replaces files with the same names.

## Metrics calculated

- **Characteristic cell size** is the cube root of the computed cell volume,
  `h = V^(1/3)`. It is a volume-based length scale, not necessarily an actual
  edge length.
- **Largest neighbour size ratio** is the largest symmetric ratio of cell
  sizes across any interior face, `max(h_i/h_j, h_j/h_i)`. A value of 1 means
  equal sizes; larger values indicate a sharper local size change.
- **Largest face angle** compares the line between two neighbouring cell
  centroids with the axis of their shared face normal. Zero degrees is aligned;
  larger angles indicate greater misalignment. Each cell receives the largest
  angle from its interior faces.

The neighbour ratio and angle are computed from **interior faces only** and
assigned to both cells sharing each face. The script recognizes Fluent cell
type code 4 as hexahedra and code 7 as polyhedra.

The reader expects the single-section Fluent HDF5 layout used by this case.
It checks face connectivity and rejects non-positive or invalid cell volumes;
other case layouts may need a different reader.

## Interpretation and limitations

These are geometric indicators calculated by this script. The face-angle
indicator is a custom centroid-to-face-normal measure and is not guaranteed to
match Fluent's native orthogonal quality or skewness. The analysis does not
report boundary-face quality, aspect ratio, wall `y+`, solver convergence, or
solution accuracy. No universal pass/fail thresholds are assigned: inspect
the distribution and locations, consult the quality guidance for the solver
and application, and use mesh-sensitivity studies to assess numerical impact.
