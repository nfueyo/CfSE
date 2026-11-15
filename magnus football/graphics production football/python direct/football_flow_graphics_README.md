# Football flow graphics

Produce **one annotated PNG per plot** from one already-solved Fluent case and its matching data file. This script does not run iterations or change the flow setup. It starts its own Fluent session, creates temporary sections, extracts their nodal fields and polygon connectivity, and closes that session without saving the case. The PNGs are rendered with Matplotlib rather than captured from Fluent's graphics window.

## Outputs

| Section | Images |
|---|---|
| Horizontal, **y = 0 m**, through the ball centre | Static pressure; total pressure; dynamic pressure; velocity magnitude; velocity vectors; velocity LIC |
| Vertical, **z = 0 m**, through the ball centre | Static pressure; velocity magnitude; velocity vectors |

The supplied CFD coordinates have **x downstream, y vertical, z transverse**. The horizontal plot shows x–z, and the vertical plot shows x–y. Incoming air is always left to right. Coordinates are plotted relative to the ball centre and divided by D, with equal geometric scales. The default view covers **−2D to +8D** downstream and **±2D** transversely: enough to show the ball and wake without wasting most of the image on the far field. Each 14 × 7.8 inch image is saved at 220 DPI by default (3080 × 1716 pixels).

These defaults mean **mid-height of the ball**, not the midpoint between the outer domain's unequal y bounds. To plot the geometric mid-height of the original outer domain, use `--horizontal-y -0.1755` instead. All sections are configurable.

## Install and run

Use a Python environment with access to an installed, licensed Fluent release and a compatible PyFluent version. Fluent 2024 R2 or later is recommended; loading an older saved case into a supported newer Fluent release is different from running an unsupported old Fluent executable. The field reader supports current request objects and the documented older scalar/surface access methods when available.

```sh
python -m pip install -r football_flow_graphics_requirements.txt
```

For example (replace both paths with the **matching solved files for one matrix case**, and specify that case's actual velocity and spin):

```sh
python football_flow_graphics.py --case "case_U25_w50.cas.h5" --data "case_U25_w50.dat.h5" --speed 25 --omega 50 --label "U25 / omega50" --output "plots_U25_w50"
```

The command works as a single line in Windows PowerShell as well as a Unix shell. `--speed` and `--omega` label the data; they **do not set the boundary conditions**. Do not put a different velocity/spin here from the one used to solve the loaded case. The spin annotation assumes the supplied case's +y axis. Speed also determines the vector-arrow reference scale. Omit `--omega` if the spin is unknown; the title then says so.

The script requires the `.dat.h5` solution, not only a case or the report `.out` files. Initialization fields are not a validated solution: use a case with convergence already assessed. For a transient solution the images represent the single loaded snapshot, not a time average.

The output directory must be **new**. This prevents mixing new images with previous exports. Use `--version 25.2` (for example) if you need to select a particular installed Fluent release, and `--processors 4` to choose the launch size.

## Important details for this football case

- The exterior fluid cell zone is misleadingly named **`solid`**. The zone named `fluid` is inside the ball. The default `--zone solid` intentionally extracts only external air. The script checks that a centre-plane cut has a hole at the ball centre, rather than silently plotting the internal ball cells.
- The density is checked against material **`air`**. The default `--rho 1.225` must match its constant density. If you have renamed the material, use `--material NAME`. This renderer's pressure definitions require constant-density flow.
- Defaults: D = **0.229 m**, ball centre = **(0,0,0) m**, static gauge-pressure key = **`pressure`**. Change these if you changed the case. `--pressure-field` must identify a **gauge static-pressure** field, not a coefficient or absolute-pressure field.
- Static pressure comes from Fluent. Dynamic pressure is calculated as **q = 0.5 ρ (Vx² + Vy² + Vz²)**. Total pressure is calculated consistently as **p_total,gauge = p_static,gauge + q** in the stationary inertial frame. These are not additional fields read from Fluent or a compressible stagnation-pressure formula. Absolute total pressure would require adding the operating pressure.
- The density check does not establish that the selected material is assigned to a renamed/generalised cell zone. If you adapt the script to a different model, confirm its material assignment as well as the zone selection.
- Shared colour limits are based on the actual minima and maxima in both visible sections, with **no percentile clipping**. Static pressure has a symmetric scale about zero gauge pressure; speed and dynamic pressure start at zero. Colours are comparable between these two planes of this run, not automatically between different runs.
- The original section connectivity is retained when interpolating to the regular image grid. The script does not create a Delaunay interpolation that bridges the ball hole. Overlapping/invalid cut topology raises a diagnostic error rather than quietly inventing field values.
- Vectors show the **in-plane components**, with equal physical arrow scaling and an annotated velocity key. Colours show full three-component speed. A 2D section therefore does not show the entire 3D flow direction.
- The LIC is a genuine, reproducible **symmetric line integral convolution**: seeded noise is integrated forwards and backwards along the projected velocity with a Hann kernel and midpoint integration. It is **not a streamline or pathline plot**, and does not indicate the sign of travel. Its hue represents full speed; brightness is modulated by the texture. It stops at section holes and domain edges. `--lic-length` controls texture length in pixels and `--pixels` controls sampling resolution, not CFD mesh resolution.
- Rotating walls can exchange mechanical energy with the fluid. The total-pressure caption therefore asks students to inspect its distribution rather than assuming every change is solely an irreversible loss.
- PNGs have explanatory captions, flow direction, coordinate labels, field units, case parameters, and a ball-section callout. The callout targets the projected ball centre, not a CFD separation point.

## View and resolution examples

Tighter view to emphasise the boundary layer and immediate wake:

```sh
python football_flow_graphics.py --case "one.cas.h5" --data "one.dat.h5" --speed 25 --omega 50 --xlim -1.5 4 --transverse-lim -1.25 1.25 --pixels 1200 --dpi 250 --output "plots_close"
```

Select a plane or change geometry inputs with `--horizontal-y`, `--vertical-z`, `--centre X Y Z`, and `--diameter`. The plotted axes are always measured from that configured centre. Section coordinates are **absolute metres**, while the view limits are in **diameters from the centre**.

## Replot without Fluent

Each export saves `horizontal_section.npz`, `vertical_section.npz`, and a `manifest.json` describing fields, coordinates, colour limits, parameters and output filenames. The NPZ files contain the extracted section geometry, pressure and three velocity components, and load without Python pickle. `available_scalar_fields.json` helps diagnose a release-specific field-name mismatch.

```sh
python football_flow_graphics.py --from-cache "plots_U25_w50" --speed 25 --pixels 1200 --dpi 250 --output "plots_U25_w50_fine"
```

For cache replotting, the original physical parameters and case label are restored from the original manifest. View limits, PNG resolution and LIC length remain configurable. The required `--speed` argument is then superseded by the cached value. Keep the manifest alongside both NPZ files.

## Validation and limits

This delivery was checked without a licensed Fluent session: syntax and command-line parsing, all nine image exports, field-cache round trips, retention of the ball hole, duplicate-partition welding, full-3D pressure calculations, and LIC alignment/reversal/masking were verified using **synthetic analytic data**. The synthetic images are conspicuously watermarked and are not football CFD results. The live Fluent extraction step has **not** been executed here; release-specific API issues may still need adjustment. The script writes a failed-status manifest with the exception if extraction or rendering fails.

Run the numerical checks with:

```sh
python check_football_flow_graphics.py
```

Optional renderer demonstration:

```sh
python football_flow_graphics.py --demo --speed 25 --pixels 500 --output "renderer_demo"
```

The demo is an analytic flow around a circular hole, not a spinning-football prediction. Its purpose is to check geometry, masks, scales and image layout.

## API references

- [PyFluent field-data requests and response models](https://fluent.docs.pyansys.com/version/stable/api/fields/field_data/requests.html)
- [PyFluent live field-data access](https://fluent.docs.pyansys.com/version/stable/api/fields/field_data/live_field_data.html)
- [Iso-surface settings, including cell-zone restriction](https://fluent.docs.pyansys.com/version/stable/api/solver/_autosummary/settings/iso_surface_child.html)
- [Fluent's explanation of LIC](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/flu_ug_sec_graphics_basic.html)
