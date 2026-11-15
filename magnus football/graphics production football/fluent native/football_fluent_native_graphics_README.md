# Football plots rendered by Fluent

This version **drives Fluent's own graphics objects and PNG exporter**. It uses native Contours, Vectors, LICs, camera controls and Save Picture. When the Fluent/PyFluent settings API exposes native Annotation objects, it uses them for captions. If not, it still exports the plots without script captions and records that limitation in `manifest.json`. Python does not render, interpolate or modify the images. It requires only `ansys-fluent-core`, a compatible installed Fluent release and a Fluent licence.

## Run

```sh
python -m pip install ansys-fluent-core
python football_fluent_native_graphics.py --case "one.cas.h5" --data "one.dat.h5" --speed 25 --omega 50 --save-case --output "native_plots_U25_w50"
```

Replace the filenames, velocity and spin with the actual values for your selected solved case. `--speed` and `--omega` are annotation values, not commands to change the boundary conditions. Omega is in rad/s, about the supplied case's +y axis. The output folder must be new. Paths may contain spaces when quoted; the command works on Windows as one line.

Use `--version 25.2` to select a particular installed Fluent release. **Fluent 2025 R2 or newer is recommended**, with matching PyFluent. The settings API must expose native LIC objects; missing LIC support is reported as an error rather than replaced with Python graphics. Native text annotations are used when available; otherwise the PNGs are exported without the script captions, with a warning in the console and manifest. An old `.cas.h5` can be read into a supported newer Fluent release.

The default launch mode is **no GUI, with graphics enabled**. LIC requires a working Fluent graphics context on the execution machine. Do not change the launch mode to `no_gui_or_graphics`. Add `--gui` to see Fluent while it exports the plots or to diagnose a machine-specific graphics-driver issue. This is a local Fluent launch, so the case/data/output paths refer to the machine running Python.

## Nine separate PNGs

| Plane | Fluent graphics |
|---|---|
| Horizontal **y=0** | Static pressure, total pressure, dynamic pressure, velocity magnitude, velocity vectors, velocity LIC |
| Vertical **z=0** | Static pressure, velocity magnitude, velocity vectors |

These are **ball-centre** sections. The default y=0 is not the geometric mid-height of the original outer domain, whose y limits are asymmetric. For that alternative use `--horizontal-y -0.1755`. Set `--vertical-z` to change the vertical cut.

All four scalar fields are Fluent's **native fields**: `pressure`, `total-pressure`, `dynamic-pressure`, and `velocity-magnitude`. No Python formula is substituted for them. The native colour key provides the field name, units and range. For this incompressible, stationary-frame case, static pressure is gauge pressure and the pressure relationships are consistent with the reference pressure in Fluent. The native total-pressure definition remains Fluent's responsibility when applied to a different model.

The provided case calls the **exterior fluid cell zone `solid`**, and the interior ball volume `fluid`. The script therefore uses `--zone solid` by default and creates coordinate iso-surfaces restricted to that zone for contours and vectors. Change the option if the zone was renamed. In the tested Fluent 2026 R1 interface, native LIC accepts `plane_surface` objects rather than those zone-restricted iso-surfaces, so the script creates matching whole-domain planes for LIC. If another zone has velocity data, it may also appear in the LIC image.

## Framing, scales and annotations

- Physical orthographic views show the incoming +x flow **left to right** in both cuts. Horizontal: camera on the −y side, screen-up +z. Vertical: camera on the +z side, screen-up +y. The native axis triad remains visible.
- The requested view surrounds the ball/wake from **−2D to +8D**, with transverse bounds **−2D to +2D**, using D=0.229 m and centre (0,0,0). The camera expands that rectangle to match the image aspect and leave annotation margins; it does not stretch the geometry or show exactly cropped bounds. Example tighter view: `--xlim -1.5 4 --transverse-lim -1.25 1.25`.
- `--diameter` and `--centre X Y Z` adjust the reference geometry. Section coordinates are absolute metres; view bounds are in diameters relative to that centre.
- Output is **3200 × 1800 pixels**, configurable with `--width` and `--height`. A PNG signature and resolution check follows each Fluent export.
- Fluent computes field ranges over **both complete selected planes**. Each field then receives one custom range shared across planes, rather than changing its colour scale between pictures. These ranges include any off-camera part of those planes. Static pressure is symmetric around zero, and dynamic pressure/speed start at zero. Shared speed limits are used for speed contours, vectors and LIC. Ranges are not automatically shared between separate runs.
- The camera is reapplied after each object's `display()` so a saved display state cannot unexpectedly reset the requested orientation/zoom before export.
- When native annotation objects are available, Fluent annotations give the field/plane, case label, actual U/spin, ball diameter, flow direction, screen-up direction and one teaching caption. If unavailable, the script omits those captions and marks `native_annotations.enabled` as `false` in the manifest.
- Vectors use **in-plane projection**, variable lengths and Fluent's scale/skip controls. Colour represents full velocity magnitude. Start with the defaults; use `--vector-scale 0.015 --vector-skip 20` for shorter/sparser arrows. Native vector scaling depends on the cut's mesh, so the defaults cannot guarantee ideal density for every mesh. Fixed-length vectors are deliberately disabled so length retains magnitude information.
- LIC is Fluent's **native velocity LIC**, coloured by velocity magnitude when that release exposes the option. Default texture size is 16 (hundreds of pixels) and the optional native step control is 30. Adjust with `--lic-texture` and `--lic-steps`. Additional LIC display controls vary by release; unavailable optional controls retain Fluent's defaults. Symmetric LIC illustrates alignment; vectors are the clearer sign-of-motion reference.

## Reusable objects and records

Add `--save-case` to save **a new `graphics_setup.cas.h5` in the output folder**, containing the native surfaces, plots and annotations. The input case/data files are not overwritten. Reopen this new case with the original matching `.dat.h5` to inspect or edit the plots. No new solution data is calculated or written.

Each PNG also has a corresponding `*_settings.json`. `manifest.json` records the actual native object names, field ranges, camera state, exported image dimensions and completion/failure status. A native `graphics_diagnostic.json` is attempted if an error occurs. Temporary range-probe contours are retained in the saved case as well, with names identifying their purpose.

This script does not initialize, iterate, change physical boundary conditions or create synthetic results. A saved transient data file will yield its loaded snapshot, not a time average. A case/data pair alone does not prove convergence; check that separately.

## Validation status

Syntax, argument handling, camera handedness, equal-scale framing and PNG-header checks have been tested without Fluent. The live native graphics/API/export path has **not been executed here**, because no licensed Fluent renderer is available in this environment. It may require release-specific adjustments; failures are reported with native settings diagnostics rather than falling back to another renderer.

## Official interfaces used

- [Contour object](https://fluent.docs.pyansys.com/version/dev/api/solver/_autosummary/settings/contour_child.html)
- [Vector object](https://fluent.docs.pyansys.com/version/stable/api/solver/_autosummary/settings/vector_1_child.html)
- [Native LIC object](https://fluent.docs.pyansys.com/version/stable/api/solver/_autosummary/settings/lic_child.html)
- [Native annotations](https://fluent.docs.pyansys.com/version/dev/api/solver/_autosummary/settings/annotation_child.html)
- [Camera properties](https://fluent.docs.pyansys.com/version/stable/api/solver/_autosummary/settings/camera.html)
- [PNG hardcopy driver](https://fluent.docs.pyansys.com/version/stable/api/solver/tui/display/set/picture/driver/driver_contents.html)
