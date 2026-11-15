"""Make nine PNGs with Fluent's own graphics renderer; annotate when supported.

Requires installed/licensed Fluent with native Contour, Vector and LIC settings
objects (2025 R2 or newer recommended) and compatible PyFluent. Native Annotation
settings are optional; without them, exported images omit the script captions.
No Matplotlib/PyVista, external rendering, interpolation or Python-computed LIC.
No iterations or initialization. See football_fluent_native_graphics_README.md.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import json
import math
from pathlib import Path
import struct
import uuid


FIELDS = {'pressure': 'pressure', 'total': 'total-pressure',
          'dynamic': 'dynamic-pressure', 'speed': 'velocity-magnitude'}
CAPTIONS = {
    'pressure': ('Static pressure', 'Pressure differences around the ball contribute to aerodynamic force.'),
    'total': ('Total pressure', 'Inspect the wake and stagnation-pressure distribution; a rotating wall can exchange energy.'),
    'dynamic': ('Dynamic pressure', 'Local kinetic energy per unit volume: q = 0.5 rho |V|^2.'),
    'speed': ('Velocity magnitude', 'Inspect acceleration beside the ball and the downstream speed deficit.'),
    'vectors': ('Velocity vectors', 'Arrows: in-plane velocity, variable length. Colour: full three-component speed.'),
    'lic': ('Velocity LIC', 'Native Fluent line integral convolution of velocity on the horizontal plane.'),
}


def names(group):
    return set(getattr(group, 'child_names', ()))


def child(group, *candidates):
    for name in candidates:
        if name in names(group):
            return getattr(group, name)
    raise RuntimeError(f'Expected {candidates} in {getattr(group, "path", type(group).__name__)}; '
                       f'available children: {sorted(names(group))}')


def annotation_collection(graphics):
    for name in ('annotation', 'annotations'):
        if name in names(graphics):
            return getattr(graphics, name)
    return None


def set_value(group, name, value):
    child(group, name).set_state(value)


def optional(group, name, value):
    if name in names(group):
        setting = getattr(group, name)
        if setting.is_active():
            setting.set_state(value)
            return True
    return False


def choice(setting, *wanted):
    """Use this release's actual spelling; fail rather than guess an enum."""
    allowed = setting.allowed_values()
    if allowed is None:
        setting.set_state(wanted[0])
        return wanted[0]
    normalize = lambda s: str(s).lower().replace('_', '-').replace(' ', '-')
    for desired in wanted:
        for actual in allowed:
            if normalize(actual) == normalize(desired):
                setting.set_state(actual)
                return actual
    raise RuntimeError(f'Cannot select {wanted}; allowed values: {allowed}')


def select_range_mode(obj, *, automatic):
    setting = child(obj, 'range_option')
    if 'option' in names(setting):
        mode = 'auto-range-on' if automatic else 'auto-range-off'
        actual = choice(child(setting, 'option'), mode)
        if automatic and 'auto_range_on' in names(setting):
            optional(child(setting, 'auto_range_on'), 'global_range', False)
        return actual
    return choice(setting, 'local' if automatic else 'custom')


def new_object(collection, name):
    collection.create(name=name)
    return collection[name]


def camera_spec(args, plane):
    """Physical orthographic camera; right is +x in BOTH views.

    Horizontal: look from -y, screen-up +z. Vertical: look from +z, up +y.
    The right vector is (target-position) cross up, so neither view mirrors x.
    Width/height preserve geometry; margins leave space for native annotations.
    """
    D = args.diameter
    xmid = args.centre[0] + .5*sum(args.xlim)*D
    width = (args.xlim[1]-args.xlim[0])*D
    height = (args.transverse_lim[1]-args.transverse_lim[0])*D
    transverse = 2 if plane == 'horizontal' else 1
    target = list(args.centre)
    target[0] = xmid
    target[transverse] += .5*sum(args.transverse_lim)*D
    target[1 if plane == 'horizontal' else 2] = args.horizontal_y if plane == 'horizontal' else args.vertical_z
    position = target.copy()
    position[1 if plane == 'horizontal' else 2] += (-1 if plane == 'horizontal' else 1)*20*D
    # Fit the requested rectangle to the image aspect without distorting it.
    # This can reveal extra space beyond the requested rectangle, never crop it.
    aspect = args.width/args.height
    field_width = max(width, height*aspect)*1.32
    return dict(position=position, target=target,
                up_vector=[0., 0., 1.] if plane == 'horizontal' else [0., 1., 0.],
                field=[field_width, field_width/aspect])


def set_camera(graphics, args, plane):
    cam = child(child(graphics, 'views'), 'camera')
    choice(child(cam, 'projection'), 'orthographic', 'parallel')
    spec = camera_spec(args, plane)
    for name, value in spec.items():
        set_value(cam, name, value)
    return spec


def configure_picture(solver, graphics, args):
    picture = child(graphics, 'picture')
    optional(picture, 'use_window_resolution', False)
    set_value(picture, 'x_resolution', args.width)
    set_value(picture, 'y_resolution', args.height)
    # TUI driver selection is explicit: do not rely on a .png suffix alone.
    solver.tui.display.set.picture.driver.png()
    optional(picture, 'landscape', True)
    optional(picture, 'invert_background', True)
    # Fluent retains its native colour key, captions, units and axis triad.
    return picture


def png_size(path):
    with path.open('rb') as f:
        header = f.read(24)
    if len(header) != 24 or header[:8] != b'\x89PNG\r\n\x1a\n':
        raise RuntimeError(f'Fluent did not write a valid PNG: {path}')
    return list(struct.unpack('>II', header[16:24]))


def annotation(collection, name, text, position, size=20):
    obj = new_object(collection, name)
    set_value(obj, 'text', text)
    set_value(obj, 'position', position)
    set_value(obj, 'font_size', size)
    if 'font_color' in names(obj):
        choice(child(obj, 'font_color'), 'black')
    return name


def attach_annotations(collection, obj, prefix, plane, plot, args):
    if collection is None:
        return False
    coordinate = args.horizontal_y if plane == 'horizontal' else args.vertical_z
    normal = 'y' if plane == 'horizontal' else 'z'
    up = '+z' if plane == 'horizontal' else '+y'
    title, description = CAPTIONS[plot]
    label = args.label or args.case.stem
    spin = '' if args.omega is None else f' | spin = {args.omega:g} rad/s about +y'
    lines = [
        (f'{title} | {plane.capitalize()} section: {normal} = {coordinate:g} m', [.17, .95], 28),
        (f'{label} | U_in = {args.speed:g} m/s{spin} | D = {args.diameter:g} m', [.17, .90], 20),
        (description, [.17, .12], 20),
        (f'Incoming air +x: left to right | Screen up: {up} | External air only | Orthographic view', [.17, .07], 18),
    ]
    annotations = [annotation(collection, f'{prefix}-a{i}', text, position, size)
                   for i, (text, position, size) in enumerate(lines)]
    set_value(obj, 'annotations_list', annotations)
    return True


def set_field(obj, field, surfaces):
    prop = child(obj, 'field')
    allowed = prop.allowed_values()
    if allowed is not None and field not in allowed:
        raise RuntimeError(f'Native field {field!r} unavailable. Allowed fields: {allowed}')
    prop.set_state(field)
    surface_setting = child(obj, 'surfaces_list')
    allowed = surface_setting.allowed_values()
    if allowed is not None and any(s not in allowed for s in surfaces):
        raise RuntimeError(f'This native graphics object does not accept surface(s) {surfaces}. '
                           f'Allowed surfaces: {allowed}.')
    surface_setting.set_state(surfaces)
    optional(obj, 'draw_mesh', False)


def range_values(obj):
    """Compute with Fluent on selected surfaces, including both planes."""
    select_range_mode(obj, automatic=True)
    opts = child(obj, 'range_options')
    optional(opts, 'global_range', False)
    optional(opts, 'auto_range', True)
    opts.compute()
    lo = float(child(opts, 'minimum').get_state())
    hi = float(child(opts, 'maximum').get_state())
    if not math.isfinite(lo) or not math.isfinite(hi) or hi < lo:
        raise RuntimeError(f'Invalid native field range {lo}, {hi}; check loaded solution data.')
    return lo, hi


def custom_range(obj, limits):
    select_range_mode(obj, automatic=False)
    opts = child(obj, 'range_options')
    optional(opts, 'auto_range', False)
    # Some generated schemas put custom limits in an auto_range_off subgroup.
    if 'auto_range_off' in names(opts):
        opts = child(opts, 'auto_range_off')
    set_value(opts, 'minimum', limits[0])
    set_value(opts, 'maximum', limits[1])


def vector_options(obj, args):
    choice(child(obj, 'vector_field'), 'velocity')
    vopt = child(obj, 'vector_opt')
    set_value(vopt, 'in_plane', True)
    set_value(vopt, 'fixed_length', False)
    optional(vopt, 'scale_head', .15)
    opts = child(obj, 'options')
    optional(opts, 'auto_scale', False)
    if 'scale' in names(opts):
        set_value(opts, 'scale', args.vector_scale)
        set_value(opts, 'skip', args.vector_skip)
    else:
        set_value(obj, 'scale', args.vector_scale)
        set_value(obj, 'skip', args.vector_skip)
    # "arrow" is a common value; retain the release's existing style otherwise.
    for group, prop in ((opts, 'vector_style'), (obj, 'style')):
        if prop in names(group):
            setting = child(group, prop)
            allowed = setting.allowed_values()
            for value in allowed or ():
                if str(value).lower() in ('arrow', 'arrows'):
                    setting.set_state(value)
                    return


def lic_options(obj, args):
    choice(child(obj, 'vector_field'), 'velocity')
    optional(obj, 'lic_color_by_field', True)
    optional(obj, 'lic_normalize', True)
    optional(obj, 'lic_oriented', False)
    optional(obj, 'texture_size', args.lic_texture)
    optional(obj, 'lic_max_steps', args.lic_steps)


def run(solver, args, folder, manifest):
    s = solver.settings
    s.file.read_case(file_name=str(args.case.resolve()))
    s.file.read_data(file_name=str(args.data.resolve()))
    manifest['fluent_version'] = str(solver.get_fluent_version())
    graphics = s.results.graphics
    # Require real native objects. There is deliberately no external-renderer fallback.
    for collection in ('contour', 'vector', 'lic'):
        child(graphics, collection)
    annotation_group = annotation_collection(graphics)
    manifest['native_annotations'] = {'enabled': annotation_group is not None}
    if annotation_group is None:
        manifest['native_annotations']['reason'] = (
            'This Fluent/PyFluent version does not expose native annotation objects; '
            'PNG exports omit the script captions.')
        print('WARNING: Native annotation objects are unavailable; exporting PNGs without script captions.',
              flush=True)
    manifest['graphics_children'] = sorted(names(graphics))
    zones = s.setup.cell_zone_conditions.fluid.get_object_names()
    if args.zone not in zones:
        raise ValueError(f'External fluid zone {args.zone!r} missing. Available fluid zones: {zones}')
    prefix = 'ffl-' + uuid.uuid4().hex[:8]
    surfaces = {}
    for plane, axis, coordinate in [('horizontal', 'y', args.horizontal_y),
                                     ('vertical', 'z', args.vertical_z)]:
        name = f'{prefix}-{plane}'
        obj = new_object(s.results.surfaces.iso_surface, name)
        obj.field = axis+'-coordinate'
        obj.zones = [args.zone]
        obj.iso_values = [coordinate]
        surfaces[plane] = name
    manifest['surfaces'] = surfaces
    lic_surfaces = {}
    for plane, method, axis, coordinate in [
        ('horizontal', 'zx-plane', 'y', args.horizontal_y),
        ('vertical', 'xy-plane', 'z', args.vertical_z),
    ]:
        name = f'{prefix}-{plane}-lic'
        obj = new_object(s.results.surfaces.plane_surface, name)
        set_value(obj, 'method', method)
        set_value(obj, axis, coordinate)
        lic_surfaces[plane] = name
    manifest['native_lic_surfaces'] = lic_surfaces
    native_contours = child(graphics, 'contour')
    ranges = {}
    for plot, field in FIELDS.items():
        probe = new_object(native_contours, f'{prefix}-range-{plot}')
        set_field(probe, field, list(surfaces.values()))
        lo, hi = range_values(probe)
        if plot == 'pressure':
            bound = max(abs(lo), abs(hi), 1e-8)
            lo, hi = -bound, bound
        elif plot in ('dynamic', 'speed'):
            lo = 0.
        if hi <= lo:
            hi = lo + max(abs(lo)*1e-6, 1e-6)
        ranges[plot] = [lo, hi]
    manifest['shared_native_ranges'] = ranges
    picture = configure_picture(solver, graphics, args)
    manifest['plots'] = []
    for plane, surface in surfaces.items():
        plots = ('pressure', 'total', 'dynamic', 'speed', 'vectors', 'lic') if plane == 'horizontal' else ('pressure', 'speed', 'vectors')
        for plot in plots:
            name = f'{prefix}-{plane}-{plot}'
            collection = child(graphics, 'vector' if plot == 'vectors' else 'lic' if plot == 'lic' else 'contour')
            obj = new_object(collection, name)
            field_plot = 'speed' if plot in ('vectors', 'lic') else plot
            plot_surface = lic_surfaces[plane] if plot == 'lic' else surface
            set_field(obj, FIELDS[field_plot], [plot_surface])
            if plot == 'vectors':
                vector_options(obj, args)
            elif plot == 'lic':
                lic_options(obj, args)
            else:
                set_value(obj, 'filled', True)
                set_value(obj, 'node_values', True)
                optional(obj, 'contour_lines', False)
            custom_range(obj, ranges[field_plot])
            attach_annotations(annotation_group, obj, name, plane, plot, args)
            # Display first: a graphics object's saved display state can reset
            # the camera. Apply our camera AFTER display, BEFORE save-picture.
            obj.display()
            view = set_camera(graphics, args, plane)
            path = folder / f'{plane}_{plot}.png'
            picture.save_picture(file_name=str(path))
            resolution = png_size(path)
            if resolution != [args.width, args.height]:
                raise RuntimeError(f'Unexpected PNG resolution {resolution}; requested {args.width}x{args.height}.')
            state = obj.get_state()
            (folder / f'{plane}_{plot}_settings.json').write_text(
                json.dumps(state, indent=2, default=str), encoding='utf-8')
            manifest['plots'].append(dict(file=path.name, native_object=name,
                                          field=FIELDS[field_plot], resolution=resolution, camera=view))
            print(path, flush=True)
    if args.save_case:
        saved = folder / 'graphics_setup.cas.h5'
        s.file.write_case(file_name=str(saved))
        manifest['graphics_case'] = str(saved)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--case', type=Path, required=True)
    p.add_argument('--data', type=Path, required=True, help='Matching solved data file')
    p.add_argument('--speed', type=float, required=True, help='Actual solved inlet speed, m/s (label only)')
    p.add_argument('--omega', type=float, help='Actual solved spin, rad/s about +y (label only)')
    p.add_argument('--label', help='Short label for annotations')
    p.add_argument('--zone', default='solid', help='EXTERIOR fluid cell zone; supplied case calls it solid')
    p.add_argument('--diameter', type=float, default=.229)
    p.add_argument('--centre', type=float, nargs=3, default=(0., 0., 0.))
    p.add_argument('--horizontal-y', type=float, default=0.)
    p.add_argument('--vertical-z', type=float, default=0.)
    p.add_argument('--xlim', type=float, nargs=2, default=(-2., 8.), help='Requested streamwise view bounds in D')
    p.add_argument('--transverse-lim', type=float, nargs=2, default=(-2., 2.), help='Requested transverse bounds in D')
    p.add_argument('--width', type=int, default=3200)
    p.add_argument('--height', type=int, default=1800)
    p.add_argument('--vector-scale', type=float, default=.025, help='Fluent native vector scale factor; tune if needed')
    p.add_argument('--vector-skip', type=int, default=12, help='Fluent native vector skip factor')
    p.add_argument('--lic-texture', type=int, default=16, help='Fluent texture size (hundreds of pixels)')
    p.add_argument('--lic-steps', type=int, default=30)
    p.add_argument('--processors', type=int, default=4)
    p.add_argument('--version', help='Installed Fluent release, e.g. 25.2')
    p.add_argument('--gui', action='store_true', help='Show Fluent GUI while exporting')
    p.add_argument('--save-case', action='store_true', help='Save a NEW case copy with reusable native graphics objects')
    p.add_argument('--output', type=Path, default=Path('football_native_graphics'))
    return p


def main(argv=None):
    p = parser(); args = p.parse_args(argv)
    if not args.case.is_file() or not args.data.is_file():
        p.error('Case and matching solved data files must exist.')
    if args.speed <= 0 or args.diameter <= 0 or args.vector_scale <= 0:
        p.error('Speed, diameter and vector scale must be positive.')
    if args.xlim[0] >= args.xlim[1] or args.transverse_lim[0] >= args.transverse_lim[1]:
        p.error('View bounds must increase.')
    if args.width < 640 or args.height < 480 or args.processors < 1 or args.vector_skip < 0:
        p.error('Use width >=640, height >=480, processors >=1 and skip >=0.')
    if not 5 <= args.lic_steps <= 50 or args.lic_texture < 1:
        p.error('LIC steps must be 5–50 and texture size positive.')
    folder = args.output.resolve(); folder.mkdir(parents=True, exist_ok=False)
    import ansys.fluent.core as pyfluent
    manifest = {'created': datetime.now().isoformat(), 'configuration': vars(args),
                'renderer': 'Fluent native graphics only', 'pyfluent_version': pyfluent.__version__}
    solver = None
    try:
        kw = dict(mode='solver', dimension=3, precision='double',
                  processor_count=args.processors, cwd=str(folder),
                  ui_mode='gui' if args.gui else 'no_gui')
        if args.version:
            kw['product_version'] = args.version
        # no_gui retains graphics. Do not use no_gui_or_graphics for LIC/PNG.
        solver = pyfluent.launch_fluent(**kw)
        run(solver, args, folder, manifest)
        manifest['status'] = 'complete'
    except Exception as exc:
        manifest['status'] = 'failed'; manifest['error'] = repr(exc)
        if solver is not None:
            try:
                (folder / 'graphics_diagnostic.json').write_text(json.dumps(
                    solver.settings.results.graphics.get_state(), indent=2, default=str))
            except Exception as diagnostic_error:
                manifest['diagnostic_error'] = repr(diagnostic_error)
        raise
    finally:
        (folder / 'manifest.json').write_text(json.dumps(manifest, indent=2, default=str), encoding='utf-8')
        if solver is not None:
            solver.exit()
    print(f'\nSaved {len(manifest["plots"])} Fluent-rendered PNGs to {folder}')


if __name__ == '__main__':
    main()
