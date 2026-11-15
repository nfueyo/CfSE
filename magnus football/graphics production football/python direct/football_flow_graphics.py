"""Annotated plane plots from one solved Fluent football case.

PyFluent extracts nodal fields and the original section connectivity. Matplotlib
renders PNGs; LIC is an actual noise convolution along the in-plane velocity.
No solution iterations, initialization, boundary edits, or case writes occur.
See football_flow_graphics_README.md. --help needs no third-party imports.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
import importlib
import json
from pathlib import Path
import sys
import uuid


@dataclass
class Section:
    name: str
    axis: int
    coordinate: float
    vertices: object
    faces: list
    pressure: object
    velocity: object


def load_libraries():
    global np, plt, mtri, Normalize, TwoSlopeNorm, LineCollection
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from matplotlib.colors import Normalize, TwoSlopeNorm
    from matplotlib.collections import LineCollection
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})


def unwrap(value):
    """Old PyFluent containers have .data; current responses are arrays."""
    if isinstance(value, np.ndarray):
        return value
    if hasattr(value, 'data') and not isinstance(value, (dict, list, tuple)):
        return value.data
    return value


def one_surface(response, surface_id, name):
    if not isinstance(response, dict):
        return response
    for key in (surface_id, str(surface_id), name):
        if key in response:
            return response[key]
    if len(response) == 1:
        return next(iter(response.values()))
    raise RuntimeError(f'Cannot identify {name}: response keys {list(response)}')


def member(obj, *names):
    for name in names:
        if isinstance(obj, dict) and name in obj:
            return unwrap(obj[name])
        if hasattr(obj, name):
            return unwrap(getattr(obj, name))
    if isinstance(obj, dict):
        for key, value in obj.items():
            if getattr(key, 'value', None) in names:
                return unwrap(value)
    raise RuntimeError(f'Missing geometry member {names}; type={type(obj)}')


def decode_faces(raw):
    """Accept documented per-face arrays or count-prefixed flat connectivity."""
    if isinstance(raw, (list, tuple)) and len(raw) and np.ndim(raw[0]) > 0:
        return [np.asarray(f, dtype=np.int64).ravel() for f in raw]
    if isinstance(raw, np.ndarray) and raw.ndim == 2:
        return [np.asarray(f, dtype=np.int64) for f in raw]
    flat = np.asarray(unwrap(raw), dtype=np.int64).ravel()
    faces, pos = [], 0
    while pos < len(flat):
        n = int(flat[pos]); pos += 1
        if n < 3 or pos + n > len(flat):
            raise ValueError('Invalid count-prefixed face connectivity.')
        faces.append(flat[pos:pos+n]); pos += n
    return faces


class FieldReader:
    """Modern request API, with the documented legacy API for older releases."""
    def __init__(self, solver):
        self.fd = solver.fields.field_data
        self.info = solver.fields.field_info
        self.requests = None
        for module in ('ansys.fluent.core.fields.field_data',
                       'ansys.fluent.core.fields.field_data_interfaces',
                       'ansys.fluent.core.field_data_interfaces'):
            try:
                mod = importlib.import_module(module)
                if hasattr(mod, 'ScalarFieldDataRequest'):
                    self.requests = mod
                    break
            except ImportError:
                pass
        if not hasattr(self.fd, 'get_field_data'):
            self.requests = None

    def scalar(self, field, sid, name):
        if self.requests:
            request = self.requests.ScalarFieldDataRequest(
                field_name=field, surfaces=[sid], node_value=True, boundary_value=True)
            response = self.fd.get_field_data(request)
        else:
            response = self.fd.get_scalar_field_data(
                field_name=field, surfaces=[sid], node_value=True, boundary_value=True)
        values = np.asarray(unwrap(one_surface(response, sid, name)), dtype=float).ravel()
        if not values.size or not np.isfinite(values).all():
            raise RuntimeError(f'{name}: empty/non-finite {field}; load solved data.')
        return values

    def geometry(self, sid, name):
        if self.requests:
            request = self.requests.SurfaceFieldDataRequest(
                surfaces=[sid], data_types=['vertices', 'faces'],
                flatten_connectivity=False)
            response = self.fd.get_field_data(request)
        else:
            response = self.fd.get_surface_data(
                surfaces=[sid], data_types=['vertices', 'faces'])
        geometry = one_surface(response, sid, name)
        vertices = np.asarray(member(geometry, 'vertices'), dtype=float).reshape(-1, 3)
        faces = decode_faces(member(geometry, 'connectivity', 'faces'))
        if not vertices.size or not faces:
            raise RuntimeError(f'Empty section {name}. Check its coordinate and cell zone.')
        if min(f.min() for f in faces) < 0 or max(f.max() for f in faces) >= len(vertices):
            raise ValueError('Connectivity is not zero-based vertex indexing.')
        return vertices, faces


def extract(args, folder):
    import ansys.fluent.core as pyfluent
    kwargs = dict(mode='solver', dimension=3, precision='double',
                  processor_count=args.processors, ui_mode='no_gui', cwd=str(folder))
    if args.version:
        kwargs['product_version'] = args.version
    solver = pyfluent.launch_fluent(**kwargs)
    manifest = {'pyfluent_version': pyfluent.__version__, 'fields': {}}
    try:
        settings = solver.settings
        settings.file.read_case(file_name=str(args.case.resolve()))
        settings.file.read_data(file_name=str(args.data.resolve()))
        zones = settings.setup.cell_zone_conditions.fluid.get_object_names()
        if args.zone not in zones:
            raise ValueError(f'External fluid zone {args.zone!r} not found. Fluid zones: {zones}')
        density = settings.setup.materials.fluid[args.material].density.get_state()
        manifest['material_density_setting'] = density
        if density.get('option') != 'constant':
            raise ValueError('These pressure definitions require constant-density air. '
                             f'Material density setting: {density}')
        if not np.isclose(float(density['value']), args.rho, rtol=1e-6):
            raise ValueError(f'--rho {args.rho} differs from material density {density["value"]}.')
        reader = FieldReader(solver)
        available = reader.info.get_scalar_fields_info()
        (folder / 'available_scalar_fields.json').write_text(
            json.dumps(available, indent=2, default=str), encoding='utf-8')
        needed = [args.pressure_field, 'x-velocity', 'y-velocity', 'z-velocity']
        absent = [f for f in needed if f not in available]
        if absent:
            raise ValueError(f'Missing fields {absent}. See available_scalar_fields.json.')
        surfaces = settings.results.surfaces.iso_surface
        sections = []
        specs = [('horizontal', 1, args.horizontal_y), ('vertical', 2, args.vertical_z)]
        for label, axis, coordinate in specs:
            name = 'ffl-' + label + '-' + uuid.uuid4().hex[:8]
            surfaces.create(name=name)
            surface = surfaces[name]
            # Restrict the section to the EXTERIOR air. In the supplied case its
            # misleading saved name is "solid"; "fluid" is the ball interior.
            surface.field = ['x-coordinate', 'y-coordinate', 'z-coordinate'][axis]
            surface.zones = [args.zone]
            surface.iso_values = [coordinate]
            metadata = reader.info.get_surfaces_info()
            if name not in metadata:
                raise RuntimeError(f'Section {name} was not created: {surface.get_state()}')
            ids = metadata[name]['surface_id']
            if isinstance(ids, int):
                ids = [ids]
            all_vertices, all_faces, all_pressure, all_velocity = [], [], [], []
            offset = 0
            for sid in ids:
                xyz, faces = reader.geometry(sid, name)
                p = reader.scalar(args.pressure_field, sid, name)
                velocity = np.column_stack([
                    reader.scalar(field, sid, name)
                    for field in ('x-velocity', 'y-velocity', 'z-velocity')])
                if len(p) != len(xyz) or len(velocity) != len(xyz):
                    raise RuntimeError('Expected nodal scalar data aligned with section vertices.')
                all_vertices.append(xyz); all_pressure.append(p); all_velocity.append(velocity)
                all_faces.extend(f + offset for f in faces)
                offset += len(xyz)
            section = Section(label, axis, coordinate, np.vstack(all_vertices),
                              all_faces, np.concatenate(all_pressure), np.vstack(all_velocity))
            if np.max(np.abs(section.vertices[:, axis] - coordinate)) > args.diameter * 1e-5:
                raise RuntimeError('Section vertices are not on the requested coordinate plane.')
            sections.append(section)
            manifest['fields'][label] = {'surface': name, 'axis': axis,
                'coordinate_m': coordinate, 'vertices': len(section.vertices),
                'polygons': len(section.faces), 'external_cell_zone': args.zone}
        return sections, manifest
    finally:
        # Closing our unsaved, script-owned session also discards temporary cuts.
        # No pre-existing user session is connected to or closed.
        solver.exit()


def save_section(section, path):
    flat = np.concatenate([np.r_[len(f), f] for f in section.faces])
    np.savez_compressed(path, name=section.name, axis=section.axis,
                        coordinate=section.coordinate, vertices=section.vertices,
                        connectivity=flat, pressure=section.pressure, velocity=section.velocity)


def load_section(path):
    with np.load(path, allow_pickle=False) as d:
        return Section(str(d['name']), int(d['axis']), float(d['coordinate']),
                       d['vertices'], decode_faces(d['connectivity']),
                       d['pressure'], d['velocity'])


def section_grid(section, args):
    """Use the section's own polygons, never a Delaunay mesh across the ball hole."""
    transverse_axis = 2 if section.axis == 1 else 1
    xyz = np.asarray(section.vertices)
    xy = (xyz[:, [0, transverse_axis]] -
          np.array(args.centre)[[0, transverse_axis]]) / args.diameter
    # Welding shared vertices avoids overlapping duplicate edges in parallel data.
    _, inverse = np.unique(np.round(xy, 10), axis=0, return_inverse=True)
    counts = np.bincount(inverse)
    def weld(a):
        return np.bincount(inverse, weights=a) / counts
    points = np.column_stack([weld(xy[:, k]) for k in range(2)])
    p = weld(np.asarray(section.pressure))
    vel = np.column_stack([weld(np.asarray(section.velocity)[:, k]) for k in range(3)])
    triangles = []
    for face in section.faces:
        f = inverse[face]
        f = f[np.r_[True, np.diff(f) != 0]]
        if len(f) > 1 and f[0] == f[-1]:
            f = f[:-1]
        for j in range(1, len(f)-1):
            triangles.append([f[0], f[j], f[j+1]])
    triangles = np.asarray(triangles, dtype=int).reshape(-1, 3)
    a, b, c = [points[triangles[:, k]] for k in range(3)]
    area2 = (b[:, 0]-a[:, 0])*(c[:, 1]-a[:, 1]) - (b[:, 1]-a[:, 1])*(c[:, 0]-a[:, 0])
    triangles = triangles[np.abs(area2) > 1e-16]
    negative = area2[np.abs(area2) > 1e-16] < 0
    triangles[negative] = triangles[negative][:, [0, 2, 1]]
    _, unique = np.unique(np.sort(triangles, axis=1), axis=0, return_index=True)
    triangles = triangles[np.sort(unique)]
    if not len(triangles):
        raise RuntimeError('No valid plane triangles.')
    tri = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    # Do not silently fill invalid geometry with a new unconstrained triangulation.
    try:
        finder = tri.get_trifinder()
    except RuntimeError as exc:
        raise RuntimeError('Invalid/overlapping cut polygons; inspect cached section geometry.') from exc
    if abs(section.coordinate - args.centre[section.axis]) < args.diameter*.45:
        if int(finder(0., 0.)) >= 0:
            raise RuntimeError('The section includes the ball centre. Select the exterior air '
                               'cell zone only (the supplied case calls it "solid").')
    x = np.linspace(*args.xlim, args.pixels)
    nt = max(80, int(round(args.pixels * np.diff(args.transverse_lim)[0] / np.diff(args.xlim)[0])))
    t = np.linspace(*args.transverse_lim, nt)
    X, T = np.meshgrid(x, t)
    def interp(a):
        return mtri.LinearTriInterpolator(tri, a)(X, T).filled(np.nan)
    # Derived at the grid points from ALL THREE velocity components, not just
    # the in-plane components used for vector/LIC orientation.
    U = np.stack([interp(vel[:, k]) for k in range(3)], axis=-1)
    speed = np.linalg.norm(U, axis=-1)
    pressure = interp(p)
    q = .5 * args.rho * speed**2
    if np.isfinite(speed).sum() < .1 * speed.size:
        raise ValueError('Too little section data inside view. Check the view limits and plane.')
    edges = np.sort(np.vstack([triangles[:, [0, 1]], triangles[:, [1, 2]],
                              triangles[:, [2, 0]]]), axis=1)
    edges, multiplicity = np.unique(edges, axis=0, return_counts=True)
    boundary = points[edges[multiplicity == 1]]
    return dict(x=x, t=t, X=X, T=T, U=U, speed=speed, pressure=pressure,
                dynamic=q, total=pressure+q, boundary=boundary,
                a=U[..., 0], b=U[..., transverse_axis], section=section,
                transverse_label='z' if section.axis == 1 else 'y')


def bilinear(a, x, y):
    """Masked interpolation: stop at a hole/domain edge, without wrapping."""
    h, w = a.shape
    finite = np.isfinite(x) & np.isfinite(y)
    inside = finite & (x >= 0) & (x < w-1) & (y >= 0) & (y < h-1)
    xx = np.where(inside, x, 0); yy = np.where(inside, y, 0)
    i = xx.astype(int); j = yy.astype(int)
    fx = xx-i; fy = yy-j
    v = (a[j, i]*(1-fx)*(1-fy) + a[j, i+1]*fx*(1-fy) +
         a[j+1, i]*(1-fx)*fy + a[j+1, i+1]*fx*fy)
    return np.where(inside, v, np.nan)


def lic_texture(a, b, dx, dy, half_length=24, seed=7):
    """Symmetric Hann-kernel LIC of seeded noise along a steady 2D field.

    Midpoint integration steps are one pixel in arc length. A symmetric LIC
    depicts alignment, not direction of travel. It is not a pathline plot.
    """
    valid = np.isfinite(a) & np.isfinite(b)
    px, py = a/dx, b/dy
    mag = np.hypot(px, py)
    valid &= mag > 1e-10
    px = np.where(valid, px / np.maximum(mag, 1e-30), np.nan)
    py = np.where(valid, py / np.maximum(mag, 1e-30), np.nan)
    Y, X = np.mgrid[:a.shape[0], :a.shape[1]].astype(float)
    noise = np.random.default_rng(seed).random(a.shape)
    noise[~valid] = np.nan
    accum = np.nan_to_num(noise); weight = valid.astype(float)
    for sign in (-1, 1):
        x, y, alive = X.copy(), Y.copy(), valid.copy()
        for step in range(1, half_length+1):
            vx, vy = bilinear(px, x, y), bilinear(py, x, y)
            mx, my = x + sign*.5*vx, y + sign*.5*vy
            vx, vy = bilinear(px, mx, my), bilinear(py, mx, my)
            norm = np.hypot(vx, vy)
            good = alive & np.isfinite(norm) & (norm > 1e-12)
            x = x + sign*vx/np.maximum(norm, 1e-30)
            y = y + sign*vy/np.maximum(norm, 1e-30)
            sample = bilinear(noise, x, y)
            alive = good & np.isfinite(sample)
            w = .5*(1+np.cos(np.pi*step/(half_length+1)))
            accum += np.where(alive, sample, 0)*w
            weight += alive*w
    texture = np.where(valid, accum/np.maximum(weight, 1e-30), np.nan)
    finite = texture[np.isfinite(texture)]
    if not finite.size:
        raise ValueError('No nonzero in-plane velocity for LIC.')
    low, high = np.quantile(finite, [.01, .99])
    return np.clip((texture-low)/max(high-low, 1e-12), 0, 1)


DESCRIPTIONS = {
    'pressure': ('Static pressure', 'Gauge static pressure [Pa]', 'RdBu_r',
        'Inspect pressure differences around the ball: they contribute to aerodynamic force.'),
    'total': ('Total pressure', 'p + ½ρ|V|² [Pa, gauge]', 'cividis',
        'Inspect stagnation-pressure variations and the wake; the rotating wall can exchange mechanical energy.'),
    'dynamic': ('Dynamic pressure', '½ρ|V|² [Pa]', 'magma',
        'Dynamic pressure is kinetic energy per unit volume; it uses all three velocity components.'),
    'speed': ('Velocity magnitude', '|V| [m/s]', 'viridis',
        'Compare accelerated flow beside the ball with the downstream speed deficit.'),
    'vectors': ('Velocity vectors', '|V| [m/s]', 'viridis',
        'Arrows show in-plane velocity direction and magnitude; colour shows full 3D speed.'),
    'lic': ('Velocity line integral convolution', '|V| [m/s]', 'viridis',
        'Texture follows the projected velocity field. Symmetric LIC shows alignment, not the sign of motion.'),
}


def colour_scales(grids):
    scales = {}
    for field in ('pressure', 'total', 'dynamic', 'speed'):
        v = np.concatenate([g[field][np.isfinite(g[field])] for g in grids])
        lo, hi = float(v.min()), float(v.max())
        if field == 'pressure':
            limit = max(abs(lo), abs(hi), 1e-12)
            scales[field] = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
        else:
            lo = 0. if field in ('speed', 'dynamic') else lo
            scales[field] = Normalize(lo, max(hi, lo+max(abs(lo)*1e-8, 1e-12)))
    return scales


def render_one(g, plot, scales, args, folder):
    title, units, cmap, description = DESCRIPTIONS[plot]
    scalar = 'speed' if plot in ('vectors', 'lic') else plot
    norm = scales[scalar]
    fig = plt.figure(figsize=(14, 7.8))
    ax = fig.add_axes([.075, .235, .77, .575])
    cax = fig.add_axes([.875, .25, .017, .54])
    data = np.ma.masked_invalid(g[scalar])
    dx = g['x'][1]-g['x'][0]; dy = g['t'][1]-g['t'][0]
    extent = [g['x'][0]-dx/2, g['x'][-1]+dx/2,
              g['t'][0]-dy/2, g['t'][-1]+dy/2]
    if plot == 'lic':
        texture = lic_texture(g['a'], g['b'], g['x'][1]-g['x'][0],
                              g['t'][1]-g['t'][0], args.lic_length)
        rgba = plt.get_cmap(cmap)(norm(data.filled(np.nan)))
        rgba[..., :3] *= (.45+.55*np.nan_to_num(texture, nan=.5))[..., None]
        rgba[..., 3] = np.isfinite(g['speed']).astype(float)
        ax.imshow(rgba, origin='lower', extent=extent, interpolation='nearest')
        image = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    else:
        image = ax.imshow(data, norm=norm, cmap=cmap, origin='lower',
                          extent=extent, interpolation='nearest')
        if plot != 'vectors' and np.ptp(data.compressed()) > max(1e-9, abs(norm.vmax)*1e-8):
            ax.contour(g['X'], g['T'], data, levels=np.linspace(norm.vmin, norm.vmax, 13)[1:-1],
                       colors='#243b53', linewidths=.3, alpha=.23)
    if plot == 'vectors':
        sy = max(1, len(g['t'])//18); sx = max(1, len(g['x'])//44)
        cut = (slice(sy//2, None, sy), slice(sx//2, None, sx))
        # One U_infinity arrow is 0.25D long, on an equal-aspect x/D plane.
        qv = ax.quiver(g['X'][cut], g['T'][cut],
                       np.ma.masked_invalid(g['a'][cut]), np.ma.masked_invalid(g['b'][cut]),
                       color='#172b3a', angles='xy', scale_units='xy',
                       scale=args.speed/.25, pivot='mid', width=.0018,
                       headwidth=3.2, headlength=4.2)
        ax.quiverkey(qv, .68, 1.06, args.speed, f'{args.speed:g} m/s (in-plane)',
                     coordinates='axes', labelpos='E', fontproperties={'size': 10})
    ax.add_collection(LineCollection(g['boundary'], colors='#243b53', linewidths=.7))
    ax.set(xlim=args.xlim, ylim=args.transverse_lim,
           xlabel='(x − ball-centre x)/D',
           ylabel=f'({g["transverse_label"]} − ball-centre {g["transverse_label"]})/D')
    ax.set_aspect('equal', adjustable='box')
    ax.annotate('Incoming air: +x', xy=(.23, 1.035), xytext=(.02, 1.035),
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops={'arrowstyle': '->', 'color': '#176b61', 'lw': 1.7},
                fontsize=10, color='#176b61', annotation_clip=False)
    ax.annotate('Ball section', xy=(0, 0), xytext=(.14, .94),
                textcoords='axes fraction', arrowprops={'arrowstyle': '->', 'color': '#384b60'},
                fontsize=10, color='#243b53', bbox={'fc': 'white', 'ec': 'none', 'alpha': .83})
    fig.colorbar(image, cax=cax, label=units)
    plane = 'y' if g['section'].axis == 1 else 'z'
    label = f'{g["section"].name.capitalize()} section · {plane} = {g["section"].coordinate:g} m'
    fig.text(.075, .94, f'{title} | {label}', fontsize=18, weight='bold', color='#244667')
    spin = 'spin not specified' if args.omega is None else f'Ω = {args.omega:g} rad/s about +y'
    case_label = args.label or args.case.stem
    fig.text(.075, .885, f'{case_label}   ·   U∞ = {args.speed:g} m/s   ·   {spin}   ·   D = {args.diameter:g} m',
             fontsize=11, color='#526172')
    fig.text(.075, .125, description, fontsize=11, color='#244667', wrap=True)
    fig.text(.075, .068, 'Exterior air only · Equal geometric scales · Shared colour limits between planes · '
             f'ρ = {args.rho:g} kg/m³', fontsize=10, color='#657486')
    if args.demo:
        fig.text(.5, .49, 'SYNTHETIC DEMO — NOT CFD RESULTS', ha='center',
                 color='#aa3300', alpha=.65, fontsize=22, rotation=10)
    dest = folder / f'{g["section"].name}_{plot}.png'
    fig.savefig(dest, dpi=args.dpi, facecolor='white')
    plt.close(fig)
    print(dest, flush=True)
    return dest.name


def demo_sections(args):
    """Analytic potential-flow field, ONLY for renderer checks; never CFD data."""
    # An annular structured mesh has an exact topology hole and enough extent
    # for the requested rectangular view, without a Delaunay bridge across it.
    sections = []
    angles = np.linspace(0, 2*np.pi, 241)[:-1]
    radii = np.geomspace(.501, 14, 130)
    rr, aa = np.meshgrid(radii, angles, indexing='ij')
    xx, tt = rr*np.cos(aa), rr*np.sin(aa)
    vx = args.speed*(1-.25/rr**2*np.cos(2*aa))
    vt = -args.speed*.25/rr**2*np.sin(2*aa)
    p = .5*args.rho*(args.speed**2-vx**2-vt**2)
    faces = []
    n = len(angles)
    for j in range(len(radii)-1):
        for i in range(n):
            k = (i+1)%n
            faces.append(np.array([j*n+i, (j+1)*n+i, (j+1)*n+k, j*n+k]))
    for label, axis, coord in [('horizontal', 1, args.horizontal_y), ('vertical', 2, args.vertical_z)]:
        transverse = 2 if axis == 1 else 1
        xyz = np.tile(args.centre, (xx.size, 1)).astype(float)
        xyz[:, 0] += xx.ravel()*args.diameter
        xyz[:, transverse] += tt.ravel()*args.diameter
        xyz[:, axis] = coord
        vel = np.zeros_like(xyz); vel[:, 0] = vx.ravel(); vel[:, transverse] = vt.ravel()
        sections.append(Section(label, axis, coord, xyz, faces, p.ravel(), vel))
    return sections


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--case', type=Path, help='Solved .cas.h5 file')
    p.add_argument('--data', type=Path, help='Matching solved .dat.h5 file')
    p.add_argument('--speed', type=float, required=True, help='Actual inlet speed, m/s; annotation/vector scale only')
    p.add_argument('--omega', type=float, help='Actual spin, rad/s about +y; annotation only')
    p.add_argument('--label', help='Short case label, e.g. U25_spin50')
    p.add_argument('--zone', default='solid', help='Exterior FLUID cell-zone name (supplied case: solid)')
    p.add_argument('--material', default='air', help='Constant-density fluid material to verify')
    p.add_argument('--diameter', type=float, default=.229)
    p.add_argument('--rho', type=float, default=1.225, help='Constant air density used for q and p+q')
    p.add_argument('--centre', type=float, nargs=3, default=(0, 0, 0), metavar=('X', 'Y', 'Z'))
    p.add_argument('--horizontal-y', type=float, default=0., help='Horizontal section coordinate, m')
    p.add_argument('--vertical-z', type=float, default=0., help='Vertical section coordinate, m')
    p.add_argument('--xlim', type=float, nargs=2, default=(-2, 8), metavar=('MIN_D', 'MAX_D'))
    p.add_argument('--transverse-lim', type=float, nargs=2, default=(-2, 2), metavar=('MIN_D', 'MAX_D'))
    p.add_argument('--pressure-field', default='pressure', help='Fluent gauge static-pressure field key')
    p.add_argument('--pixels', type=int, default=900, help='Sampling grid width; increase for finer LIC detail')
    p.add_argument('--lic-length', type=int, default=24, help='Convolution half-length in pixels')
    p.add_argument('--dpi', type=int, default=220)
    p.add_argument('--processors', type=int, default=4)
    p.add_argument('--version', help='Fluent release, e.g. 25.2; otherwise use installed default')
    p.add_argument('--output', type=Path, default=Path('football_flow_graphics'))
    p.add_argument('--from-cache', type=Path, help='Replot an extracted sections folder without Fluent')
    p.add_argument('--demo', action='store_true', help='Synthetic renderer check; outputs visibly watermarked')
    return p


def main(argv=None):
    p = parser(); args = p.parse_args(argv)
    if not args.demo and not args.from_cache:
        if not args.case or not args.data:
            p.error('--case and matching --data are required (or use --from-cache).')
        if not args.case.is_file() or not args.data.is_file():
            p.error('Case or data file does not exist.')
    if args.speed <= 0 or args.diameter <= 0 or args.rho <= 0:
        p.error('Speed, diameter and density must be positive.')
    if args.pixels < 100 or args.lic_length < 1 or args.processors < 1 or args.dpi < 50:
        p.error('Use pixels >=100, LIC length >=1, processors >=1 and DPI >=50.')
    if args.xlim[0] >= args.xlim[1] or args.transverse_lim[0] >= args.transverse_lim[1]:
        p.error('View limits must increase.')
    if args.demo and args.from_cache:
        p.error('Choose --demo or --from-cache, not both.')
    # Fresh directory keeps the current export from mixing with older plot sets.
    folder = args.output.resolve(); folder.mkdir(parents=True, exist_ok=False)
    load_libraries()
    manifest = {'created': datetime.now().isoformat(), 'configuration': vars(args).copy(),
                'pressure_convention': 'gauge static p; q=0.5*rho*|V|^2; total=p+q',
                'velocity_convention': 'absolute velocity in the stationary CFD frame',
                'sampling': 'original section polygons; linear nodal interpolation; no hole filling',
                'lic': 'Python symmetric Hann-kernel LIC of in-plane velocity, not a Fluent screenshot',
                'synthetic': args.demo}
    try:
        if args.demo:
            sections = demo_sections(args)
        elif args.from_cache:
            metadata = json.loads((args.from_cache / 'manifest.json').read_text())
            manifest['synthetic'] = metadata.get('synthetic', False)
            args.demo = manifest['synthetic']
            # Metadata labels must follow the original extracted case.
            original = metadata['configuration']
            for key in ('speed', 'omega', 'diameter', 'rho', 'centre', 'label', 'case', 'data', 'zone'):
                if key in original:
                    setattr(args, key, original[key])
            sections = [load_section(args.from_cache / f'{s}_section.npz')
                        for s in ('horizontal', 'vertical')]
        else:
            sections, extraction = extract(args, folder)
            manifest.update(extraction)
        for s in sections:
            save_section(s, folder / f'{s.name}_section.npz')
        if args.case is None:
            args.case = Path('synthetic-demo')
        else:
            args.case = Path(args.case)
        manifest['configuration'] = vars(args).copy()
        grids = [section_grid(s, args) for s in sections]
        scales = colour_scales(grids)
        manifest['colour_limits'] = {k: [float(v.vmin), float(v.vmax)] for k, v in scales.items()}
        outputs = []
        for g in grids:
            plots = ('pressure', 'total', 'dynamic', 'speed', 'vectors', 'lic') if g['section'].name == 'horizontal' else ('pressure', 'speed', 'vectors')
            outputs.extend(render_one(g, plot, scales, args, folder) for plot in plots)
        manifest['outputs'] = outputs
        manifest['status'] = 'complete'
    except Exception as exc:
        manifest['status'] = 'failed'
        manifest['error'] = repr(exc)
        raise
    finally:
        (folder / 'manifest.json').write_text(json.dumps(manifest, indent=2, default=str), encoding='utf-8')
    print(f'\nSaved {len(outputs)} PNGs to {folder}')


if __name__ == '__main__':
    main()
