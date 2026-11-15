"""Offline numerical checks; no Fluent licence or actual CFD results required."""
from pathlib import Path
import tempfile
from types import SimpleNamespace
import football_flow_graphics as f

f.load_libraries()
np = f.np

# The LIC kernel is symmetric: reversing all vectors must preserve the texture.
a = np.ones((51, 91)); b = np.zeros_like(a)
a[20:30, 40:50] = np.nan; b[20:30, 40:50] = np.nan
forward = f.lic_texture(a, b, 1, 1, half_length=12)
reverse = f.lic_texture(-a, -b, 1, 1, half_length=12)
np.testing.assert_allclose(forward, reverse, atol=1e-12, equal_nan=True)
assert np.isnan(forward[20:30, 40:50]).all()
horizontal_variation = np.nanmean(np.abs(np.diff(forward, axis=1)))
vertical_variation = np.nanmean(np.abs(np.diff(forward, axis=0)))
assert horizontal_variation < vertical_variation*.5  # true flow-aligned texture

args = SimpleNamespace(centre=(0, 0, 0), diameter=.229, speed=25., rho=1.225,
                       horizontal_y=0., vertical_z=0., pixels=240,
                       xlim=(-2, 8), transverse_lim=(-2, 2))
sections = f.demo_sections(args)
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp)/'plane.npz'
    f.save_section(sections[0], path)
    section = f.load_section(path)
    np.testing.assert_array_equal(section.vertices, sections[0].vertices)
    for original, restored in zip(section.faces, sections[0].faces):
        np.testing.assert_array_equal(original, restored)
    # A nonzero out-of-plane component must contribute to q and speed.
    section.velocity[:, 1] = 3.5
    g = f.section_grid(section, args)
    np.testing.assert_allclose(g['dynamic'], .5*args.rho*np.sum(g['U']**2, axis=-1),
                               equal_nan=True)
    np.testing.assert_allclose(g['total'], g['pressure']+g['dynamic'], equal_nan=True)
    centre = np.hypot(g['X'], g['T']) < .48
    assert np.isnan(g['speed'][centre]).all()  # no interpolation across the ball
    assert np.isfinite(g['speed'][np.hypot(g['X'], g['T']) > .55]).all()
    # Welding identical per-partition vertices must preserve the same result.
    n = len(section.vertices)
    duplicated = f.Section(section.name, section.axis, section.coordinate,
                          np.vstack([section.vertices]*2),
                          section.faces + [face+n for face in section.faces],
                          np.tile(section.pressure, 2), np.vstack([section.velocity]*2))
    again = f.section_grid(duplicated, args)
    np.testing.assert_allclose(again['speed'], g['speed'], equal_nan=True)
print('Passed: LIC orientation, reversal symmetry, mask, original topology, cache round-trip, '
      '3D pressure definitions and duplicated partition geometry.')
