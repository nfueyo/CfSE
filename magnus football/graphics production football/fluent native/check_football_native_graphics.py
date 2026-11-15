"""Offline checks of scientific framing and image format; no mock CFD output."""
import math
from pathlib import Path
from types import SimpleNamespace
import tempfile
import struct
from football_fluent_native_graphics import (
    annotation_collection, attach_annotations, camera_spec, choice, png_size,
    select_range_mode,
)

def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

args = SimpleNamespace(diameter=.229, centre=(1., 2., 3.), xlim=(-2., 8.),
                       transverse_lim=(-2., 2.), horizontal_y=2., vertical_z=3.,
                       width=3200, height=1800)
for plane in ('horizontal', 'vertical'):
    c = camera_spec(args, plane)
    look = [t-p for p, t in zip(c['position'], c['target'])]
    right = cross(look, c['up_vector'])
    assert right[0] > 0 and right[1:] == [0., 0.], 'View must not mirror the downstream axis.'
    assert sum(a*b for a,b in zip(look, c['up_vector'])) == 0
    w, h = c['field']
    assert math.isclose(w/h, args.width/args.height), 'No geometric stretching.'
    assert w >= 10*args.diameter and h >= 4*args.diameter
    transverse = 2 if plane == 'horizontal' else 1
    assert math.isclose(c['target'][transverse], args.centre[transverse])
    assert math.isclose(c['target'][0], 1.+3.*args.diameter)

with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp)/'header.png'
    # Exercise header validation only; this is not an image or CFD result.
    path.write_bytes(b'\x89PNG\r\n\x1a\n' + b'\x00\x00\x00\rIHDR' + struct.pack('>II', 3200, 1800))
    assert png_size(path) == [3200, 1800]
    path.write_bytes(b'not a png')
    try:
        png_size(path)
    except RuntimeError:
        pass
    else:
        raise AssertionError('Invalid image signature was accepted.')

class Setting:
    def __init__(self, values=None):
        self.values = values
        self.active = True
    def allowed_values(self): return self.values
    def set_state(self, value): self.value = value
    def is_active(self): return self.active
setting = Setting(['Orthographic', 'Perspective'])
assert choice(setting, 'orthographic') == 'Orthographic'
try:
    choice(setting, 'unsupported')
except RuntimeError:
    pass
else:
    raise AssertionError('Unsupported native option must fail clearly.')

auto_range_on = SimpleNamespace(child_names=['global_range'], global_range=Setting())
auto_range_off = SimpleNamespace(child_names=[],)
range_setting = SimpleNamespace(
    child_names=['option', 'auto_range_on', 'auto_range_off'],
    option=Setting(['auto-range-on', 'auto-range-off']),
    auto_range_on=auto_range_on,
    auto_range_off=auto_range_off,
)
range_obj = SimpleNamespace(child_names=['range_option'], range_option=range_setting)
assert select_range_mode(range_obj, automatic=True) == 'auto-range-on'
assert auto_range_on.global_range.value is False
assert select_range_mode(range_obj, automatic=False) == 'auto-range-off'

legacy_range_obj = SimpleNamespace(
    child_names=['range_option'],
    range_option=Setting(['local', 'custom']),
)
assert select_range_mode(legacy_range_obj, automatic=True) == 'local'

annotations = object()
graphics = SimpleNamespace(child_names=['annotations'], annotations=annotations)
assert annotation_collection(graphics) is annotations
assert annotation_collection(SimpleNamespace(child_names=['contour'])) is None
assert not attach_annotations(None, object(), 'test', 'horizontal', 'pressure', args)
print('Passed: cameras, PNG headers, option selection, optional annotations, and range-mode compatibility.')
