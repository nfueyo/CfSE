#!/usr/bin/env python3
"""Analyse the supplied Fluent elbow case using h5py, NumPy and Matplotlib.

    python3 -m pip install h5py numpy matplotlib
    python3 elbow_mesh_quality.py elbow_h2_v11Dec21.msh.cas.h5

The reader expects the single-section Fluent HDF5 layout used by this case.
The output includes a per-cell CSV and three PNG figures.
"""

import argparse
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter, MaxNLocator, NullFormatter


def read_mesh(path):
    """Read face loops and align adjacency arrays using Fluent's face IDs."""
    with h5py.File(path, "r") as case:
        mesh = case["meshes/1"]
        nodes = mesh["nodes/coords/1"][:]
        cell_types = mesh["cells/ctype/1/cell-types"][:]
        sizes = mesh["faces/nodes/1/nnodes"][:].astype(np.int64)
        vertices = mesh["faces/nodes/1/nodes"][:].astype(np.int64) - 1
        offsets = np.r_[0, np.cumsum(sizes)]
        owner = mesh["faces/c0/1"][:].astype(np.int64) - 1
        neighbours = mesh["faces/c1/1"]
        first = int(neighbours.attrs["minId"][0]) - 1
        last = int(neighbours.attrs["maxId"][0])
        other = np.full(len(owner), -1, dtype=np.int64)
        other[first:last] = neighbours[:].astype(np.int64) - 1

    if (
        len(sizes) != len(owner)
        or offsets[-1] != len(vertices)
        or np.any(sizes < 3)
        or np.any((vertices < 0) | (vertices >= len(nodes)))
        or np.any((owner < 0) | (owner >= len(cell_types)))
        or np.any((other < -1) | (other >= len(cell_types)))
        or np.any(owner == other)
    ):
        raise ValueError("Invalid connectivity for the supported Fluent layout.")
    return nodes, cell_types, vertices, offsets, owner, other


def cell_geometry(nodes, vertices, offsets, owner, other, cell_count):
    """Integrate signed tetrahedra from an origin to oriented face triangles."""
    origin = nodes.mean(axis=0)
    nodes = nodes - origin
    normals = np.empty((len(owner), 3))
    face_volume = np.empty(len(owner))
    face_moment = np.empty((len(owner), 3))

    for f, (start, stop) in enumerate(zip(offsets[:-1], offsets[1:])):
        ring = nodes[vertices[start:stop]]
        a = ring[0]
        b, c = ring[1:-1], ring[2:]
        cross = np.cross(b - a, c - a)
        normals[f] = cross.sum(axis=0)
        signed = cross @ a / 6
        face_volume[f] = signed.sum()
        face_moment[f] = signed @ (a + b + c) / 4

    interior = other >= 0
    volume = np.zeros(cell_count)
    moment = np.zeros((cell_count, 3))
    np.add.at(volume, owner, face_volume)
    np.add.at(volume, other[interior], -face_volume[interior])
    np.add.at(moment, owner, face_moment)
    np.add.at(moment, other[interior], -face_moment[interior])

    # The global winding convention may be reversed. Volume and moment
    # change sign together, so their ratio still gives the centroid.
    if np.all(volume < 0):
        volume = -volume
        moment = -moment
    if np.any(volume <= 0) or not np.isfinite(volume).all():
        raise ValueError("Found non-positive or invalid cell volumes. Check face winding.")
    return volume, moment / volume[:, None] + origin, normals


def neighbour_metrics(size, centres, normals, owner, other):
    """Attribute each interior face's metrics to both adjoining cells."""
    interior = other >= 0
    i, j = owner[interior], other[interior]
    delta = centres[j] - centres[i]
    normal = normals[interior]
    scale = np.linalg.norm(delta, axis=1) * np.linalg.norm(normal, axis=1)
    if np.any(scale <= 0):
        raise ValueError("Found a zero-area face or coincident neighbouring centroids.")

    # Absolute cosine compares an undirected normal axis with the line
    # joining the centroids, independently of face winding.
    cosine = np.abs(np.einsum("ij,ij->i", delta, normal)) / scale
    face_angle = np.degrees(np.arccos(np.clip(cosine, 0, 1)))
    face_jump = np.maximum(size[i] / size[j], size[j] / size[i])

    results = []
    for values in (face_jump, face_angle):
        worst = np.full(len(size), -np.inf)
        np.maximum.at(worst, i, values)
        np.maximum.at(worst, j, values)
        worst[~np.isfinite(worst)] = np.nan
        results.append(worst)
    return results


def plot_results(output, centres, cell_types, metrics, slab_centres, half_width):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), layout="constrained")
    for index, (ax, (label, values)) in enumerate(zip(axes, metrics)):
        valid = np.isfinite(values)
        bins = np.histogram_bin_edges(values[valid], bins=30)
        if index == 0 and bins[0] > 0:
            bins = np.geomspace(bins[0], bins[-1], 31)
            ax.set_xscale("log")
            ax.xaxis.set_major_formatter(LogFormatter())
            ax.xaxis.set_minor_formatter(NullFormatter())
        for code, name in [(4, "Hexahedra"), (7, "Polyhedra")]:
            sample = values[valid & (cell_types == code)]
            if sample.size:
                ax.hist(sample, bins=bins, weights=np.full(sample.size, 100 / sample.size),
                        histtype="step", linewidth=1.6, label=name)
        ax.set(xlabel=label, ylabel="Cells of each type (%)")
    axes[0].legend()
    fig.savefig(output / "histograms.png", dpi=180)
    plt.close(fig)

    x, y, z = (centres * 1000).T
    fig, axes = plt.subplots(len(slab_centres), 3, squeeze=False,
                             figsize=(12, 3.5 * len(slab_centres)), layout="constrained")
    for column, (label, values) in enumerate(metrics):
        valid = np.isfinite(values)
        lo, hi = np.min(values[valid]), np.max(values[valid])
        for ax, target in zip(axes[:, column], slab_centres):
            selected = valid & (np.abs(z - target) <= half_width)
            points = ax.scatter(x[selected], y[selected], c=values[selected], s=7,
                                vmin=lo, vmax=hi, cmap="viridis", linewidths=0)
            ax.set(xlabel="x (mm)", ylabel="y (mm)", xlim=(x.min(), x.max()),
                   ylim=(y.min(), y.max()), aspect="equal",
                   title=f"z = {target:g} ± {half_width:g} mm · {selected.sum():,} cells")
        fig.colorbar(points, ax=axes[:, column], label=label, shrink=0.85)
    fig.savefig(output / "centroid_slabs.png", dpi=180)
    plt.close(fig)

    angles = metrics[-1][1]
    valid = np.flatnonzero(np.isfinite(angles))
    count = max(1, int(np.ceil(0.01 * len(valid))))
    worst = valid[np.argsort(angles[valid], kind="stable")[-count:]]
    fig = plt.figure(figsize=(8, 6), layout="constrained")
    ax = fig.add_subplot(projection="3d", computed_zorder=False)
    ax.scatter(x, y, z, c="#64748b", s=2, alpha=0.08, depthshade=False)
    points = ax.scatter(x[worst], y[worst], z[worst], c=angles[worst], cmap="plasma",
                        s=13, depthshade=False)
    ax.set(xlabel="x (mm)", ylabel="y (mm)", zlabel="z (mm)",
           title=f"Largest non-orthogonality angles · {count:,} cells (top 1%)")
    ax.set_box_aspect(np.ptp(centres, axis=0))
    ax.zaxis.set_major_locator(MaxNLocator(2))
    fig.colorbar(points, ax=ax, label=metrics[-1][0], shrink=0.65)
    fig.savefig(output / "worst_cells_3d.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path, help="Fluent elbow .cas.h5 file")
    parser.add_argument("--output", type=Path, default=Path("elbow_mesh_results"))
    parser.add_argument("--slab-z-mm", type=float, nargs="+", default=[0, 15])
    parser.add_argument("--slab-half-width-mm", type=float, default=2)
    args = parser.parse_args()
    if args.slab_half_width_mm <= 0:
        parser.error("--slab-half-width-mm must be positive")

    nodes, types, vertices, offsets, owner, other = read_mesh(args.case)
    volume, centres, normals = cell_geometry(nodes, vertices, offsets, owner, other, len(types))
    size = np.cbrt(volume)
    jump, angle = neighbour_metrics(size, centres, normals, owner, other)
    metrics = [
        ("Characteristic cell size (mm)", size * 1000),
        ("Largest neighbour size ratio", jump),
        ("Largest face angle (°)", angle),
    ]

    args.output.mkdir(parents=True, exist_ok=True)
    table = np.column_stack((np.arange(1, len(types) + 1), types, centres, volume,
                             size * 1000, jump, angle))
    np.savetxt(args.output / "cell_metrics.csv", table, delimiter=",",
               fmt=["%d", "%d"] + ["%.10g"] * 7, comments="",
               header="cell_id,fluent_cell_type,x_m,y_m,z_m,volume_m3,size_mm,max_size_ratio,max_angle_deg")
    plot_results(args.output, centres, types, metrics, args.slab_z_mm, args.slab_half_width_mm)

    print(f"{len(types):,} cells, {len(owner):,} faces, {len(nodes):,} nodes")
    print(f"{np.count_nonzero(types == 4):,} hexahedra, {np.count_nonzero(types == 7):,} polyhedra")
    for label, values in metrics:
        low, median, p95, high = np.nanpercentile(values, [0, 50, 95, 100])
        print(f"{label}: min {low:.3g}, median {median:.3g}, 95th percentile {p95:.3g}, max {high:.3g}")
    print(f"Results saved in {args.output.resolve()}")


if __name__ == "__main__":
    main()
