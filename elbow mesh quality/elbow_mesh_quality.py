#!/usr/bin/env python3
"""Analyse the supplied Fluent elbow case using h5py, NumPy and Matplotlib.

    python3 -m pip install h5py numpy matplotlib
    python3 elbow_mesh_quality.py elbow_h2_v11Dec21.msh.cas.h5

The reader expects the single-section Fluent HDF5 layout used by this case.
The output includes a per-cell CSV, three PNG figures, and a student report PDF.
"""

import argparse
import textwrap
from datetime import datetime
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
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


def write_report(output, centres, cell_types, metrics, slab_centres, half_width):
    size, jump, angle = (values for _, values in metrics)
    cell_count = len(cell_types)
    type_counts = {code: np.count_nonzero(cell_types == code) for code in (4, 7)}
    summaries = {
        name: np.nanpercentile(values, [0, 50, 95, 100])
        for name, values in (("size", size), ("ratio", jump), ("angle", angle))
    }
    valid_angles = np.flatnonzero(np.isfinite(angle))
    worst_count = max(1, int(np.ceil(0.01 * len(valid_angles))))
    worst = valid_angles[np.argsort(angle[valid_angles], kind="stable")[-worst_count:]][::-1]
    centres_mm = centres * 1000

    navy, blue, teal = "#17324d", "#2677a8", "#148a8a"
    ink, muted, rule = "#213547", "#526579", "#c9d5de"
    page_number = 0

    def new_page(title, subtitle):
        nonlocal page_number
        page_number += 1
        fig = plt.figure(figsize=(11.69, 8.27), facecolor="white")
        fig.text(0.055, 0.948, title, fontsize=21, weight="bold",
                 color=navy, va="top")
        fig.text(0.056, 0.901, subtitle, fontsize=9.5, color=muted, va="top")
        fig.add_artist(Line2D([0.055, 0.945], [0.878, 0.878],
                              color=rule, linewidth=1))
        fig.add_artist(Line2D([0.055, 0.945], [0.064, 0.064],
                              color=rule, linewidth=0.8))
        fig.text(0.055, 0.034, "ELBOW MESH QUALITY  |  STUDENT GUIDE",
                 fontsize=7.5, color=muted, va="center")
        fig.text(0.945, 0.034, f"{page_number:02d}", fontsize=8,
                 color=muted, ha="right", va="center")
        return fig

    def wrapped(text, width):
        return "\n".join(
            textwrap.fill(paragraph, width=width) if paragraph else ""
            for paragraph in text.split("\n")
        )

    def add_table(fig, bbox, headers, rows, widths, fontsize=8.5):
        ax = fig.add_axes(bbox)
        ax.axis("off")
        table = ax.table(
            cellText=rows, colLabels=headers, colWidths=widths,
            cellLoc="left", colLoc="left", bbox=[0, 0, 1, 1],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(fontsize)
        for (row, _), cell in table.get_celld().items():
            cell.set_edgecolor("white")
            cell.set_linewidth(0.6)
            if row == 0:
                cell.set_facecolor(navy)
                cell.get_text().set_color("white")
                cell.get_text().set_weight("bold")
            else:
                cell.set_facecolor("#f0f5f8" if row % 2 else "white")
                cell.get_text().set_color(ink)
                cell.PAD = 0.08

    def add_image(fig, filename, bbox):
        ax = fig.add_axes(bbox)
        ax.imshow(plt.imread(output / filename))
        ax.set_axis_off()

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "text.color": ink,
        "axes.labelcolor": ink,
        "xtick.color": "#425466",
        "ytick.color": "#425466",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })

    with PdfPages(
        output / "mesh_quality_report.pdf",
        metadata={"Title": "Elbow mesh quality: student guide"},
    ) as pdf:
        fig = new_page(
            "Elbow mesh quality",
            "A student guide to the terminology, metrics and evidence in the mesh analysis",
        )
        fig.text(0.06, 0.83, "What this report is for", fontsize=13,
                 weight="bold", color=blue)
        intro = (
            "A computational mesh divides a fluid region into small control "
            "volumes. This report explains three cell-level indicators: "
            "volume-based size, the largest size change across a cell's "
            "neighbours, and the largest angle between a face-normal axis "
            "and the line joining neighbouring cell centres. The indicators "
            "complement one another; none is a single all-purpose "
            "mesh-quality score."
        )
        fig.text(0.06, 0.775, wrapped(intro, 65), fontsize=10.5, color=ink,
                 va="top", linespacing=1.55)
        fig.text(0.06, 0.515, "Dataset snapshot", fontsize=13,
                 weight="bold", color=blue)
        fig.text(0.06, 0.47, f"{cell_count:,} cells in the per-cell table",
                 fontsize=11, weight="bold", color=navy)
        fig.text(
            0.06, 0.425,
            f"{type_counts[4]:,} hexahedra "
            f"({100 * type_counts[4] / cell_count:.2f}%)  |  "
            f"{type_counts[7]:,} polyhedra "
            f"({100 * type_counts[7] / cell_count:.2f}%)",
            fontsize=9.5, color=ink,
        )
        fig.text(0.06, 0.36, "Main lesson", fontsize=13, weight="bold",
                 color=blue)
        lesson = (
            "Most cells have modest neighbour-size changes and low computed "
            "angles, but the upper tail contains a small, spatially clustered "
            "set of larger angles. Use the location and cell IDs to inspect "
            "the mesh; do not label it pass/fail from one maximum alone."
        )
        fig.text(0.06, 0.315, wrapped(lesson, 64), fontsize=10.2,
                 color=ink, va="top", linespacing=1.5)
        fig.text(0.56, 0.83, "Three metrics at a glance", fontsize=13,
                 weight="bold", color=blue)
        overview = [
            ["Characteristic size", f"{summaries['size'][1]:.3f} mm",
             f"{summaries['size'][2]:.3f} mm",
             f"{summaries['size'][3]:.3f} mm"],
            ["Neighbour size ratio", f"{summaries['ratio'][1]:.3f}",
             f"{summaries['ratio'][2]:.3f}",
             f"{summaries['ratio'][3]:.3f}"],
            ["Largest face angle", f"{summaries['angle'][1]:.2f} deg",
             f"{summaries['angle'][2]:.2f} deg",
             f"{summaries['angle'][3]:.2f} deg"],
        ]
        add_table(fig, [0.55, 0.50, 0.40, 0.26],
                  ["Metric", "Median", "95th pct.", "Maximum"], overview,
                  [0.38, 0.20, 0.21, 0.21], 8.3)
        fig.text(0.56, 0.435, "How to use the table", fontsize=11,
                 weight="bold", color=navy)
        percentile_note = (
            "The median describes a typical cell; the 95th percentile "
            "describes the upper tail; the maximum identifies the single "
            "worst reported value. A maximum should trigger inspection, "
            "not an automatic verdict."
        )
        fig.text(0.56, 0.395, wrapped(percentile_note, 53), fontsize=9.5,
                 color=ink, va="top", linespacing=1.45)
        fig.text(0.56, 0.25, "Scope", fontsize=11, weight="bold", color=navy)
        scope = (
            "This report summarizes the current run's CSV and graphics. "
            "It does not replace Fluent's solver-native quality checks or "
            "a mesh-sensitivity study."
        )
        fig.text(0.56, 0.21, wrapped(scope, 53), fontsize=9.2,
                 color=muted, va="top", linespacing=1.4)
        pdf.savefig(fig)
        plt.close(fig)

        fig = new_page(
            "1 | The mesh: vocabulary first",
            "The metrics make sense only when the mesh objects and their relationships are clear",
        )
        glossary = [
            ["Node / vertex", "A point used to describe cell and face geometry. A node is not a cell."],
            ["Face", "A polygonal boundary shared by two interior cells, or separating a cell from the outside/boundary."],
            ["Cell / control volume", "A 3-D region over which the finite-volume method balances mass, momentum and energy."],
            ["Cell centroid", "The volume-weighted geometric centre of a cell, computed here from face geometry."],
            ["Owner and neighbour", "The two cells attached to an interior face in Fluent connectivity."],
            ["Face-area normal", "A vector perpendicular to a face, with magnitude proportional to face area; its sign depends on winding."],
            ["Hexahedron / polyhedron", "Cell topologies recorded here as Fluent type codes 4 and 7."],
            ["Interior-only metric", "A value formed from faces with two adjacent cells. Boundary faces do not contribute to the two neighbour metrics."],
        ]
        add_table(fig, [0.055, 0.32, 0.89, 0.52],
                  ["Term", "Meaning in this analysis"], glossary,
                  [0.22, 0.78], 8.2)
        ax = fig.add_axes([0.06, 0.095, 0.43, 0.19])
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 5)
        ax.axis("off")
        ax.add_patch(Polygon(
            [(1, 1), (4.1, 1.25), (4.55, 4.15), (1, 4)],
            closed=True, facecolor="#e7f1f7", edgecolor=blue, linewidth=1.5,
        ))
        ax.add_patch(Polygon(
            [(4.1, 1.25), (8.5, 1.65), (8.1, 4.5), (4.55, 4.15)],
            closed=True, facecolor="#e9f5f2", edgecolor=teal, linewidth=1.5,
        ))
        c1, c2, face_mid = (2.5, 2.5), (6.4, 2.95), (4.325, 2.7)
        ax.plot([c1[0], c2[0]], [c1[1], c2[1]], color=navy, linewidth=1.5)
        ax.scatter([c1[0], c2[0]], [c1[1], c2[1]], s=24, color=navy, zorder=4)
        ax.text(2.05, 1.95, "centroid i", fontsize=8, color=navy)
        ax.text(6.1, 2.4, "centroid j", fontsize=8, color=navy)
        ax.annotate(
            "", xy=(6.0, 2.89), xytext=(4.65, 2.73),
            arrowprops=dict(arrowstyle="->", color=navy, linewidth=1.5),
        )
        ax.text(4.75, 4.62, "centroid line (delta)", fontsize=7.2, color=navy)
        ax.annotate(
            "", xy=(5.55, 2.51), xytext=face_mid,
            arrowprops=dict(arrowstyle="->", color=teal, linewidth=1.8),
        )
        ax.text(4.65, 1.6, "face normal n", fontsize=7.5, color=teal)
        ax.text(
            0.12, 0.06,
            "Concept sketch (not to scale): angle compares the centroid line with the face-normal axis.",
            transform=ax.transAxes, fontsize=7.4, color=muted,
        )
        fig.text(0.53, 0.255, "Why the geometry matters", fontsize=11,
                 weight="bold", color=navy)
        geometry_note = (
            "For a well-aligned face, the normal points along the line between "
            "the two adjacent cell centroids. The angle metric quantifies "
            "their misalignment. It is not the angle between cell edges, nor "
            "a visual judgement of how 'crooked' a cell looks."
        )
        fig.text(0.53, 0.215, wrapped(geometry_note, 57), fontsize=9.1,
                 color=ink, va="top", linespacing=1.45)
        pdf.savefig(fig)
        plt.close(fig)

        fig = new_page(
            "2 | What each metric means",
            "Definitions follow the calculations in elbow_mesh_quality.py; they are not interchangeable",
        )
        cards = [
            (
                "A. Characteristic cell size",
                "h_i = cube_root(V_i)",
                "V_i is the computed cell volume. The cube root is the side "
                "length of a cube with the same volume. It is a convenient "
                "length scale (reported in millimetres), not an actual edge "
                "length for a general polyhedron. Smaller values mean smaller "
                "control volumes, not automatically better geometry.",
            ),
            (
                "B. Largest neighbour size ratio",
                "R_i = max over adjacent j of max(h_i/h_j, h_j/h_i)",
                "R_i is at least 1: a value of 1 means equal size for that cell "
                "and neighbour. Larger values mean a stronger local change "
                "across an interior face. It does not measure cell shape.",
            ),
            (
                "C. Largest face angle (computed non-orthogonality proxy)",
                "theta_f = arccos(abs(delta dot n) / (|delta| |n|)); "
                "theta_i = max over interior faces",
                "Delta joins the two cell centroids; n is the face-area "
                "normal. Absolute cosine ignores the normal's sign. Zero "
                "degrees is ideal alignment; angles approaching 90 degrees "
                "indicate stronger misalignment. The cell value is its "
                "worst interior face.",
            ),
        ]
        y = 0.82
        for title, formula, explanation in cards:
            fig.text(0.07, y, title, fontsize=12, weight="bold",
                     color=blue, va="top")
            fig.text(0.075, y - 0.048, formula, fontsize=10.5,
                     weight="bold", color=navy, va="top")
            fig.text(0.075, y - 0.09, wrapped(explanation, 119),
                     fontsize=9.1, color=ink, va="top", linespacing=1.42)
            y -= 0.24
        fig.text(0.07, 0.17, "Percentiles in this report", fontsize=10.5,
                 weight="bold", color=navy)
        percentile_text = (
            "Median = 50th percentile (half the cells are at or below it); "
            "95th percentile = upper-tail summary (95% are at or below it); "
            "maximum = largest cell value. Read these together."
        )
        fig.text(0.07, 0.135, wrapped(percentile_text, 125), fontsize=8.8,
                 color=ink, va="top")
        pdf.savefig(fig)
        plt.close(fig)

        fig = new_page(
            "3 | Results for this mesh",
            f"The reported statistics summarize {cell_count:,} rows in the exported per-cell CSV",
        )
        fig.text(0.06, 0.835, "Whole-mesh distribution", fontsize=12,
                 weight="bold", color=blue)
        metric_rows = [
            ["Characteristic cell size", *(f"{v:.3f}" for v in summaries["size"]), "mm"],
            ["Largest neighbour size ratio", *(f"{v:.3f}" for v in summaries["ratio"]), "dimensionless"],
            ["Largest face angle", *(f"{v:.2f}" for v in summaries["angle"]), "degrees"],
        ]
        add_table(
            fig, [0.06, 0.615, 0.88, 0.18],
            ["Metric", "Minimum", "Median", "95th pct.", "Maximum", "Unit"],
            metric_rows, [0.30, 0.12, 0.12, 0.13, 0.12, 0.21], 8.2,
        )
        fig.text(0.06, 0.565, "By cell type", fontsize=12,
                 weight="bold", color=blue)
        type_rows = []
        for code, name in ((4, "Hexahedra"), (7, "Polyhedra")):
            selected = cell_types == code
            type_rows.append([
                name, f"{type_counts[code]:,}",
                f"{100 * type_counts[code] / cell_count:.2f}%",
                f"{np.nanmedian(size[selected]):.3f}",
                f"{np.nanpercentile(size[selected], 95):.3f}",
                f"{np.nanmedian(jump[selected]):.3f}",
                f"{np.nanmedian(angle[selected]):.2f}",
                f"{np.nanmax(angle[selected]):.2f}",
            ])
        add_table(
            fig, [0.06, 0.37, 0.88, 0.15],
            ["Type", "Cells", "Share", "Size med (mm)", "Size p95 (mm)",
             "Ratio med", "Angle med (deg)", "Angle max (deg)"],
            type_rows, [0.15, 0.10, 0.10, 0.13, 0.13, 0.12, 0.15, 0.12],
            7.3,
        )
        fig.text(0.06, 0.315, "What stands out", fontsize=12,
                 weight="bold", color=blue)
        high_ratio = np.count_nonzero(jump >= 2.0)
        high_angle = np.count_nonzero(angle >= 45.0)
        observations = [
            f"Polyhedra make up {100 * type_counts[7] / cell_count:.1f}% "
            "of the listed cells and have a wider angle tail than the "
            "hexahedra in this export.",
            f"{high_ratio:,} cells ({100 * high_ratio / cell_count:.2f}%) "
            "have a neighbour-size ratio of at least 2.0; "
            f"{high_angle:,} cells ({100 * high_angle / cell_count:.2f}%) "
            "have an angle of at least 45 degrees.",
            "These cut points are descriptive examples, not universal "
            "accept/reject thresholds. Use the solver, physics and intended "
            "discretization to decide what to inspect.",
        ]
        bullets = "\n".join(
            "• " + textwrap.fill(text, width=123, subsequent_indent="  ")
            for text in observations
        )
        fig.text(0.07, 0.275, bullets, fontsize=8.9, color=ink, va="top",
                 linespacing=1.45)
        pdf.savefig(fig)
        plt.close(fig)

        fig = new_page(
            "4 | Read the distributions, not just the worst cell",
            "Per-cell histograms, shown separately for each cell topology",
        )
        add_image(fig, "histograms.png", [0.055, 0.40, 0.89, 0.42])
        fig.text(0.07, 0.335, "How to read this figure", fontsize=11.5,
                 weight="bold", color=blue)
        histogram_note = (
            "Each outline is normalized within its own cell type: the y-axis "
            "is the percentage of that type in a bin. Their shapes are "
            "comparable, but bar heights do not show how many cells of each "
            "type exist. Cell size uses a logarithmic x-axis; ratio and angle "
            "use linear axes. A long right tail flags cells for spatial "
            "inspection; a narrow peak describes the common pattern."
        )
        fig.text(0.07, 0.29, wrapped(histogram_note, 132), fontsize=9,
                 color=ink, va="top", linespacing=1.45)
        fig.text(
            0.07, 0.145,
            "In this mesh, the size distributions differ by topology; "
            "polyhedra also carry the largest angle and size-ratio tails.",
            fontsize=9.1, weight="bold", color=navy, va="top", wrap=True,
        )
        pdf.savefig(fig)
        plt.close(fig)

        fig = new_page(
            "5 | Where are the values in space?",
            "Centroid maps for the requested z-slabs; each dot is a cell centroid",
        )
        add_image(fig, "centroid_slabs.png", [0.095, 0.18, 0.81, 0.64])
        slab_descriptions = []
        angle_valid = np.isfinite(angle)
        for target in slab_centres:
            selected = angle_valid & (
                np.abs(centres_mm[:, 2] - target) <= half_width
            )
            slab_descriptions.append(
                f"z = {target:g} ± {half_width:g} mm "
                f"({np.count_nonzero(selected):,} centroids)"
            )
        caption = (
            "Rows show " + "; ".join(slab_descriptions) + ". Columns are "
            "size, neighbour ratio and angle. Within each column, the colour "
            "scale is shared between slabs. These are sampled slices, not "
            "every cell in the 3-D mesh; blank regions are outside the "
            "selected slabs or fluid domain."
        )
        fig.text(0.06, 0.115, wrapped(caption, 140), fontsize=8.6,
                 color=ink, va="top", linespacing=1.35)
        pdf.savefig(fig)
        plt.close(fig)

        fig = new_page(
            "6 | Investigate the upper tail",
            "The 3-D plot highlights the top 1% of cells by computed largest face angle",
        )
        add_image(fig, "worst_cells_3d.png", [0.035, 0.19, 0.62, 0.63])
        fig.text(0.68, 0.82, "Five largest values", fontsize=11,
                 weight="bold", color=blue, va="top")
        cell_ids = np.arange(1, cell_count + 1)
        top_rows = [
            [
                str(cell_ids[index]),
                f"{centres_mm[index, 0]:.1f}",
                f"{centres_mm[index, 1]:.1f}",
                f"{centres_mm[index, 2]:.1f}",
                f"{angle[index]:.1f}",
                f"{jump[index]:.2f}",
            ]
            for index in worst[:5]
        ]
        add_table(
            fig, [0.675, 0.54, 0.29, 0.23],
            ["Cell ID", "x mm", "y mm", "z mm", "Angle deg", "Ratio"],
            top_rows, [0.20, 0.15, 0.15, 0.15, 0.18, 0.17], 7,
        )
        fig.text(0.68, 0.48, "How to use this view", fontsize=10.5,
                 weight="bold", color=navy, va="top")
        outlier_note = (
            f"The highlighted {worst_count:,} cells are the largest 1% of "
            f"the {len(valid_angles):,} valid angle values. Use the Fluent "
            "cell IDs to inspect local faces, neighbouring cells and "
            "boundary proximity before deciding whether refinement or "
            "repair is appropriate."
        )
        fig.text(0.68, 0.44, wrapped(outlier_note, 42), fontsize=8.6,
                 color=ink, va="top", linespacing=1.4)
        fig.text(
            0.68, 0.205,
            "This geometry plot does not establish solver convergence or "
            "solution accuracy.",
            fontsize=8.8, weight="bold", color=navy, va="top", wrap=True,
        )
        pdf.savefig(fig)
        plt.close(fig)

        fig = new_page(
            "7 | A practical student workflow",
            "Turn a metric into a defensible mesh decision",
        )
        steps = [
            ("1. Check the definition",
             "Confirm what was measured, the units, whether a low or high "
             "value is desirable, and which faces are included."),
            ("2. Read the distribution",
             "Compare median, upper percentile and maximum. Ask whether the "
             "issue is widespread or concentrated in a small tail."),
            ("3. Find the cells",
             "Use the spatial maps and cell IDs to inspect local geometry, "
             "topology, boundaries and the flow region that matters."),
            ("4. Use solver-native checks",
             "Compare these custom indicators with the quality measures and "
             "guidance for the Fluent version and discretization you will use."),
            ("5. Test numerical impact",
             "Run a mesh-sensitivity study and compare engineering outputs; "
             "geometric metrics alone cannot establish accuracy."),
        ]
        y = 0.82
        for heading, explanation in steps:
            fig.text(0.07, y, heading, fontsize=11, weight="bold",
                     color=blue, va="top")
            fig.text(0.075, y - 0.04, wrapped(explanation, 110),
                     fontsize=8.9, color=ink, va="top", linespacing=1.35)
            y -= 0.115
        fig.text(0.07, 0.235, "Important limits of this report",
                 fontsize=11.5, weight="bold", color=navy)
        limitations = (
            "The neighbour ratio and angle use interior faces only. "
            "Characteristic size is cube-root volume, not the longest edge "
            "or a measure of anisotropy. The angle is a custom centroid-to-"
            "face-normal calculation, not necessarily identical to Fluent's "
            "native orthogonal quality or skewness. This analysis does not "
            "report aspect ratio, boundary-face quality, wall y+, solver "
            "residuals, convergence or discretization error. Thresholds "
            "depend on the solver setup and application; no universal "
            "pass/fail label is assigned."
        )
        fig.text(0.075, 0.195, wrapped(limitations, 127), fontsize=8.7,
                 color=ink, va="top", linespacing=1.4)
        fig.text(
            0.07, 0.095,
            f"Source: current run's CSV and PNG figures  |  "
            f"Report date: {datetime.now():%Y-%m-%d}",
            fontsize=7.8, color=muted, va="top",
        )
        pdf.savefig(fig)
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
    write_report(args.output, centres, types, metrics, args.slab_z_mm,
                 args.slab_half_width_mm)

    print(f"{len(types):,} cells, {len(owner):,} faces, {len(nodes):,} nodes")
    print(f"{np.count_nonzero(types == 4):,} hexahedra, {np.count_nonzero(types == 7):,} polyhedra")
    for label, values in metrics:
        low, median, p95, high = np.nanpercentile(values, [0, 50, 95, 100])
        print(f"{label}: min {low:.3g}, median {median:.3g}, 95th percentile {p95:.3g}, max {high:.3g}")
    print(f"Results saved in {args.output.resolve()}")
    print(f"Report saved to {(args.output / 'mesh_quality_report.pdf').resolve()}")


if __name__ == "__main__":
    main()
