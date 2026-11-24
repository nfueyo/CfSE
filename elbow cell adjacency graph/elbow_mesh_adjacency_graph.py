"""Read a Fluent .cas.h5 mesh with PyFluent and visualize cell adjacency.

The full graph is exported as a face-level CSV and an adjacency-matrix view.
A node-link view shows a selected cell and its one- or two-ring neighbours.
"""
from __future__ import annotations

import argparse
import csv
from collections import deque
from pathlib import Path

import h5py
import numpy as np


def read_adjacency(case_path: Path):
    """Return cell types and unique undirected cell pairs from Fluent HDF5."""
    with h5py.File(case_path, "r") as case:
        mesh = case["meshes/1"]
        cell_types = mesh["cells/ctype/1/cell-types"][:]
        owner = mesh["faces/c0/1"][:].astype(np.int64) - 1
        neighbour_dataset = mesh["faces/c1/1"]
        neighbour_start = int(neighbour_dataset.attrs["minId"][0]) - 1
        neighbour_stop = int(neighbour_dataset.attrs["maxId"][0])
        neighbour_values = neighbour_dataset[:].astype(np.int64) - 1
        if (
            neighbour_start < 0
            or neighbour_stop > len(owner)
            or neighbour_start > neighbour_stop
            or len(neighbour_values) != neighbour_stop - neighbour_start
        ):
            raise ValueError("Invalid neighbour-face range in the Fluent case.")
        neighbours = np.full(len(owner), -1, dtype=np.int64)
        neighbours[neighbour_start:neighbour_stop] = neighbour_values

    if (
        len(owner) != len(neighbours)
        or len(cell_types) == 0
        or np.any(owner < 0)
        or np.any(owner >= len(cell_types))
        or np.any(neighbours < -1)
        or np.any(neighbours >= len(cell_types))
        or np.any(owner == neighbours)
    ):
        raise ValueError("Invalid connectivity for the supported Fluent HDF5 layout.")

    interior = neighbours >= 0
    face_ids = np.flatnonzero(interior) + 1
    face_pairs = np.column_stack((owner[interior], neighbours[interior]))
    unique_pairs = np.unique(np.sort(face_pairs, axis=1), axis=0)
    if not len(unique_pairs):
        raise ValueError("The case contains no interior faces for an adjacency graph.")
    return cell_types, face_ids, face_pairs, unique_pairs


def build_adjacency(cell_count, unique_pairs):
    adjacency = [set() for _ in range(cell_count)]
    for cell_a, cell_b in unique_pairs:
        adjacency[cell_a].add(int(cell_b))
        adjacency[cell_b].add(int(cell_a))
    return adjacency


def neighborhood(adjacency, seed, rings):
    distances = {seed: 0}
    pending = deque([seed])
    while pending:
        cell = pending.popleft()
        if distances[cell] == rings:
            continue
        for neighbour in sorted(adjacency[cell]):
            if neighbour not in distances:
                distances[neighbour] = distances[cell] + 1
                pending.append(neighbour)
    return distances


def write_edge_list(output, face_ids, face_pairs):
    with (output / "cell_adjacency.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(["face_id", "cell_id_a", "cell_id_b"])
        for face_id, (cell_a, cell_b) in zip(face_ids, face_pairs):
            writer.writerow([face_id, cell_a + 1, cell_b + 1])


def plot_graph(output, unique_pairs, adjacency, cell_types, seed, rings, show):
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    node_distances = neighborhood(adjacency, seed, rings)
    if len(node_distances) > 250:
        raise ValueError(
            "The selected neighborhood has more than 250 cells; use --rings 1."
        )
    nodes = sorted(node_distances, key=lambda cell: (node_distances[cell], cell))
    selected = set(nodes)
    edges = [
        (int(cell_a), int(cell_b))
        for cell_a, cell_b in unique_pairs
        if cell_a in selected and cell_b in selected
    ]

    positions = {seed: (0.0, 0.0)}
    for ring in range(1, rings + 1):
        ring_nodes = [cell for cell in nodes if node_distances[cell] == ring]
        for index, cell in enumerate(ring_nodes):
            angle = 2 * np.pi * index / max(1, len(ring_nodes))
            positions[cell] = (ring * np.cos(angle), ring * np.sin(angle))

    figure, (matrix_axis, local_axis) = plt.subplots(
        1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [1.1, 1]}
    )
    matrix_axis.scatter(
        np.concatenate((unique_pairs[:, 0] + 1, unique_pairs[:, 1] + 1)),
        np.concatenate((unique_pairs[:, 1] + 1, unique_pairs[:, 0] + 1)),
        s=0.35, marker=".", color="#64748b", alpha=0.6,
    )
    matrix_axis.set(
        xlabel="Cell ID", ylabel="Adjacent cell ID",
        title=f"Full adjacency matrix ({len(adjacency):,} cells)",
        xlim=(1, len(adjacency)), ylim=(1, len(adjacency)),
    )
    matrix_axis.set_aspect("equal", adjustable="box")
    matrix_axis.grid(False)

    segments = [(positions[a], positions[b]) for a, b in edges]
    if segments:
        local_axis.add_collection(
            LineCollection(segments, colors="#94a3b8", linewidths=1.2, zorder=1)
        )
    colors = {4: "#168a78", 7: "#3975a8"}
    for cell in nodes:
        x, y = positions[cell]
        is_seed = cell == seed
        local_axis.scatter(
            [x], [y], s=170 if is_seed else 95,
            c="#d9473f" if is_seed else colors.get(int(cell_types[cell]), "#7c3aed"),
            edgecolors="white", linewidths=0.8, zorder=3,
        )
        if node_distances[cell] < rings or is_seed or len(nodes) <= 60:
            local_axis.annotate(
                str(cell + 1), (x, y), xytext=(0, 7),
                textcoords="offset points", ha="center", fontsize=7.5,
            )
    local_axis.set_aspect("equal")
    local_axis.set_xlim(-rings - 0.8, rings + 0.8)
    local_axis.set_ylim(-rings - 0.8, rings + 0.8)
    local_axis.axis("off")
    local_axis.set_title(
        f"Cell {seed + 1} neighborhood: {rings} ring(s), "
        f"{len(nodes):,} cells, {len(edges):,} edges"
    )
    legend = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#d9473f",
               markersize=8, label=f"Selected cell {seed + 1}"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[4],
               markersize=7, label="Hexahedron (type 4)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[7],
               markersize=7, label="Polyhedron (type 7)"),
    ]
    local_axis.legend(handles=legend, loc="lower center",
                      bbox_to_anchor=(0.5, 0.02), ncol=1, frameon=False)
    figure.suptitle("Fluent cell-face adjacency", fontsize=15, weight="bold")
    figure.tight_layout(rect=(0, 0.03, 1, 0.94))
    figure.savefig(output / "adjacency_graph.png", dpi=180)
    if show:
        plt.show()
    plt.close(figure)
    return len(nodes), len(edges)


def load_case_with_pyfluent(case_path, output, processors, version):
    import ansys.fluent.core as pyfluent

    kwargs = {
        "mode": "solver",
        "precision": "double",
        "dimension": 3,
        "processor_count": processors,
        "cwd": str(output.resolve()),
    }
    if version:
        kwargs["product_version"] = version

    solver = pyfluent.launch_fluent(**kwargs)
    try:
        solver.settings.file.read_case(file_name=str(case_path.resolve()))
        return solver.get_fluent_version(), pyfluent.__version__
    finally:
        solver.exit()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True,
                        help="Fluent HDF5 case file (.cas.h5)")
    parser.add_argument("--output", type=Path, default=Path("adjacency_graph_results"))
    parser.add_argument("--cell-id", type=int, default=1,
                        help="One-based Fluent cell ID for the local graph (default: 1)")
    parser.add_argument("--rings", type=int, choices=(1, 2), default=1,
                        help="Neighbourhood depth around --cell-id (default: 1)")
    parser.add_argument("--processors", type=int, default=1)
    parser.add_argument("--version", help="Fluent version, e.g. 26.1")
    parser.add_argument("--show", action="store_true",
                        help="Open the graph window after saving the PNG")
    args = parser.parse_args(argv)

    if not args.case.is_file():
        parser.error(f"Case file does not exist: {args.case}")
    if not args.case.name.lower().endswith(".cas.h5"):
        parser.error("This reader requires a Fluent HDF5 case ending in .cas.h5.")
    if args.cell_id < 1 or args.processors < 1:
        parser.error("--cell-id and --processors must be positive.")

    args.output.mkdir(parents=True, exist_ok=True)
    fluent_version, pyfluent_version = load_case_with_pyfluent(
        args.case, args.output, args.processors, args.version
    )
    cell_types, face_ids, face_pairs, unique_pairs = read_adjacency(args.case)
    if args.cell_id > len(cell_types):
        parser.error(
            f"--cell-id must be between 1 and {len(cell_types):,}."
        )
    adjacency = build_adjacency(len(cell_types), unique_pairs)
    write_edge_list(args.output, face_ids, face_pairs)
    shown_cells, shown_edges = plot_graph(
        args.output, unique_pairs, adjacency, cell_types,
        args.cell_id - 1, args.rings, args.show,
    )

    degrees = np.fromiter((len(neighbours) for neighbours in adjacency), dtype=int)
    print(f"Loaded {args.case.resolve()} with Fluent {fluent_version} "
          f"(PyFluent {pyfluent_version})")
    print(f"{len(cell_types):,} cells; {len(face_ids):,} interior faces; "
          f"{len(unique_pairs):,} unique adjacency edges")
    print(f"Cell degree: min {degrees.min()}, median {np.median(degrees):g}, "
          f"max {degrees.max()}")
    print(f"Cell {args.cell_id} view: {shown_cells:,} cells and "
          f"{shown_edges:,} edges")
    print(f"Graph saved to {(args.output / 'adjacency_graph.png').resolve()}")
    print(f"Edge list saved to {(args.output / 'cell_adjacency.csv').resolve()}")


if __name__ == "__main__":
    main()
