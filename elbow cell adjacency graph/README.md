# Elbow cell adjacency graph

`elbow_mesh_adjacency_graph.py` reads a Fluent `.cas.h5` case through PyFluent
and builds a cell adjacency graph from the case's face connectivity. Cells are
graph nodes; an undirected edge joins two cells that share an interior face.
The script does not run solver iterations or modify the case.

## Requirements

The Python environment needs `ansys-fluent-core`, `h5py`, NumPy, and
Matplotlib, as well as a working Fluent installation and license:

```powershell
python -m pip install ansys-fluent-core h5py numpy matplotlib
```

## Run

Open a terminal in this directory and run:

```powershell
python .\elbow_mesh_adjacency_graph.py `
  --case "C:\path\to\elbow_no_h2.msh.cas.h5" `
  --cell-id 5364 --rings 1 --show
```

`--cell-id` is a one-based Fluent cell ID (default: `1`). `--rings` may be `1`
or `2` (default: `1`); it sets how far the local graph expands from that cell.
The selected cell is the starting point (distance 0). A cell is one ring away
if it shares an interior face with the starting cell; a cell is two rings away
if the shortest path through face-sharing cells is two edges. The local view
includes the starting cell and all cells up to the requested distance. Thus,
`--rings 1` shows the selected cell and its immediate neighbors, while
`--rings 2` also includes their previously unseen neighbors. The rings describe
topological distance, not physical distance or mesh layers.

`--show` opens the saved graph in a Matplotlib window. Use `--help` for the
Fluent version, processor-count, and output-directory options.

## Outputs

By default, files are written to `adjacency_graph_results` under the current
terminal directory. The supplied graph preview and edge list are already in
that folder:

- `adjacency_graph.png` shows the full adjacency matrix and the selected cell's
  local node-link graph. The circular arrangement in the local view represents
  graph rings, not physical cell coordinates.
- `cell_adjacency.csv` lists interior faces with their face ID and the two
  adjacent one-based cell IDs. The plotted matrix collapses repeated pairs
  into one undirected edge.

The graph is read from the `.cas.h5` face-owner (`c0`) and face-neighbour
(`c1`) connectivity. Boundary faces have no second cell and are not graph
edges. The reader supports the single-section Fluent HDF5 layout used by the
supplied case.
