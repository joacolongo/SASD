# Solvent Accessible Surface Distance Calculator

    Joaquin Algorta - Max Planck Institute of Molecular Plant Physiology

This tool computes the shortest accessible path between two points around a molecular surface, avoiding the interior of the mesh. It is designed for structural biology applications, such as analyzing protein-protein interfaces or ligand accessibility.

## Features
Surface Mesh Generation: Uses MSMS to generate a molecular surface mesh from a PDB file.
Grid Sampling: Samples points outside the mesh to allow pathfinding around the surface.
Shortest Path Calculation: Computes the shortest path between two user-defined points using A* search, avoiding the mesh interior.
Path Refinement: Refines the path to minimize unnecessary detours.
Visualization: Optionally visualizes the mesh, sampled points, and computed path (requires Plotly).
PyMOL Export: Generates a PyMOL script to visualize the path and structure.
No-Contact Mode: Optionally computes the path with separated chains.

## Usage 
    python run_shortest_path.py <input_pdb> <start_x> <start_y> <start_z> <end_x> <end_y> <end_z> --output_dir <output_folder>

Example
    python run_shortest_path.py ../ranked_1.pdb -10.87 2.89 -11.47 0 39 25 --output_dir prueba

## Arguments
Mandatory positional arguments:
- `<input_pdb>`:  
  Path to the input PDB file for which the surface and path will be calculated.
- `<start_x> <start_y> <start_z>`:  
  The X, Y, and Z coordinates of the starting point (in the same coordinate system as the PDB).
- `<end_x> <end_y> <end_z>`:  
  The X, Y, and Z coordinates of the ending point.
  
Optional Arguments
* --output_dir: Directory to save outputs (default: output)
* --max_distance: Maximum distance to connect nodes in the graph (default: 1.9)
* --resolution: Grid resolution for sampling outside points (default: 1.0)
* --just_dist: Only compute and save distances, skip visualization
* --no_contacts: Use no-contact mode (separate chains)

Optional Arguments
* --output_dir: Directory to save outputs (default: output)
* --max_distance: Maximum distance to connect nodes in the graph (default: 1.9)
* --resolution: Grid resolution for sampling outside points (default: 1.0)
* --just_dist: Only compute and save distances, skip visualization
* --no_contacts: Use no-contact mode (separate chains)

Output
* SASD.txt: Shortest accessible surface distance
* euclid.txt: Euclidean distance between the two points
* SASD.pml: PyMOL script for visualization


## Requirements
* Python 3.x
* trimesh
* plotly (for visualization)
* MSMS and pdb_to_xyzr binaries in your PATH
