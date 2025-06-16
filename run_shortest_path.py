import subprocess
import numpy as np
import os
import plotly.graph_objects as go
import trimesh
from scipy.spatial import KDTree
import heapq
from itertools import combinations

from functions import *  # Import all functions from functions.py



def main(input_pdb, start_point, end_point, max_distance=1.9, resolution=1, output_dir=None, just_dist=False):
    # Step 1: Generate MSMS files
    mesh = run_msms(input_pdb, output_dir)
    print("Number of vertices:", mesh["num_vertexs"])
    print("Number of triangles:", mesh["num_triangles"])
    print("Number of spheres:", mesh["num_spheres"])

    # Create the mesh using trimesh
    vertices = np.array(mesh['vertices'])
    faces = np.array(mesh['faces'])
    your_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

    # Step 2: Define grid parameters
    min_bounds = vertices.min(axis=0) - 1
    max_bounds = vertices.max(axis=0) + 1

    points_outside_mesh = generate_3d_grid_outside_mesh2(your_mesh, min_bounds, max_bounds, resolution)
    
    print("Points outside mesh:", points_outside_mesh.shape)

    # Combine mesh vertices and points outside
    combined_points = np.vstack([vertices, points_outside_mesh])

    # Step 3: Build adjacency list
    adjacency_list = build_adjacency_list(combined_points, max_distance)

    print(start_point,end_point)

    # Find indices of start and end points
    start_idx = np.argmin(np.linalg.norm(combined_points - start_point, axis=1))
    end_idx = np.argmin(np.linalg.norm(combined_points - end_point, axis=1))

    euclidean_distance = np.linalg.norm(combined_points[start_idx] - combined_points[end_idx])

    # Step 4: Find shortest path using Dijkstra's algorithm
    path_indices = astar_shortest_path(start_idx, end_idx, adjacency_list, combined_points)
    shortest_path_points = combined_points[path_indices]

    # Step 5, refine the path in the forward and reversedirection
    refined_forward_path = refine_path_forward(shortest_path_points, your_mesh)
    refined_reverse_path = refine_path_reverse(refined_forward_path, your_mesh)

    #refined_reverse_path = best_refine_path(refined_reverse_path, your_mesh) # ACTIVAR ESTO CUNADO QUERAMOS EL MEJOR CAMINO

    min_dist = calculate_path_distance(refined_reverse_path)

    print(f"Final distance: {min_dist}")

    # save the final distance to a text file
    with open(os.path.join(output_dir, 'SASD.txt'), 'w') as f:
        f.write(str(min_dist)+'\n')
    with open(os.path.join(output_dir, 'euclid.txt'), 'w') as f:
        f.write(str(euclidean_distance)+'\n')

    # Step 6: Visualize the results
    if just_dist==False:
        visualize_mesh_outside_points_results(points_outside_mesh, your_mesh, refined_reverse_path)

    # Step 7: Write the PML file
    write_pml_file(refined_reverse_path, input_pdb, output_dir=output_dir)

    #if just_dist==True:
    #   rm_cmd = f"rm -r {output_dir}"
    #   subprocess.run(rm_cmd, shell=True, check=True)
        
    return min_dist, euclidean_distance

def main_no_contacts(input_pdb, start_point, end_point, max_distance=1.9, resolution=1, output_dir=None, just_dist=False):
    # Step 1: Generate MSMS files
    mesh = run_msms_separate_chains(input_pdb, output_dir)
    print("Number of vertices:", mesh["num_vertexs"])
    print("Number of triangles:", mesh["num_triangles"])

    # Create the mesh using trimesh
    vertices = np.array(mesh['vertices'])
    faces = np.array(mesh['faces'])
    your_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

    # Step 2: Define grid parameters
    min_bounds = vertices.min(axis=0) - 1
    max_bounds = vertices.max(axis=0) + 1

    points_outside_mesh = generate_3d_grid_outside_mesh2(your_mesh, min_bounds, max_bounds, resolution)
    
    print("Points outside mesh:", points_outside_mesh.shape)

    # Combine mesh vertices and points outside
    combined_points = np.vstack([vertices, points_outside_mesh])

    # Step 3: Build adjacency list
    adjacency_list = build_adjacency_list(combined_points, max_distance)

    print(start_point,end_point)

    # Find indices of start and end points
    start_idx = np.argmin(np.linalg.norm(combined_points - start_point, axis=1))
    end_idx = np.argmin(np.linalg.norm(combined_points - end_point, axis=1))

    euclidean_distance = np.linalg.norm(combined_points[start_idx] - combined_points[end_idx])

    # Step 4: Find shortest path using Dijkstra's algorithm
    path_indices = astar_shortest_path(start_idx, end_idx, adjacency_list, combined_points)
    shortest_path_points = combined_points[path_indices]

    # Step 5, refine the path in the forward and reversedirection
    refined_forward_path = refine_path_forward(shortest_path_points, your_mesh)
    refined_reverse_path = refine_path_reverse(refined_forward_path, your_mesh)

    #refined_reverse_path = best_refine_path(refined_reverse_path, your_mesh) # ACTIVAR ESTO CUNADO QUERAMOS EL MEJOR CAMINO

    min_dist = calculate_path_distance(refined_reverse_path)

    print(f"Final distance: {min_dist}")

    # save the final distance to a text file
    with open(os.path.join(output_dir, 'SASD.txt'), 'w') as f:
        f.write(str(min_dist)+'\n')
    with open(os.path.join(output_dir, 'euclid.txt'), 'w') as f:
        f.write(str(euclidean_distance)+'\n')

    # Step 6: Visualize the results
    if just_dist==False:
        visualize_mesh_outside_points_results(points_outside_mesh, your_mesh, refined_reverse_path)

    # Step 7: Write the PML file
    write_pml_file(refined_reverse_path, input_pdb, output_dir=output_dir)

    #if just_dist==True:
    #   rm_cmd = f"rm -r {output_dir}"
    #   subprocess.run(rm_cmd, shell=True, check=True)
        
    return min_dist, euclidean_distance


'''
# Example usage
input_pdb = "trimmed_AF2_outputs/P0A796_P0AB71/ranked_1.pdb"
start_point = np.array([-10.87, 2.89, -11.47])  # Example start point
end_point = np.array([0, 39, 25])          # Example end point

main(input_pdb, start_point, end_point,output_dir='prueba',just_dist=True)

'''
'''
import json

# Define the directory
directory = "../ProximityReactions/Ecoli/dock2_bueno/P00509_P00805"

# Open and read the JSON file
with open(os.path.join(directory, 'distances.json'), 'r') as f:
    data = json.load(f)

# Extract values for keys COM1 and COM2
com1_value = data.get('COM1')
com2_value = data.get('COM2')


print(f"COM1: {com1_value}, COM2: {com2_value}")

main(f"{directory}/model_1.pdb", com1_value, com2_value, output_dir='P00509_P00805')
'''

#### NEW ADJACENCY LIST WITH POINTS IN THE SURFACE THAT DO NOT GO INSIDE THE MESH
def build_adjacency_list_visibility(vertices, mesh, num_samples=5):
    """
    Build an adjacency list for the mesh surface vertices.
    Two points are connected only if the direct path between them does not enter the protein.

    Parameters:
        vertices (np.ndarray): Array of vertex coordinates (N, 3).
        mesh (trimesh.Trimesh): The mesh object representing the protein surface.
        num_samples (int): Number of intermediate points to sample along the path.

    Returns:
        dict: Adjacency list where each vertex index maps to a list of connected vertex indices.
    """
    adjacency_list = {i: [] for i in range(len(vertices))}

    for i, j in combinations(range(len(vertices)), 2):
        # Get the two points
        point_A = vertices[i]
        point_B = vertices[j]

        # Generate evenly spaced points along the segment
        segment_points = np.array([np.linspace(point_A[k], point_B[k], num_samples) for k in range(3)]).T

        # Check if any point along the path lies inside the mesh
        try:
            inside_status = mesh.contains(segment_points)
        except Exception as e:
            print(f"Error in mesh.contains for pair ({i}, {j}): {e}")
            continue

        # If no points are inside the mesh, connect the vertices
        if not any(inside_status):
            adjacency_list[i].append(j)
            adjacency_list[j].append(i)  # Ensure the graph is undirected

    return adjacency_list

def dijkstra_with_dynamic_graph(vertices, mesh, start, end, num_samples=5, distance_threshold=2.0):
    # Initialize Dijkstra structures
    queue = [(0, start)]  # (distance, vertex)
    distances = {start: 0}
    predecessors = {start: None}
    visited = set()  # To avoid processing the same vertex multiple times

    # Priority Queue for processing vertices
    while queue:
        current_dist, current_vertex = heapq.heappop(queue)

        if current_vertex == end:
            break

        if current_vertex in visited:
            continue
        visited.add(current_vertex)

        # For each vertex, dynamically build its neighbors on demand
        point_A = vertices[current_vertex]
        # Find nearby vertices within the threshold distance using KDTree (or any method)
        for j, point_B in enumerate(vertices):
            if j == current_vertex:
                continue

            distance = np.linalg.norm(point_A - point_B)

            # Check if the distance is within the threshold
            if distance <= distance_threshold:
                # Generate points between the current vertex and the neighbor
                segment_points = np.linspace(point_A, point_B, num_samples)

                # Check if any point along the path is inside the mesh
                inside_status = mesh.contains(segment_points)

                # If no points are inside the mesh, consider this as a valid neighbor
                if not any(inside_status):
                    if j not in distances or current_dist + distance < distances[j]:
                        distances[j] = current_dist + distance
                        predecessors[j] = current_vertex
                        heapq.heappush(queue, (distances[j], j))

    # Reconstruct the shortest path
    path = []
    current_vertex = end
    while current_vertex is not None:
        path.append(current_vertex)
        current_vertex = predecessors[current_vertex]
    path.reverse()

    return path, distances.get(end, float('inf'))  # Return path and the distance to the end node



def main_visibility(input_pdb, start_point, end_point, max_distance=1.9, resolution=1, output_dir=None, just_dist=False):
    # Step 1: Generate MSMS files
    mesh = run_msms(input_pdb, output_dir)
    print("Number of vertices:", mesh["num_vertexs"])
    print("Number of triangles:", mesh["num_triangles"])
    print("Number of spheres:", mesh["num_spheres"])

    # Create the mesh using trimesh
    vertices = np.array(mesh['vertices'])
    faces = np.array(mesh['faces'])
    your_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

    print("Vertices shape:", vertices.shape)
    print("Vertices sample:", vertices[:5])
    print("Mesh bounds:", your_mesh.bounds)

    # Step 3: Build adjacency list
    adjacency_list = build_adjacency_list_visibility(vertices[:100], your_mesh, num_samples=5)
    print(len(adjacency_list))
    print(start_point, end_point)

    # Find indices of start and end points
    start_idx = np.argmin(np.linalg.norm(vertices - start_point, axis=1))
    end_idx = np.argmin(np.linalg.norm(vertices - end_point, axis=1))

    start_idx=1
    end_idx=99

    # Step 4: Find shortest path using Dijkstra's algorithm
    path_indices = astar_shortest_path(start_idx, end_idx, adjacency_list, vertices)
    shortest_path_points = vertices[path_indices]
    
    distance = calculate_path_distance(shortest_path_points)
    print(f"Final distance: {distance}")

    # save the final distance to a text file
    with open(os.path.join(output_dir, 'SASD.txt'), 'w') as f:
        f.write(str(distance)+'\n')

    # Step 6: Visualize the results
    if just_dist == False:
        visualize_mesh_outside_points_results(vertices, your_mesh, shortest_path_points)

    # Step 7: Write the PML file
    write_pml_file(shortest_path_points, input_pdb, output_dir=output_dir)

    if just_dist == True:
        rm_cmd = f"rm -r {output_dir}"
        subprocess.run(rm_cmd, shell=True, check=True)
        
    return distance

'''
# SUPER COMPUTATIONALLY EXPENSIVE
# Example usage
input_pdb = "AF2_outputs/A0A0H3JGH6_P00509/ranked_1.pdb"
start_point = np.array([-10.87, 2.89, -11.47])  # Example start point
end_point = np.array([0, 39, 25])          # Example end point

main_visibility(input_pdb, start_point, end_point,output_dir='prueba')
'''



main("../ranked_1.pdb", np.array([-10.87, 2.89, -11.47]), np.array([0, 39, 25]), output_dir='prueba')

'''
# Example with no contacts
start_point = np.array([-13,30,13])  # Example start point
end_point = np.array([-5,-52,-20])          # Example end point

main_no_contacts(input_pdb="AF2_outputs/P0ABF6_P12758/ranked_1.pdb", output_dir='prueba',start_point=start_point, end_point=end_point)


'''