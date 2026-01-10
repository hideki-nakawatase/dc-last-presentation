from mesh import Node, Element, Mesh
from calc_matrix import Material, Solver
from mesh_plot import create_mesh

MESH_SIZE = 5
LENGTH = 30.0
HEIGHT = 30.0
THICKNESS = 1.0
E = 208000
NU = 0.27

nodes_dict, elems = create_mesh(LENGTH, HEIGHT, MESH_SIZE, THICKNESS)
material = Material(E, NU)
mesh = Mesh(nodes_dict, elems)

num_mesh_len = int(LENGTH / MESH_SIZE)
num_mesh_hei = int(HEIGHT / MESH_SIZE)
if num_mesh_len == 0:
    num_mesh_len = 1
if num_mesh_hei == 0:
    num_mesh_hei = 1

for k, v in mesh.nodes.items():
    if v.x == 0.0:
        v.bc_dist[0] = 0.0
    if v.x == 0.0:
        v.bc_dist[1] = 0.0
    if v.x == LENGTH and (v.y == 0.0 or v.y == HEIGHT):
        v.bc_force = [10000 / num_mesh_len / 2, 0.0]
    elif v.x == LENGTH:
        v.bc_force = [10000 / num_mesh_len, 0.0]

solver = Solver(mesh, material)
solver.run()
solver.plot_stress()
