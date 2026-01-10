import numpy as np


class Node:
    def __init__(self, idx, x, y):
        self.id = idx
        self.x = x
        self.y = y
        self.bc_dist = [None, None]
        self.bc_force = [0, 0]

    def __str__(self):
        return "Node {}: x{}, y{}".format(self.id, self.x, self.y)


class Element:
    def __init__(self, idx, node_idxs, thickness):
        self.id = idx
        self.node = node_idxs
        self.thickness = thickness
        self.xy = None

    def __str__(self):
        return "Element {}: Node {}, Thickness: {}".format(
            self.id, self.node, self.thickness
        )

    def get_coordination(self, node_dict):
        res = []
        for node_idx in self.node:
            res.append([node_dict[node_idx].x, node_dict[node_idx].y])
        self.xy = np.array(res)


class Mesh:
    def __init__(self, nodes_dict, elements):
        self.nodes = nodes_dict
        self.elements = elements
        self.get_element_coord()

    def get_element_coord(self):
        for elm in self.elements:
            elm.get_coordination(self.nodes)


mesh_size = 5.0
length = 10.0
height = 10.0

# node作成
num_mesh_len = int(length / mesh_size)
num_mesh_hei = int(height / mesh_size)
if num_mesh_len == 0:
    num_mesh_len = 1
if num_mesh_hei == 0:
    num_mesh_hei = 1
x = np.linspace(0, length, num_mesh_len + 1)
y = np.linspace(0, height, num_mesh_hei + 1)
X, Y = np.meshgrid(x, y)
X = X.ravel()
Y = Y.ravel()

nodes_dict = {}
for i, coord in enumerate(zip(X, Y)):
    nodes_dict[i] = Node(i, coord[0], coord[1])

for node in nodes_dict.values():
    print(node)

thickness = 1

# element作成
node_idx = 0
elem_idx = 0
elems = []
for i in range(num_mesh_hei):
    for j in range(num_mesh_len + 1):
        if (node_idx + 1) % (num_mesh_len + 1) == 0:
            node_idx += 1
            continue
        else:
            node_idxs = [
                node_idx,
                node_idx + 1,
                node_idx + num_mesh_len + 2,
                node_idx + num_mesh_len + 1,
            ]
            elems.append(Element(elem_idx, node_idxs, thickness))
            node_idx += 1
            elem_idx += 1

for elem in elems:
    print(elem)
