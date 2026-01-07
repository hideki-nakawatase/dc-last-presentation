import numpy as np


class Node:
    def __init__(self, idx, x, y):
        self.id = idx
        self.x = x
        self.y = y

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
            elm.get_coodination(self.nodes)
