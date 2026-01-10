import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from mesh import Node, Element, Mesh
from calc_matrix import Material, Solver


def create_mesh(length, height, mesh_size, thickness):
    # Node作成
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

        # Element作成
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
    return nodes_dict, elems


def plot_mesh(mesh):
    fig, ax = plt.subplots()
    for elem in mesh.elements:
        patch = patches.Polygon(xy=elem.xy, ec="black")
        ax.add_patch(patch)
        text_xy = np.mean(elem.xy, axis=0)
        ax.text(text_xy[0], text_xy[1], elem.id, fontsize=12, va="center", ha="center")
    for node in mesh.nodes.values():
        ax.scatter(node.x, node.y, fc="black", s=100)
        ax.text(
            node.x, node.y, node.id, fontsize=8, color="white", va="center", ha="center"
        )
    ax.autoscale()
    ax.set_aspect("equal", "box")
    plt.show()


# モデル、要素サイズ定義
# MESH_SIZE = 5.0
# LENGTH = 10.0
# HEIGHT = 10.0
# THICKNESS = 1.0

# nodes_dict, elems = create_mesh(LENGTH, HEIGHT, MESH_SIZE, THICKNESS)
# mesh = Mesh(nodes_dict, elems)
# for node in mesh.nodes.values():
#     print(node)
# for elem in mesh.elements:
#     print(elem)

# plot_mesh(mesh)
