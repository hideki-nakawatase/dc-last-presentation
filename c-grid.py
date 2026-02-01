import numpy as np
import matplotlib.pyplot as plt


def create_c_grid(radius, n_surface=60, n_layers=10, wake_length=40, x_shift=30):
    theta = np.linspace(np.pi / 2, 3 * np.pi / 2, int(n_surface / 2))
    x_wake = np.linspace(0, wake_length, int(n_surface / 2))

    x_front_new = radius * np.cos(theta) + x_shift
    y_front_new = radius * np.sin(theta)

    y_upper_new = -radius * (x_wake[1:] / wake_length) ** 2 + radius
    y_lower_new = radius * (x_wake[1:] / wake_length) ** 2 - radius
    x_wake_actual = x_wake[1:] + x_shift

    x_line = np.concatenate([x_wake_actual[::-1], x_front_new[::-1], x_wake_actual])
    y_line = np.concatenate([y_lower_new[::-1], y_front_new[::-1], y_upper_new])

    dx = np.gradient(x_line)
    dy = np.gradient(y_line)
    nx = -dy
    ny = dx
    mag = np.sqrt(nx**2 + ny**2)
    nx, ny = nx / mag, ny / mag

    is_wake = x_line >= x_shift

    box_x_min = 0
    box_x_max = x_shift + wake_length + 20
    box_y_max = 60

    grid_x = []
    grid_y = []
    for i in range(n_layers):
        w = i / (n_layers - 1)

        dist = (1.3**i - 1) * 2.0
        inner_x = x_line + nx * dist
        inner_y = y_line + ny * dist

        inner_x[is_wake] = x_line[is_wake]
        inner_y[is_wake] = y_line[is_wake] + np.sign(y_line[is_wake]) * dist

        outer_x = np.clip(inner_x, box_x_min, box_x_max)
        outer_y = np.clip(inner_y, -box_y_max, box_y_max)

        x_layer = (1 - w**2) * inner_x + w**2 * outer_x
        y_layer = (1 - w**2) * inner_y + w**2 * outer_y

        x_back = np.linspace(x_shift + wake_length, box_x_max, 10)
        y_back_upper = np.ones_like(x_back) * y_layer[-2]
        y_back_lower = np.ones_like(x_back) * (-y_layer[-2])

        x = np.concatenate([x_layer, x_back, x_back[::-1]])
        y = np.concatenate([y_layer, y_back_upper, y_back_lower])

        grid_x.append(x)
        grid_y.append(y)

    return np.array(grid_x), np.array(grid_y)


X, Y = create_c_grid(10)
n_layers, n_points = X.shape

nodes = []
node_id = 0
for i in range(n_layers):
    for j in range(n_points):
        nodes.append([node_id, X[i, j], Y[i, j]])
        node_id += 1
nodes = np.array(nodes)

elements = []
elem_id = 0
for i in range(n_layers - 1):
    for j in range(n_points - 1):
        n1 = i * n_points + j
        n2 = i * n_points + j + 1
        n3 = (i + 1) * n_points + j + 1
        n4 = (i + 1) * n_points + j

        elements.append([elem_id, n1, n2, n3, n4])
        elem_id += 1
elements = np.array(elements)

filename = "c_grid_data/c_grid_mesh_data.dat"

with open(filename, "w") as f:
    f.write(f"# Mesh Data: Nodes={len(nodes)}, Elements={len(elements)}\n")

    f.write("NODES\n")
    for node in nodes:
        f.write(f"{int(node[0])} {node[1]:.8f} {node[2]:.8f}\n")

    f.write("ELEMENTS\n")
    for elem in elements:
        f.write(
            f"{int(elem[0])} {int(elem[1])} {int(elem[2])} {int(elem[3])} {int(elem[4])}\n"
        )

print(f"Successfully saved to {filename}")