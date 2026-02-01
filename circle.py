import numpy as np
import matplotlib.pyplot as plt

plt.switch_backend("agg")


def create_circle_mesh(radius=1, center=[0, 0], layers=20, n_points=30, max_dist=5):
    theta = np.linspace(0, 2 * np.pi, n_points, endpoint=True)

    box_w = 6.0
    box_h = 6.0

    x_mesh = []
    y_mesh = []

    for i in range(layers):
        w = i / (layers - 1)

        offset = i * max_dist / layers
        inner_x = center[0] + (radius + offset) * np.cos(theta)
        inner_y = center[1] + (radius + offset) * np.sin(theta)

        scale = np.minimum(
            box_w / (np.abs(np.cos(theta)) + 1e-10),
            box_h / (np.abs(np.sin(theta)) + 1e-10),
        )
        outer_x = center[0] + scale * np.cos(theta)
        outer_y = center[1] + scale * np.sin(theta)

        x_layer = (1 - w) * inner_x + (w) * outer_x
        y_layer = (1 - w) * inner_y + (w) * outer_y

        x_mesh.append(x_layer)
        y_mesh.append(y_layer)

    return np.array(x_mesh), np.array(y_mesh)


X, Y = create_circle_mesh()
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

filename = "circle_data_u=0/circle_mesh_data.dat"

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

bc_wall = []
bc_in = []  
bc_out = []  
bc_side = []  

for i in range(n_layers):
    for j in range(n_points):
        node_id = i * n_points + j
        x_val = nodes[node_id, 1]
        if i == 0:
            bc_wall.append(node_id)
        elif i == n_layers - 1:
            if x_val < -5.0:
                bc_in.append(node_id)
            elif x_val > 5.0:
                bc_out.append(node_id)
            else:
                bc_side.append(node_id)

with open(filename, "a") as f:
    f.write("BOUNDARY_CONDITIONS\n")
    f.write(f"WALL {' '.join(map(str, bc_wall))}\n")
    f.write(f"INLET {' '.join(map(str, bc_in))}\n")
    f.write(f"OUTLET {' '.join(map(str, bc_out))}\n")
    f.write(f"SIDE {' '.join(map(str, bc_side))}\n")

print(f"Successfully updated with endpoint=True. Total Nodes: {len(nodes)}")

u = np.ones(len(nodes))
u[bc_wall] = 0
v = np.zeros(len(nodes))
p = np.zeros(len(nodes))


for node_id in bc_wall:
    u[node_id] = 0.0
    v[node_id] = 0.0

for node_id in bc_in:
    u[node_id] = 1.0
    v[node_id] = 0.0


element_centers = []
element_areas = []

for elem in elements:
    points = nodes[elem[1:].astype(int), 1:]
    center = np.mean(points, axis=0)
    element_centers.append(center)

    area = 0.5 * abs(
        (points[0, 0] - points[2, 0]) * (points[1, 1] - points[3, 1])
        - (points[1, 0] - points[3, 0]) * (points[0, 1] - points[2, 1])
    )
    element_areas.append(area)

element_centers = np.array(element_centers)
element_areas = np.array(element_areas)

fig, ax = plt.subplots(figsize=(8, 8))
for elem in elements:
    idx = elem[1:].astype(int)
    idx_closed = np.append(idx, idx[0])
    ax.plot(nodes[idx_closed, 1], nodes[idx_closed, 2], "k-", lw=0.3, alpha=0.3)

step = 1
ax.quiver(
    nodes[::step, 1],
    nodes[::step, 2],
    u[::step],
    v[::step],
    color="blue",
    scale=20,
    width=0.005,
    pivot="mid",
    label="Velocity (Initial)",
)

ax.scatter(nodes[bc_wall, 1], nodes[bc_wall, 2], c="red", s=10, label="Wall (Fixed 0)")
ax.scatter(nodes[bc_in, 1], nodes[bc_in, 2], c="green", s=10, label="In (Fixed 1.0)")
ax.scatter(nodes[bc_out, 1], nodes[bc_out, 2], c="orange", s=10, label="Out")

ax.set_aspect("equal")
ax.set_title("Initial Velocity Field & Boundary Condition Check")
ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.0))
plt.tight_layout()
plt.savefig("circle_data_u=1/circle_initial.png")

nu = 0.5  
dt = 0.5  
n_step = 5000  
rho = 1.0  
p_iterations = 20
omega = 0.1
friction_coeff = 0.2

p = np.zeros(len(nodes))

for step in range(n_step):
    u_old = u.copy()
    v_old = v.copy()

    u_star = u.copy()
    v_star = v.copy()

    for i in range(1, n_layers - 1):
        for j in range(n_points):
            idx = i * n_points + j
            up, down = (i + 1) * n_points + j, (i - 1) * n_points + j
            left = i * n_points + (j - 1 if j > 0 else n_points - 2)
            right = i * n_points + (j + 1 if j < n_points - 1 else 1)

            # if u_old[idx] > 0:
            #     df_dx_u = u_old[idx] - u_old[left]
            #     df_dx_v = v_old[idx] - v_old[left]
            # else:
            #     df_dx_u = u_old[right] - u_old[idx]
            #     df_dx_v = v_old[right] - v_old[idx]

            # if v_old[idx] > 0:
            #     df_dy_u = u_old[idx] - u_old[down]
            #     df_dy_v = v_old[idx] - v_old[down]
            # else:
            #     df_dy_u = u_old[up] - u_old[idx]
            #     df_dy_v = v_old[up] - v_old[idx]

            # # 譛邨ら噪縺ｪ蟇ｾ豬���
            # adv_u = u_old[idx] * df_dx_u + v_old[idx] * df_dy_u
            # adv_v = u_old[idx] * df_dx_v + v_old[idx] * df_dy_v

            # adv_u = 0
            # adv_v = 0

            lap_u = (
                u_old[up] + u_old[down] + u_old[left] + u_old[right] - 4 * u_old[idx]
            )
            lap_v = (
                v_old[up] + v_old[down] + v_old[left] + v_old[right] - 4 * v_old[idx]
            )

            u_star[idx] = u_old[idx] + dt * (nu * lap_u)
            v_star[idx] = v_old[idx] + dt * (nu * lap_v)
    for _ in range(p_iterations):
        p_new = p.copy()
        for i in range(1, n_layers - 1):
            for j in range(n_points):
                idx = i * n_points + j
                up, down = (i + 1) * n_points + j, (i - 1) * n_points + j
                left = i * n_points + (j - 1 if j > 0 else n_points - 2)
                right = i * n_points + (j + 1 if j < n_points - 1 else 1)

                b = (rho / dt) * (
                    (u_star[right] - u_star[left]) + (v_star[up] - v_star[down])
                )

                p_update = 0.25 * (p[up] + p[down] + p[left] + p[right] - b)
                p_new[idx] = (1 - omega) * p[idx] + omega * p_update
                p_new[idx] = 0

        p = p_new.copy()
        p[np.array(bc_in)] = p[np.array(bc_in) - n_points]
        p[bc_out] = 0.0

    for i in range(1, n_layers - 1):
        for j in range(n_points):
            idx = i * n_points + j
            up, down = (i + 1) * n_points + j, (i - 1) * n_points + j
            left = i * n_points + (j - 1 if j > 0 else n_points - 2)
            right = i * n_points + (j + 1 if j < n_points - 1 else 1)

            u[idx] = u_star[idx] - (dt / rho) * (p[right] - p[left])
            v[idx] = v_star[idx] - (dt / rho) * (p[up] - p[down])
    near_wall_idx = np.array(bc_wall) + n_points
    u[bc_wall] = u[near_wall_idx] * (1.0 - friction_coeff)
    v[bc_wall] = v[near_wall_idx] * (1.0 - friction_coeff)
    u[bc_in], v[bc_in] = 1.0, 0.0

    if step % 100 == 0:
        print(f"Step {step} finished. Pressure and Velocity updated.")
        plt.figure(figsize=(10, 8))
        speed = np.sqrt(u**2 + v**2)
        q = plt.quiver(
            nodes[:, 1],
            nodes[:, 2],
            u,
            v,
            speed,
            cmap="viridis",
            angles="xy",
            scale_units="xy",
            scale=0.5,
            width=0.003,
            minshaft=2,
        )
        plt.scatter(
            nodes[bc_wall, 1], nodes[bc_wall, 2], c="red", s=5, label="Cylinder Wall"
        )

        plt.colorbar(q, label="Flow Speed")
        plt.title(f"Velocity Field after {n_step} steps (FDM on O-grid)")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.axis("equal")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.savefig(f"circle_data_u=1/friction_az_{step}.png")