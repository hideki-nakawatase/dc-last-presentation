import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import tri


class Material:
    def __init__(self, E, nu):
        self.E = E
        self.nu = nu


class Solver:
    def __init__(self, mesh, material):
        self.mesh = mesh
        self.material = material
        self.D = self.calc_D()

    def run(self):
        self.f, self.bc_dist_dict = self.read_boundary_cond()
        self.K = self.calc_K()
        self.d = self.solve()
        self.stress_dict = self.calc_stress()

    def calc_D(self):
        D = (
            np.array(
                [
                    [1, self.material.nu, 0],
                    [self.material.nu, 1, 0],
                    [0, 0, (1 - self.material.nu) / 2],
                ]
            )
            * self.material.E
            / (1 - self.material.nu**2)
        )
        return D

    def calc_B(self, elem, xi, eta):
        print(f"Element {elem.id} has {len(elem.xy)} nodes.")
        dndxi = np.array([-1 + eta, 1 - eta, 1 + eta, -1 - eta]) / 4
        dndeta = np.array([-1 + xi, -1 - xi, 1 + xi, 1 - xi]) / 4
        x = elem.xy[:, 0]
        y = elem.xy[:, 1]
        dxdxi = np.dot(dndxi, x)
        dydxi = np.dot(dndxi, y)
        dxdeta = np.dot(dndeta, x)
        dydeta = np.dot(dndeta, y)
        J = np.array([[dxdxi, dydxi], [dxdeta, dydeta]])
        B = np.zeros((3, 8))
        for i in range(4):
            Bi = np.dot(np.linalg.inv(J), np.array([dndxi[i], dndeta[i]]).T)
            B[0, 2 * i] = Bi[0]
            B[1, 2 * i + 1] = Bi[1]
            B[2, 2 * i] = Bi[1]
            B[2, 2 * i + 1] = Bi[0]
        return B, J

    def calc_Ke(self, elem):
        gps = [
            [-1 / np.sqrt(3), -1 / np.sqrt(3), 1, 1],
            [1 / np.sqrt(3), -1 / np.sqrt(3), 1, 1],
            [1 / np.sqrt(3), 1 / np.sqrt(3), 1, 1],
            [-1 / np.sqrt(3), 1 / np.sqrt(3), 1, 1],
        ]
        Ke = np.zeros((8, 8))
        for xi, eta, wi, wj in gps:
            B, J = self.calc_B(elem, xi, eta)
            Ke += (
                wi
                * wj
                * np.dot(B.T, np.dot(self.D, B))
                * np.linalg.det(J)
                * elem.thickness
            )
        return Ke

    def calc_K(self):
        K = np.zeros((len(self.mesh.nodes) * 2, len(self.mesh.nodes) * 2))
        for elem in self.mesh.elements:
            Ke = self.calc_Ke(elem)
            for i in range(4):
                for j in range(4):
                    K[
                        2 * elem.node[i] : 2 * elem.node[i] + 2,
                        2 * elem.node[j] : 2 * elem.node[j] + 2,
                    ] += Ke[2 * i : 2 * i + 2, 2 * j : 2 * j + 2]
        return K

    def read_boundary_cond(self):
        bc_dist_dict = {}
        f = []
        for k, v in self.mesh.nodes.items():
            f += v.bc_force
            if v.bc_dist[0] is not None:
                bc_dist_dict[k * 2] = v.bc_dist[0]
            if v.bc_dist[1] is not None:
                bc_dist_dict[k * 2 + 1] = v.bc_dist[1]
        return np.array(f), bc_dist_dict

    def solve(self):
        for i in range(len(self.mesh.nodes) * 2):
            if i in self.bc_dist_dict.keys():
                self.K[i, :] = np.zeros(len(self.mesh.nodes) * 2)
                self.K[:, i] = np.zeros(len(self.mesh.nodes) * 2)
                self.K[i, i] = 1.0
        d = np.dot(np.linalg.inv(self.K), self.f.T)
        return d

    def plot_deform(self, ratio=50):
        x = np.array([[node.x, node.y] for idx, node in self.mesh.nodes.items()])
        x_new = x + self.d.reshape(len(self.mesh.nodes), 2) * ratio
        fig, ax = plt.subplots()
        for v in self.mesh.elements:
            patch = patches.Polygon(xy=v.xy, ec="black", alpha=0.3, fill=False)
            ax.add_patch(patch)
        for v in self.mesh.elements:
            xy_new = [(x_new[idx, 0], x_new[idx, 1]) for idx in v.node]
            patch = patches.Polygon(
                xy=xy_new, fc="red", ec="red", alpha=0.3, fill=False
            )
            ax.add_patch(patch)
        ax.autoscale()
        ax.set_aspect("equal", "box")
        plt.show()

    def calc_stress(self):
        stress_dict = {i: np.zeros(3) for i in range(len(self.mesh.nodes))}
        counts_dict = {i: 0 for i in range(len(self.mesh.nodes))}

        gps = (
            (-1 / np.sqrt(3), -1 / np.sqrt(3), 1, 1),
            (1 / np.sqrt(3), -1 / np.sqrt(3), 1, 1),
            (1 / np.sqrt(3), 1 / np.sqrt(3), 1, 1),
            (-1 / np.sqrt(3), 1 / np.sqrt(3), 1, 1),
        )

        for elem in self.mesh.elements:
            node_ind = []
            for ind in elem.node:
                node_ind.append(ind * 2)
                node_ind.append(ind * 2 + 1)
                d_elem = self.d[node_ind]

            stress_e = {}
            ave_stress_e = np.zeros(3)
            for (
                i,
                gp,
            ) in enumerate(gps):
                xi, eta = gp[0], gp[1]
                B, J = self.calc_B(elem, xi, eta)
                stress_e[i] = np.dot(self.D, np.dot(B, d_elem).T)
                ave_stress_e += stress_e[i]
            ave_stress_e /= 4

            for i in range(4):
                stress_dict[elem.node[i]] += (stress_e[i] - ave_stress_e) * np.sqrt(
                    3
                ) + ave_stress_e
                counts_dict[elem.node[i]] += 1

        for i in range(len(self.mesh.nodes)):
            stress_dict[i] = stress_dict[i] / counts_dict[i]
            res = stress_dict[i]
            misses = np.sqrt(
                (res[0] - res[1]) ** 2 + res[0] ** 2 + res[1] ** 2 + 6 * res[2] ** 2
            ) / np.sqrt(2)
            stress_dict[i] = np.append(stress_dict[i], misses)
        return stress_dict

    def plot_stress(self, stress_idx=0, ratio=50):
        x = np.array([[node.x, node.y] for idx, node in self.mesh.nodes.items()])
        x_new = x + self.d.reshape(len(self.mesh.nodes), 2) * ratio
        nodal_values = []
        for i in range(len(self.mesh.nodes)):
            nodal_values.append(self.stress_dict[i][stress_idx])

        nodes_x = x_new[:, 0]
        nodes_y = x_new[:, 1]
        elements_tris = []
        for v in self.mesh.elements:
            elements_tris.append([v.node[0], v.node[1], v.node[2]])
            elements_tris.append([v.node[0], v.node[2], v.node[3]])
        triangulation = tri.Triangulation(nodes_x, nodes_y, elements_tris)
        fig, ax = plt.subplots()
        result = ax.tricontourf(triangulation, nodal_values, cmap="jet")
        for v in self.mesh.elements:
            xy_new = [(x_new[idx, 0], x_new[idx, 1]) for idx in v.node]
            patch = patches.Polygon(xy=xy_new, ec="black", fill=False)
            ax.add_patch(patch)
        ax.autoscale()
        ax.set_aspect("equal", "box")
        fig.colorbar(result, ax=ax)
        plt.show()
