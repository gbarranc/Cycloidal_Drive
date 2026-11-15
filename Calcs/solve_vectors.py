# Contact solving, force balance and Hertzian stresses.

# Mixin class -- relies on Geometry methods (f, g, center_of_rollers,
# output_rollers) being available on the host object via inheritance.

import numpy as np
from scipy.optimize import fsolve


class Vectors:
    # contact points
    def Closest_Pair(self):
        #Contact [x, y, slope] for each roller -> (Z2, 3) array
        centers = self.center_of_rollers()
        contacts = np.zeros((self.Z2, 3))
        k_vals = np.linspace(0, 2 * np.pi, 500)
        l_vals = np.linspace(0, 2 * np.pi, 500)

        for i in range(self.Z2):
            self.Ex, self.Ey = centers[i]   # read by self.g()

            # coarse sweep for initial guess
            best, best_params = np.inf, (0, 0)
            for k in k_vals:
                fx, fy = self.f(k)
                for l in l_vals:
                    gx, gy = self.g(l)
                    d = (fx - gx) ** 2 + (fy - gy) ** 2
                    if d < best:
                        best, best_params = d, (k, l)

            def equations(params):
                k, l = params
                fx, fy = self.f(k)
                gx, gy = self.g(l)
                return [fx - gx, fy - gy]

            solution, info, ier, msg = fsolve(equations, best_params, full_output=True)
            k_sol = solution[0]

            if abs(info['fvec'][0]) <= 2:
                contacts[i, 0:2] = np.double(self.f(k_sol))
                contacts[i, 2] = self.normals(k_sol)
        return contacts

    def normals(self, k_sol, eps=1e-6):
        #Slope of the line from the contact point to the force convergence point (e*Z2, 0)
        x1, y1 = self.f(k_sol - eps)
        x2, y2 = self.e * self.Z2, 0
        return (y2 - y1) / (x2 - x1)

    def point_slope(self):
        #End points of each contact normal line, for plotting
        xyUP = np.zeros((self.Z2, 2))
        xyLOW = np.zeros((self.Z2, 2))
        for i in range(self.Z2):
            x, y, m = self.contacts[i]
            if x == 0:
                continue
            b = y - x * m
            xyUP[i] = (x, m * x + b)
            xyLOW[i] = (self.e * self.Z2, 0)
        return np.array([xyUP, xyLOW])

    # output side
    def Output_Reactions(self):
        rad = self.OUTPUT_PIN_RADIUS
        num = self.OUTPUT_PIN_COUNT
        centers = np.zeros((num, 2))
        reactions = np.zeros((num, 2))
        deg_shift = 2 * np.pi / num
        for i in range(num):
            centers[i, 0] = rad * np.cos(i * deg_shift)
            reactions[i, 0] = centers[i, 0] - self.OP
            centers[i, 1] = rad * np.sin(i * deg_shift)
            reactions[i, 1] = centers[i, 1]
        return (centers[:, 0], reactions[:, 0]), (centers[:, 1], reactions[:, 1])

    # forces / stresses
    def force_balance(self, T_input):
        #Returns (max outer force, outer index, max inner force, inner index)
        T_output = T_input * self.Z1
        m_outer = self.contacts[:, 2]
        x_outer = self.contacts[:, 0]
        y_outer = self.contacts[:, 1]
        _, _, M_inner = self.output_rollers()

        r_outer = np.abs(m_outer * x_outer - y_outer) / np.sqrt(m_outer ** 2 + 1)
        r_inner = M_inner

        F_outer = T_input * r_outer / np.sum(r_outer ** 2) / 2    # two disks
        F_inner = T_output * r_inner / np.sum(r_inner ** 2) / 4

        return (np.max(F_outer), np.argmax(F_outer),
                np.max(F_inner), np.argmax(F_inner))

    def Contact_Stresses(self):
        F_c = np.max(self.forces)
        D1 = self.OR / 1000
        D2 = self.OP / 1000
        L = 6.35 / 1000

        E_e = 1e9 * ((1 - self.v1 ** 2) / self.E1 + (1 - self.v2 ** 2) / self.E2) ** -1
        D_e = 2 * (D1 * D2 / (D1 - D2))

        P_max = np.sqrt((2 * F_c * E_e) / (np.pi * L * D_e)) / 1e6
        Shear_max = 0.3 * P_max
        Tensile_max = (1 - 2 * self.v1) * P_max / 3
        return np.sqrt(Tensile_max ** 2 + 3 * Shear_max ** 2)