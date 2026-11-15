# Pure geometry: profile curves, roller/pin outlines and centers.

# Mixin class -- expects the host object (Cycloid) to provide:
#     R_1, Rr_1, e, Z1, Z2, B, OR, OP, x, deg_shift
# No imports from the other project modules (prevents circular imports).

import numpy as np


class Geometry:
    # Output pin circle (hard-coded)
    OUTPUT_PIN_RADIUS = 36.98
    OUTPUT_PIN_COUNT = 8

    # ---------- cycloid disc profile ----------
    def _contact_angle(self, k):
        #Contact angle Psi as a function of the sweep parameter k
        return np.arctan(np.sin(self.Z1 * k) /
                         (np.cos(self.Z1 * k) - self.Rr_1 / (self.e * (self.Z1 + 1))))

    def f(self, k):
        #Cycloid profile point (x, y) at parameter k (scalar or array)
        psi = self._contact_angle(k)
        x = (self.Rr_1 * np.cos(k) - self.R_1 * np.cos(k + psi)
             - self.e * np.cos(k * (self.Z1 + 1)))
        y = (-self.Rr_1 * np.sin(k) + self.R_1 * np.sin(k + psi)
             + self.e * np.sin(k * (self.Z1 + 1)))
        return x, y

    def geometry(self):
        #Full cycloid profile over self.x
        return self.f(self.x)

    def g(self, l):
        # Shrunken roller circle (R_1 - scale) around the current roller
        # center (self.Ex, self.Ey); used by the contact solver
        scale = 2
        return ((self.R_1 - scale) * np.cos(l) + (self.Ex - self.e),
                (self.R_1 - scale) * np.sin(l) + self.Ey)

    # housing rollers
    def center_of_rollers(self):
        centers = np.zeros((self.Z2, 2))
        for i in range(self.Z2):
            centers[i, 0] = self.Rr_1 * np.cos(i * self.deg_shift) - self.e
            centers[i, 1] = self.Rr_1 * np.sin(i * self.deg_shift)
        return centers

    def rollers(self):
        #360 outline points per roller
        size = (self.Z2 + 1) * 360
        shell = np.zeros((size, 2))
        centers = self.center_of_rollers()
        theta = np.linspace(0, 2 * np.pi, 360, endpoint=False)
        for i in range(self.Z2):
            start = i * 360
            shell[start:start + 360, 0] = self.R_1 * np.cos(theta) + centers[i, 0]
            shell[start:start + 360, 1] = self.R_1 * np.sin(theta) + centers[i, 1]
        return shell

    # bore / output pins
    def bore(self):
        return self.B * np.cos(self.x) + self.e, self.B * np.sin(self.x)

    def output_rollers(self):
        #Returns (Output_Rollers, Output_pins, M_inner moment arms)
        n, rad = self.OUTPUT_PIN_COUNT, self.OUTPUT_PIN_RADIUS
        size = (self.Z2 + 1) * 360 + 1
        rollers = np.zeros((size, 2))
        pins = np.zeros((size, 2))
        theta = np.linspace(0, 2 * np.pi, 360, endpoint=False)
        step = 2 * np.pi / n
        M_inner = np.zeros((n, 1))

        for k in range(n):
            s = k * 360
            cx, cy = np.cos(k * step) * rad, np.sin(k * step) * rad
            rollers[s:s + 360, 0] = self.OR * np.cos(theta) + cx + self.e
            rollers[s:s + 360, 1] = self.OR * np.sin(theta) + cy
            pins[s:s + 360, 0] = self.OP * np.cos(theta) + cx
            pins[s:s + 360, 1] = self.OP * np.sin(theta) + cy
            M_inner[k] = cy
        return rollers, pins, M_inner