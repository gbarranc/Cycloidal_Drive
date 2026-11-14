
# Takes fundamental cycloidal-gearbox geometry (radii, eccentricity, pin
# counts), finds the contact points with the eccentricity parallel to the
# x axis, locates where the contact normals converge, and uses that for a
# force balance and Hertzian contact stress check.

import numpy as np
import tkinter as tk
from tkinter import messagebox

from build_geometry import Geometry
from solve_vectors import Vectors
from plotter import plotter


class Cycloid(Geometry, Vectors, plotter):
    # material properties
    E1 = 69.5 # GPa
    E2 = 69.5 # GPa
    

    v1 = v2 = 0.33
    YIELD = 275          # MPa
    INPUT_TORQUE = 10 * 1000   # N-mm  (10 N-m)

    def __init__(self, R_1, Rr_1, e, Z1, Z2, B, OR, OP):
        self.R_1, self.Rr_1 = R_1, Rr_1
        self.e, self.B, self.OR, self.OP = e, B, OR, OP
        self.Z1, self.Z2 = Z1, Z2
        self.deg_shift = 2 * np.pi / self.Z2
        self.x = np.arange(0, 2 * np.pi, 0.0001)

        self.contacts = self.Closest_Pair()                       # run once
        self.forces = self.force_balance(self.INPUT_TORQUE)
        self.max_stress = self.Contact_Stresses()

    def report(self):
        passed = "FAILED" if self.max_stress / self.YIELD < 0 else "PASSED"
        root = tk.Tk()
        root.withdraw()
        messagebox.showinfo(
            "Results",
            f"Max force at inner pins: {round(self.forces[2], 2)}N located at pin: {round(self.forces[3], 2)}"
            f"\nwith a max VM stress of: {round(self.max_stress, 2)} MPa"
            f"\n{passed}"
            f"\nPin 0 starts on rightmost x intercept and spans CCW")
        print('\n\n\nMax force at outer pins:', self.forces[0],
              'N which is at pin:', self.forces[1])


if __name__ == "__main__":
    main = Cycloid(R_1=6.35 / 2, Rr_1=104.14 / 2, e=2.54 / 2,
                   Z1=25, Z2=26, B=20, OR=9.5885, OP=7.9375)
    main.report()
    main.plotted()