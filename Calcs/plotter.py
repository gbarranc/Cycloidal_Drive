#Plotting mixin. Uses Geometry + Vectors methods via the host object.
import matplotlib.pyplot as plt

class plotter:
    def plotted(self):
        (up, low) = self.point_slope()
        x_limits = (up[:, 0], low[:, 0])
        y_extended = (up[:, 1], low[:, 1])

        Cx, Cy = self.geometry()
        shell = self.rollers()
        center = self.center_of_rollers()
        Bx, By = self.bore()
        X_OR, Y_OR = self.Output_Reactions()
        OR, OP, _ = self.output_rollers()

        plt.figure(figsize=(8, 8))
        plt.xlim(-70, 70)
        plt.ylim(-70, 70)

        for xs, ys in ((Bx, By), (Cx, Cy), (OR[:, 0], OR[:, 1]), (OP[:, 0], OP[:, 1])):
            plt.scatter(xs, ys, color='green', marker='.', s=0.5)

        plt.scatter(center[:, 0], center[:, 1], color='orange', marker='.', s=0.5)
        plt.scatter(shell[:, 0], shell[:, 1], color='blue', marker='.', s=0.5)
        plt.scatter(self.contacts[:, 0], self.contacts[:, 1], color='red', marker='.', s=100)
        plt.plot(x_limits, y_extended, color='red')
        plt.plot(X_OR, Y_OR, color='red')

        plt.title("2D Cycloid Profile")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(True)
        plt.show()