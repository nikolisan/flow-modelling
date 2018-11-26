from matplotlib import pyplot as plt
from matplotlib.pyplot import cm


def contour_plot(x, y, u):
    plt.contourf(x, y, u, alpha=0.7, cmap=cm.coolwarm)
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    import numpy as np
    X = np.arange(-5, 5, 0.25)
    Y = np.arange(-5, 5, 0.25)
    X, Y = np.meshgrid(X, Y)
    R = np.sqrt(X ** 2 + Y ** 2)
    Z = np.sin(R)
    contour_plot(X, Y, Z)
