import numpy
import sys


def sin_wave(t, prop=0.0):
    global A
    return A * numpy.sin(2 * numpy.pi * t / T - prop)


def variables(Nx, Ny, Lx, Ly, c, Tin, A0):
    global A, T
    A = A0
    T = Tin
    x = numpy.linspace(0, Lx, Nx + 3)
    y = numpy.linspace(0, Ly, Ny + 2)
    dx = Lx / Nx
    dy = Ly / Ny
    dt = round(dx * dy / c / (dx + dy))
    print('Calculated Dx: {}'.format(dx))
    print('Calculated Dy: {}'.format(dy))
    print('Calculated Dt: {}'.format(dt))
    if (T / dt) < 10:  # CFL Condition
        print('The calculated Dt is not at least an order\n'
              'of magnitude less than the tidal wave period.\n'
              'Please change Nx and/or Ny')
        sys.exit()
    return x, y, dx, dy, dt


def calculate_water_level(Lx, Nx, Ny, dx, dy, dt, c, days, P=None):
    # P : tuple of a point's coord to get the time series
    print('Calculating water level profile for {} days...'.format(days))
    # Useful variables
    Cxs = (c * dt / dx) ** 2  # Courant number for X - squared
    Cys = (c * dt / dy) ** 2  # Courant number for Y - squared

    # Instantiating arrays
    # numpy.zeros(rows, columns)
    u = numpy.zeros((Ny + 2, Nx + 3))    # u(n+1)(i,j)   -> solution array
    u_1 = numpy.zeros((Ny + 2, Nx + 3))  # u(n)(i,j)     -> solution at t-1
    u_2 = numpy.zeros((Ny + 2, Nx + 3))  # u(n-1)(i,j)   -> solution at t-2

    # For time series plot
    u_ = []
    t_ = []

    t = 0
    n = 0
    # Step 1
    t = t + dt
    n = n + 1
    for i in range(1, Ny + 1):
        for j in range(1, Nx + 2):
            u_xx = u_1[i + 1, j] - 2 * u_1[i, j] + u_1[i - 1, j]
            u_yy = u_1[i, j + 1] - 2 * u_1[i, j] + u_1[i, j - 1]
            u[i, j] = u_1[i, j] + 0.5 * Cxs * u_xx + 0.5 * Cys * u_yy

    # Loop Start
    while True:
            t = t + dt
            n = n + 1
            for i in range(1, Ny + 1):
                for j in range(1, Nx + 2):
                    u_xx = u_1[i + 1, j] - 2 * u_1[i, j] + u_1[i - 1, j]
                    u_yy = u_1[i, j + 1] - 2 * u_1[i, j] + u_1[i, j - 1]
                    u[i, j] = 2 * u_1[i, j] - u_2[i, j] + Cxs * u_xx + Cys * u_yy

            # Open sea boundary (p. 73)
            for j in range(Nx + 2):
                if t > 2*(Lx / c):  # if the reflected wave have travelled back
                    # to the open sea boundary
                    u1 = u_1[0, j] - sin_wave(t - dt)  # (2.38)
                    u2 = u_1[1, j] - sin_wave(t - dt, dx / Lx)  # (2.39)
                    u1 = u1 + dt / dx * c * (u2 - u1)  # (2.42)
                else:
                    u1 = 0
                u[0, j] = u1 + sin_wave(t)

            # Closed boundary with ghost cells
            # Fully reflective boundaries (Neumann Conditions)
            # * top (i = Ny)
            u[-1, :] = u[-2, :]
            # * left (j = Nx)
            u[:, -1] = u[:, -2]
            # * right (j = 0)
            u[:, 0] = u[:, 1]

            if P is not None:
                xp = P[0]
                yp = P[1]
                u_.append(u[xp, yp])
                t_.append(t)

            u_2 = numpy.copy(u_1)
            u_1 = numpy.copy(u)

            if t > 3600 * 24 * days:
                break
    return u, u_, t_