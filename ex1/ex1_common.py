import numpy as np
import argparse

#==============================================================================
def staggered_grid(xmin, xmax, N):
    """
    Create a uniform staggered grid on a 1D periodic domain.

    Parameters
    ----------
    xmin : float
        Left boundary of 1D domain.

    xmax : float
        Right boundary of 1D domain.

    N : int
        Number of cells in domain.

    Results
    -------
    dx : float
        Cell size.

    x_prim : 1D numpy.ndarray
        Coordinates of vertices of primal grid.

    x_dual : 1D numpy.ndarray
        Coordinates of vertices of dual grid.
    """

    dx     = (xmax - xmin) / N
    x_prim = dx * np.arange(N)
    x_dual = x_prim + dx / 2

    return dx, x_prim, x_dual

#==============================================================================
def exact_solution(xmin, xmax, a):
    """
    Manufacture a sinusoidal exact solution for the fields E(t,x) and B(t,x)
    on a periodic domain [xmin, xmax].

    Parameters
    ----------
    xmin : float
        Left boundary of 1D domain.

    xmax : float
        Right boundary of 1D domain.

    a : int
        Number of wavelengths in domain.

    Results
    -------
    f : callable
        Exact solution f(t,x) for fields E and B.

    """
    k = a * (2 * np.pi) / (xmax - xmin)
    f = lambda t, x: np.sin(k * (x - t))

    return f

#==============================================================================
def energy(E, B):
    """ Compute total energy in domain given the fields E and B.
    """
    return 0.5 * (np.dot(E, E) + np.dot(B, B))

#==============================================================================
def print_error_analysis(N_list, err_list):

    template_h = '{:^3s} | {:^9s} | {:s}'
    template_0 = '{:3d} | {:9.3e} | {:s}'
    template   = '{:3d} | {:9.3e} | {:4.3f}'

    print(template_h.format('N', 'error', 'order'))
    print(template_0.format(N_list[0], err_list[0], '--'))

    for i in range(1, len(N_list)):
        N = N_list[i]
        N_old = N_list[i-1]
        e = err_list[i]
        e_old = err_list[i-1]

        order = np.log(e_old/e) / np.log(N/N_old)

        print(template.format(N, e, order))

#==============================================================================
def make_plot(ax, t, E_ex, E_num, x_prim, xlim):
    ax.plot(x_prim, E_ex , '--', label='exact')
    ax.plot(x_prim, E_num, '-' , label='numerical')
    ax.legend(loc='upper right')
    ax.grid()
    ax.set_title('Time t = {:10.3e}'.format(t))
    ax.set_xlabel('x')
    ax.set_ylabel('E', rotation='horizontal')
    ax.set_xlim(xlim)

#==============================================================================
def update_plot(ax, t, E_ex, E_num):
    ax.set_title('Time t = {:10.3e}'.format(t))
    ax.lines[0].set_ydata(E_ex )
    ax.lines[1].set_ydata(E_num)
    ax.get_figure().canvas.draw()

#==============================================================================
def common_parser():

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description     = "Setup solution of exercise 1."
    )

    parser.add_argument( 'N',
        type = int,
        help = 'Number of grid cells (elements) in domain'
    )

    parser.add_argument( '-l',
        type    = float,
        default = 1.0,
        dest    = 'L',
        metavar = 'L',
        help    = 'Length of periodic domain [0, L]'
    )

    parser.add_argument( '-a',
        type    = int,
        default = 3,
        help    = 'Number of wavelengths in sinusoidal profile (initial conditions)'
    )

    parser.add_argument( '-c',
        type    = float,
        default = 0.5,
        dest    = 'Cp',
        metavar = 'Cp',
        help    = 'Courant parameter'
    )

    parser.add_argument( '-t',
        type    = int,
        default = 1,
        dest    = 'nsteps',
        metavar = 'NSTEPS',
        help    = 'Number of time-steps to be taken'
    )

    parser.add_argument( '-p',
        type    = int,
        default = 4,
        metavar = 'I',
        dest    = 'plot_interval',
        help    = 'No. of time steps between successive plots of solution, if I=0 no plots are made'
    )

    return parser
