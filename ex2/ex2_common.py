import numpy as np

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
def print_error_analysis(N_list, err_list):
    """
    Create (and print) simple table for error convergence under mesh refinement.

    Parameters
    ----------
    N_list : list of int
        Increasing values of N, number of cells in domain.

    err_list : list of float
        Errors obtained with the given values of N.

    """

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
