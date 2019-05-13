#==============================================================================
# SOLUTION STRATEGY
#
# Method of lines (MOL)
# ---------------------
# Maxwell's equations are first discretized in the spatial variable, thus
# obtaining a system of ordinary differential equations (ODEs) for the time-
# dependent "degrees of freedom" of E and B.
# Given some initial conditions, these "semi-discrete" Maxwell equations can be
# solved approximately by a generic ODE integrator, or by a specialized time
# splitting algorithm.
#
# Operator splitting
# ------------------
# The semi-discrete Maxwell equations are decomposed into the individual Amperè
# and Faraday equations, which are integrated exactly over one time-step.
#==============================================================================

# Import modules
import numpy as np
import matplotlib.pyplot as plt

from ex1_common import *

from splitting_methods import operator_splitting
from splitting_methods import coeffs_order2_strang
from splitting_methods import coeffs_order4_triple_jump
from splitting_methods import coeffs_order6_suzuki

#==============================================================================
# SPATIAL DISCRETIZATION
#==============================================================================
def finite_difference_coeffs(order):
    """ Get coefficients of centered finite difference formula for a given
    order of accuracy. """

    if   order == 2:  return [-1, 1]
    elif order == 4:  return [1/24, -9/8, 9/8, -1/24]
    elif order == 6:  return [-3/640, 25/384, -75/64, 75/64, -25/384, 3/640]
    else:
        raise ValueError('Finite difference of order %s is not available' % order)


def dE_dx(E, dx, coeffs):
    """ Compute ∂E/∂x on the dual grid, given the values of E on the primal grid
    """
    s =  len(coeffs) // 2 - 1
    return sum((c / dx) * np.roll(E, s-k) for k, c in enumerate(coeffs))


def dB_dx(B, dx, coeffs):
    """ Compute ∂B/∂x on the primal grid, given the values of B on the dual grid
    """
    s =  len(coeffs) // 2
    return sum((c / dx) * np.roll(B, s-k) for k, c in enumerate(coeffs))

#==============================================================================
# TIME STEPPING METHOD
#==============================================================================
def operator_splitting_coeffs(order):
    """ Get coefficients for two-way operator splitting of given order of
    accuracy. """

    if   order == 2:  return coeffs_order2_strang
    elif order == 4:  return coeffs_order4_triple_jump
    elif order == 6:  return coeffs_order6_suzuki
    else:
        raise ValueError('Time splitting of order %s is not available' % order)


def step_ampere_1d(dt, E, B, dx, coeffs):
    """ Exactly integrate the semi-discrete Amperè equation over one time-step.
    """
    E -= dt * dB_dx(B, dx, coeffs)
  # B -= 0


def step_faraday_1d(dt, E, B, dx, coeffs):
    """ Exactly integrate the semi-discrete Faraday equation over one time-step.
    """
  # E -= 0
    B -= dt * dE_dx(E, dx, coeffs)

#==============================================================================
# SIMULATION
#==============================================================================
def mol_main(*, L, a, N, Cp, nsteps, plot_interval, spatial_order, splitting_order):

    #--------------------------------------------------------
    # SETUP
    #--------------------------------------------------------

    # Exact solution
    f = exact_solution(0, L, a)

    # Staggered grid (E is on primal grid, B is on dual grid)
    dx, x_prim, x_dual = staggered_grid(0, L, N)

    # Time step size
    dt = Cp * dx

    # Initial time and initial conditions
    t = 0
    E = f(t, x_prim)  # solution array (E_{i})_i
    B = f(t, x_dual)  # solution array (B_{i+1/2})_i

    # Prepare animation
    if plot_interval != 0:
        fig, ax = plt.subplots(1, 1)
        make_plot(ax, t, E, E, x_prim, [0, L] )
        fig.show()

    # Compute energy at t=0
    print()
    print('time = {:8.4f}, energy = {}'.format(t, energy(E, B)))

    # Obtain coefficients for finite difference formulas
    fd_coeffs = finite_difference_coeffs(spatial_order)

    # Obtain coefficients for operator splitting
    os_coeffs = operator_splitting_coeffs(splitting_order)

    #--------------------------------------------------------
    # SOLUTION
    #--------------------------------------------------------

    # Arguments for time stepping
    args = (E, B, dx, fd_coeffs)

    # Time loop
    for i in range(nsteps):

        # Advance both E and B by one time-step
        operator_splitting(os_coeffs, step_faraday_1d, step_ampere_1d, dt, *args)
        t += dt

        # Skip diagnostics if required
        if plot_interval == 0:
            continue

        #... Diagnostics
        if i % plot_interval == 0 or i == nsteps-1:

            # Compute energy
            print('time = {:8.4f}, energy = {}'.format(t, energy(E, B)))

            # Update plot
            exact_E = f(t, x_prim)
            update_plot(ax, t, exact_E, E)
        #...

    #--------------------------------------------------------
    # POST-PROCESSING
    #--------------------------------------------------------

    # Error at final time
    exact_E = f(t, x_prim)
    exact_B = f(t, x_dual)
    error_E = max(exact_E - E)
    error_B = max(exact_B - B)
    print()
    print('Max error on E(t,x) at final time: {}'.format(error_E))
    print('Max error on B(t,x) at final time: {}'.format(error_B))

    # Return whole namespace as dictionary
    return locals()

#==============================================================================
if __name__ == '__main__':

    # Get generic parser from module 'ex1_common.py'
    mol_parser = common_parser()

    # Customize arguments: add order of finite-difference approx. and splitting
    mol_parser.add_argument('-f', '--spatial_order',
        type    = int,
        default = 2,
        choices = [2, 4, 6],
        help    = 'Order of accuracy of finite differences.'
    )
    mol_parser.add_argument('-o', '--splitting_order',
        type    = int,
        default = 2,
        choices = [2, 4, 6],
        help    = 'Order of accuracy of operator splitting.'
    )

    # Read input arguments
    mol_args = mol_parser.parse_args()

    # Run simulation
    mol_namespace = mol_main(**vars(mol_args))
