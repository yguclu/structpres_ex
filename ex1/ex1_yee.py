# Import modules
import numpy as np
import matplotlib.pyplot as plt

from ex1_common import *

#==============================================================================
# Standard Yee's scheme:
#   . 2nd-order finite differences in space
#   . 2nd-order leapfrog in time
#
# Note that E[i] = E_{i} and B[i] = B_{i+1/2} for i = 0, 1, ..., N - 1.
#==============================================================================

def yee_step_E(dt, dx, E, B):
    """ Evolve electric field according to equation 2.18 in lecture notes.
    """
    E += (dt / dx) * (np.roll(B, +1) - B)


def yee_step_B(dt, dx, E, B):
    """ Evolve magnetic field according to equation 2.19 in lecture notes.
    """
    B += (dt / dx) * (E - np.roll(E, -1))

#==============================================================================
def yee_main(*, L, a, N, Cp, nsteps, plot_interval):

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

    #--------------------------------------------------------
    # SOLUTION
    #--------------------------------------------------------

    # Push B backward in time by half time-step
    yee_step_B(-dt/2, dx, E, B)

    # Time loop
    for i in range(nsteps):

        # Advance both E and B by one time-step
        yee_step_B(dt, dx, E, B)  # Bring B from t - ∆t/2 to t + ∆t/2
        yee_step_E(dt, dx, E, B)  # Bring E from t to t + ∆t
        t += dt

        # Skip diagnostics if required
        if plot_interval == 0:
            continue

        #... Diagnostics
        if i % plot_interval == 0 or i == nsteps-1:

            # Bring B to same time instant as E
            yee_step_B(dt/2, dx, E, B)

            # Compute energy
            print('time = {:8.4f}, energy = {}'.format(t, energy(E, B)))

            # Update plot
            exact_E = f(t, x_prim)
            update_plot(ax, t, exact_E, E)

            # Push B back by half time-step
            yee_step_B(-dt/2, dx, E, B)
        #...

    # Bring B to same time instant as E
    yee_step_B(dt/2, dx, E, B)

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

    # Read input arguments
    yee_parser = common_parser()
#    description     = "Solve the 1D Maxwell equations using Yee's scheme."
    yee_args = yee_parser.parse_args()

    # Run simulation
    yee_namespace = yee_main(**vars(yee_args))
