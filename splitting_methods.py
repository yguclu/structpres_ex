# coding: utf-8
# Copyright 2019 Yaman Güçlü, IPP Garching

"""
Module containing splitting methods for the time integration of autonomous
systems of ODEs of the form

dU/dt = A(U) + B(U),

where we assume that the two split systems

1) dU/dt = A(U)
2) dU/dt = B(U)

can be solved exactly. Under this assumptions a splitting method is completely
defined by its set of splitting coefficients, which should be given to function
'operator_splitting'. For convenience we provide here the splitting coefficients
of some famous methods of order 2, 4 and 6.

References
----------
a) E. Sonnendrücker, Chapter 3 of "Geometric methods for the physics of
   magnetised plasmas", Lecture notes, TUM winter semester 2018/2019.
   (http://www-m16.ma.tum.de/foswiki/pub/M16/Allgemeines/GemMeth18/GeomMeth.pdf)

b) R.I. McLachlan and G.R.W. Quispel. Splitting methods. Acta Numerica, Vol. 11,
   pp. 341–434, 2002. (https://doi.org/10.1017/S0962492902000053)

"""

__all__ = ['operator_splitting',
           'coeffs_order2_strang',
           'coeffs_order4_triple_jump',
           'coeffs_order6_suzuki']

#==============================================================================
def operator_splitting(coeffs, step_A, step_B, dt, *args):
    """
    Advance solution by one time-step dt by alternating the action of two
    separate operators A and B. First A is applied, then B, then A again,
    and so on.

    Please note that the solution (which is included in *args) is updated
    in-place by the two operators A and B.

    Parameters
    ----------
    coeffs : list/tuple of float
        Coefficients defining the time-step fraction for each stage.
        The sum of these must by an integer number equal to the number of split
        operators, hence 2 in this case.

    step_A : callable
        Function that updates the current solution by applying the operator A
        for a given amount of time. Signature is step_A(c*dt, *args) where 'c'
        is one of the values of 'coeffs', 'dt' is the time-step size and 'args'
        is a list of positional arguments (which includes the various components
        of the solution as well as other constant parameters or operators).

    step_B : callable
        Function that updates the current solution by applying the operator B
        for a given amount of time. Signature is step_B(c*dt, *args) where 'c'
        is one of the values of 'coeffs', 'dt' is the time-step size and 'args'
        is a list of positional arguments (which includes the various components
        of the solution as well as other constant parameters or operators).

    dt : float
        Time-step size: solution evolves from time (t) to time (t + dt).

    *args : list
        Positional arguments for functions 'step_A' and 'step_B'.

    """
    odd = True
    for c in coeffs:
        step = step_A if odd else step_B
        step(c*dt, *args)
        odd = not odd

#==============================================================================

# 2nd-order Strang splitting
coeffs_order2_strang = (0.5, 1, 0.5)

# 4th-order "triple-jump" symmetric composition, applied to Strang splitting
coeffs_order4_triple_jump = (
        0.675603595979829,
        1.351207191959658,
       -0.175603595979828,
       -1.702414383919315,
       -0.175603595979828,
        1.351207191959658,
        0.675603595979829)

# 6th-order symmetric composition by Suzuki, applied to Strang splitting
coeffs_order6_suzuki = (
        0.0414649985182624 ,
        0.123229775946271  ,
        0.198128671918067  ,
        0.290553797799558  ,
       -0.0400061921041533 ,
       -0.127049212625417  ,
        0.0752539843015807 ,
       -0.246331761062075  ,
       -0.0115113874206879 ,
        0.357208872795928  ,
        0.23666992478693111,
        0.20477705429147008,
        0.23666992478693111,
        0.357208872795928  ,
       -0.0115113874206879 ,
       -0.246331761062075  ,
        0.0752539843015807 ,
       -0.127049212625417  ,
       -0.0400061921041533 ,
        0.290553797799558  ,
        0.198128671918067  ,
        0.123229775946271  ,
        0.0414649985182624 )

