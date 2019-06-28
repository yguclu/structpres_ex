import numpy as np
from scipy.interpolate import lagrange

#==============================================================================
# LAGRANGE INTERPOLATION
#==============================================================================

class Basis_V0:

    def __init__(self, degree):

        if degree < 1:
            raise ValueError('Degree must be strictly positive.')
        elif degree % 2 != 1:
            raise ValueError('Degree must be odd.')

        p = degree // 2
        x = [*range(-p, p+2)]
        y = lambda j: [(1 if j == k else 0) for k in range(-p, p+2)]

        self._p      = p
        self._degree = degree
        self._basis  = tuple(lagrange(x, y(j)) for j in range(-p, p+2))

    # ...
    def __len__(self):
        return self._degree + 1

    # ...
    def __iter__(self):
        return iter(self._basis)

    # ...
    def __getitem__(self, j):
        p = self._p
        if not -p <= j <= p+1:
            raise ValueError('Basis function index j can only take values between -p and p+1.')
        return self._basis[p + j]

    # ...
    @property
    def degree(self):
        return self._degree

#++++++++++++++++++++++++++++++++++++++++++++++
class Interpolator_from_C0_to_V0:

    def __init__(self, xgrid, basis):
        self._xgrid = xgrid
        self._N     = len(xgrid)
        self._dx    = xgrid[1] - xgrid[0]
        self._basis = basis

    # ...
    def __call__(self, values):
        return Interpolant_V0(self, values)

#++++++++++++++++++++++++++++++++++++++++++++++
class Interpolant_V0:

    def __init__(self, interpolator, values):
        assert(interpolator._N == len(values))
        self._interpolator = interpolator
        self._values       = np.asarray(values)

    #...
    def __call__(self, x):
        x0   = self._interpolator._xgrid[0]
        dx   = self._interpolator._dx
        N    = self._interpolator._N
        l    = self._interpolator._basis
        psi0 = self._values

        X  = (x - x0) / dx
        i  = np.array(np.floor(X), dtype=int)
        xi = X - i
    
        p = l.degree // 2

        return sum(psi0[(i + j) % N] * l[j](xi) for j in range(-p, p+2))

#==============================================================================
# LAGRANGE HISTOPOLATION
#==============================================================================

class Basis_V1:

    def __init__(self, degree):

        if degree < 1:
            raise ValueError('Degree must be strictly positive.')
        elif degree % 2 != 1:
            raise ValueError('Degree must be odd.')

        p = degree // 2
        c = lambda i, j: 0.5 * float(np.sign(j - i - 0.5))
        l = Basis_V0(degree)

        self._p      = p
        self._degree = degree
        self._basis  = tuple(
            sum(c(i, j) * l[j] for j in range(-p, p+2)).deriv()
            for i in range(-p, p+1)
        )

    # ...
    def __len__(self):
        return self._degree

    # ...
    def __iter__(self):
        return iter(self._basis)

    #...
    def __getitem__(self, j):
        p = self._p
        if not -p <= j <= p:
            raise ValueError('Basis function index j can only take values between -p and p.')
        return self._basis[p + j]

    #...
    @property
    def degree(self):
        return self._degree

#++++++++++++++++++++++++++++++++++++++++++++++
class Interpolator_from_C1_to_V1:

    def __init__(self, xgrid, basis):
        self._xgrid = xgrid
        self._N     = len(xgrid)
        self._dx    = xgrid[1] - xgrid[0]
        self._basis = basis

    # ...
    def __call__(self, values):
        return Interpolant_V1(self, values)

#++++++++++++++++++++++++++++++++++++++++++++++
class Interpolant_V1:

    def __init__(self, interpolator, values):
        assert(interpolator._N == len(values))
        self._interpolator = interpolator
        self._values       = np.asarray(values)

    def __call__(self, x):
        x0   = self._interpolator._xgrid[0]
        dx   = self._interpolator._dx
        N    = self._interpolator._N
        e    = self._interpolator._basis
        psi1 = self._values

        X  = (x - x0) / dx
        i  = np.array(np.floor(X), dtype=int)
        xi = X - i
    
        p = e.degree // 2

        return sum(psi1[(i + j) % N] * e[j](xi) for j in range(-p, p+1))
