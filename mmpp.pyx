"""
Speed up Markov modulated Poisson process calculations.

The implementation is in Cython for speed
and uses python numpy arrays for speed and convenience.
For compilation instructions see
http://docs.cython.org/src/reference/compilation.html
For example:
$ cython -a mmpp.pyx
$ gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
      -I/usr/include/python2.7 -o mmpp.so mmpp.c
"""

from cython.view cimport array as cvarray
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport exp, log, sqrt

np.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def get_mmpp_block(double a, double w, double r, double t):
    """
    Compute differential equations evaluated at the given time.
    @param a: rate from off to on
    @param w: rate from on to off
    @param r: poisson event rate
    @param t: elapsed time
    @return: P
    """
    cdef double x, denom, p, q
    cdef np.ndarray[np.float64_t, ndim=2] P = np.empty(
            (2, 2), dtype=np.float64)

    # first row of the output ndarray
    x = sqrt((a + r + w)*(a + r + w) - 4*a*r)
    denom = 2 * x * exp(t * (x + a + r + w) / 2)
    p = (exp(t*x)*(x + r + w - a) + (x - r - w + a)) / denom
    q = (2 * a * (exp(t * x) - 1)) / denom
    P[0, 0] = p
    P[0, 1] = q

    # second row of the output ndarray
    x = sqrt((a + r + w)*(a + r + w) - 4*a*r)
    denom = 2 * x * exp(t * (x + a + r + w) / 2)
    p = (2 * w * (exp(t * x) - 1)) / denom
    q = (exp(t*x)*(x - r - w + a) + (x + r + w - a)) / denom
    P[1, 0] = p
    P[1, 1] = q

    return P

