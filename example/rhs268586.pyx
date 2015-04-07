# This file is generated automatically by QuTiP.
# (C) Paul D. Nation & J. R. Johansson

from numpy import *
cimport libc.math as cmath

import numpy as np
cimport numpy as np
cimport cython
from qutip.cy.spmatfuncs import spmv_csr, spmvpy

ctypedef np.complex128_t CTYPE_t
ctypedef np.float64_t DTYPE_t



@cython.boundscheck(False)
@cython.wraparound(False)

def cy_td_ode_rhs(double t, np.ndarray[CTYPE_t, ndim=1] vec, np.ndarray[CTYPE_t, ndim=1] data0, np.ndarray[int, ndim=1] idx0, np.ndarray[int, ndim=1] ptr0):
    
    cdef Py_ssize_t row
    cdef int num_rows = len(vec)
    cdef np.ndarray[CTYPE_t, ndim=1] out = np.zeros((num_rows),dtype=np.complex)
     
    spmvpy(data0, idx0, ptr0, vec, 1.0, out)
    return out
