from numpy cimport ndarray as ar
from libc.math cimport exp
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cy_fillca(ar[double, ndim=3] C_Mat, ar[double,ndim=3]A_Mat,ar[double, ndim=1] r_path,ar[double,ndim=3] m_rates,ar[double, ndim=3] bq_vecpath, ar[double, ndim=3] w_e, ar[double,ndim=3] psi_mat, double beta, double chi, double delta, double g_A, double rho, double sigma):
    cdef int I,S,T,Q,s
    I=C_Mat.shape[0]
    S=C_Mat.shape[1]
    Q=C_Mat.shape[2]
    T=Q-S
    
    for s in xrange(S-1):
        
        C_Mat[:,s+1,s+1:T+s+1] = ((beta * (1-m_rates[:,s,s:T+s]) * (1 + r_path[s+1:T+s+1] - delta)\
                * psi_mat[:,s+1,s+1:T+s+1])/psi_mat[:,s,s:T+s])**(1/sigma) * C_Mat[:,s,s:T+s]*exp(-g_A)
        
        A_Mat[:,s+1,s+1:T+s+1] = (  (w_e[:,s,s:T+s] + (1 + r_path[s:T+s] - delta)*A_Mat[:,s,s:T+s] + bq_vecpath[:,s,s:T+s])\
                -C_Mat[:,s,s:T+s]*(1+w_e[:,s,s:T+s]*(chi/w_e[:,s,s:T+s])**rho)  )*exp(-g_A)
    

