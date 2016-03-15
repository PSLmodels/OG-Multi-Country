from numpy cimport ndarray as ar
from libc cimport math
from libc.math cimport exp
cimport cython
from cython.parallel cimport prange

@cython.boundscheck(False)
@cython.wraparound(False)
def cy_fillca(ar[double, ndim=3] C_Mat, ar[double,ndim=3] Ck_Mat, ar[double,ndim=3]A_Mat,ar[double, ndim=1] r_path,ar[double,ndim=3] m_rates,ar[double, ndim=3] bq_vecpath, ar[double, ndim=3] w_e, ar[double,ndim=3] gamma_mat,ar[double,ndim=1] lbar, ar[double,ndim=3] Kids, double beta, double chi, double delta, double g_A, double rho, double sigma):
    cdef int T,s
    T=C_Mat.shape[2]-C_Mat.shape[1]
    
    for s in prange(C_Mat.shape[1]-1, nogil=True):
        
        with gil:
            Ck_Mat[:,s+1,s+1:T+s+1] = ((beta * (1-m_rates[:,s,s:T+s]) * (1 + r_path[s+1:T+s+1] - delta))**(1/sigma)\
                    * Ck_Mat[:,s,s:T+s])*exp(-g_A)
            
            C_Mat[:,s+1,s+1:T+s+1] = Ck_Mat[:,s+1,s+1:T+s+1]/gamma_mat[:,s+1,s+1:T+s+1]
            
            #Gets assets for every agents' next year using Equation 3.19
            A_Mat[:,s+1,s+1:T+s+1] = (  (w_e[:,s,s:T+s]*lbar[s:T+s] + (1 + r_path[s:T+s] - delta)*A_Mat[:,s,s:T+s] + bq_vecpath[:,s,s:T+s])\
                    -C_Mat[:,s,s:T+s]*(1+Kids[:,s,s:T+s]*gamma_mat[:,s,s:T+s]+w_e[:,s,s:T+s]*(chi/w_e[:,s,s:T+s])**(rho) )  )*exp(-g_A)
    

