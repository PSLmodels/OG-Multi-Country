from numpy cimport ndarray as ar
from libc cimport math
from libc.math cimport exp
cimport cython
from cython.parallel cimport prange
"""
Please note that numpy isn't called because in Cython exp is faster than np.exp
when using Cython.

"""

@cython.boundscheck(False)
@cython.wraparound(False) 
#These two lines disable safety checks that take time to run. MAKE SURE THE ARRAYS YOU PASS IN ARE LEGAL OR ELSE YOU'LL CRASH THE PROGRAM
def cy_fillca(ar[double, ndim=3] C_Mat, ar[double,ndim=3] Ck_Mat, ar[double,ndim=3]A_Mat,ar[double, ndim=1] r_path,ar[double,ndim=3] m_rates,ar[double, ndim=3] bq_vecpath, ar[double, ndim=3] w_e, ar[double,ndim=3] gamma_mat,ar[double,ndim=1] lbar, ar[double,ndim=3] Kids, double beta, double chi, double delta, double g_A, double rho, double sigma):
    """
        Description:
            - Cython version of the consumption and assets fill loop in AuxiliaryClass.py. It simply fills
              the given consumption and assets matricies in cython. This module is marginally faster than
              the current pure python version.

        Inputs:
            - A_Mat         = Array: [I,S+1,T+S], Savings decisions for each cohort
            - bq_vecpath    = Array: [I,S,T+S], Transition path for distribution of bequests for each country
            - C_Mat         = Array: [I,S,T+S], Consumption decisions for each cohort
            - Ck_Mat        = Array: [I,S,T+S], Kids consumption decisions for each cohort
            - gamma_mat     = Array: [I,S,T+S], Transition path of shorthand calcuation variable gamma
            - Kids          = Array: [I,S,T], Matrix of the number of kids for each country
            - lbar          = Array: [T+S], Transition path of the time endowment
            - m_rates       = Array: [I,S,T+S], Mortality rates of each country for each age cohort and year
            - r_path        = Array: [T+S], Transition path for the interest rate
            - w_e           = Array: [I,S,T+S], Matrix product of w and e
            - beta          = Scalar: Calculated overall future discount rate
            - chi           = Scalar: Leisure preference parameter
            - delta         = Scalar: Calculated overall depreciation rate
            - g_A           = Scalar: Growth rate of technology
            - rho           = Scalar: The intratemporal elasticity of substitution between consumption and leisure
            - sigma         = Scalar: Rate of Time Preference

        Variables Called From Object:
            - None

        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Outputs:
            - None
    """
    cdef int T,s #In cython, you must declare integers
    T=C_Mat.shape[2]-C_Mat.shape[1] #Extracting T from the dimensions of the inputted matrix
    
    for s in prange(C_Mat.shape[1]-1, nogil=True): #Using Python's built in parallel module
        
        with gil:
            Ck_Mat[:,s+1,s+1:T+s+1] = ((beta * (1-m_rates[:,s,s:T+s]) * (1 + r_path[s+1:T+s+1] - delta))**(1/sigma)\
                    * Ck_Mat[:,s,s:T+s])*exp(-g_A)
            
            C_Mat[:,s+1,s+1:T+s+1] = Ck_Mat[:,s+1,s+1:T+s+1]/gamma_mat[:,s+1,s+1:T+s+1]
            
            A_Mat[:,s+1,s+1:T+s+1] = (  (w_e[:,s,s:T+s]*lbar[s:T+s] + (1 + r_path[s:T+s] - delta)*A_Mat[:,s,s:T+s] + bq_vecpath[:,s,s:T+s])\
                    -C_Mat[:,s,s:T+s]*(1+Kids[:,s,s:T+s]*gamma_mat[:,s,s:T+s]+w_e[:,s,s:T+s]*(chi/w_e[:,s,s:T+s])**(rho) )  )*exp(-g_A)
    

