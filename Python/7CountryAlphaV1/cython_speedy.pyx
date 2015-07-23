from numpy cimport ndarray as ar
cdef extern from "speedy.h":
    int test(int a, int b)
cpdef testy(int a, int b):
    return test(a,b)

 
