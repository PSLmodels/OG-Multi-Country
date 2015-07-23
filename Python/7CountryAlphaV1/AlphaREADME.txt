This folder contains our initial attempts to code up this 7 Country Model. It was done, initially, by translating the Fortran code that is in this repository.

This code has the following structure:

WorldModule.py--This contains all of the arrays and variables that were supposed to be used.

AlphaFunctions.py--This contains all of the functions that would be used in calculating the steady state and transition path.

CentralCode.py--This is the central code that controls what the overall model does.

cython_speedy.pyx--This was an abandoned module that was written in C to help speed up the process. Because of the nature of the fortran code with tons of for loops, it was our intention to use C to speed up the process if needed.

demographicswithclasses -- Contains demographics simulation code using each region as a seperate region class

demographicswithclasses2 -- Similiar to demographicswithclasses except it calculates immigration/new births/deaths in a slightly different order than demographicswithclasses. This version is probably better, but it depends on nuts and bolts of how we do the populations projections

We didn't get a chance to fully debug this before changing direction to a ground up approach. This will mostly give out NaNs.:w

