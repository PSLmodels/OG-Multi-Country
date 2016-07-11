This folder contains the versions of code of the S country OLG model.
These folders are:

7CountryAlphaV1— Contains all the code of the Fortran-Python translation. This was the 
                 original way that this model was coded in fortran. Refer to the README in
                 that folder for what the individual files do. This was abandoned because of
                 the difficulties of translating from Fortran to Python, namely, the
                 translation of indices. This version doesn't fully work properly.

7CountryAlphaV2— Contains the first attempt at coding the model from the ground up in Python
                 this was abandoned because of difficulties around iterating around 
                 consumption rather than iterating around savings, which is much simpler.
                 It does successfully calculate the steady state, however, it's missing most
                 of the complex features in the fortran code.

7CountryAlphaV3— Contains the current efforts on the code, which is currently adding skill classes. In this version,
		  we transition to using a python-class to more easily organize all of the components of the model
		  Adding skill classes hasn’t progressed beyond the steady state because of the necessary
		  change in the solution algorithm.


7CountryElliptical— A temporary offshoot of 7CountryAlphaV3, instead of the utility function we start with
		     in stages 1-4, we try an Elliptical utility function based on Phillips and Evans (2015).
		     In addition to a new utility function, this folder includes new demographics that can
		     be easily adjusted based on the number of generations. Also, this hasn’t gone beyond the
		     steady state calculations.

Archive- Contains the archive of COMPLETED stages of code. We follow the stages as outlined
         in OG-Multi-Country/Notes. As of 7/11/16, this archive contains stages 1-4. Stage
	  1 is a simple S-period I-Country OLG Model, Stage 2 adds demographics and growth
	  3 adds a Leisure decision and Stage 4 adds children to the model.
