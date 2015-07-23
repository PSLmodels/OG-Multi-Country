This folder contains the versions of code of the S country OLG model.
These folders are:

7CountryAlphaV1-- Contains all the code of the Fortran-Python translation. This was the 
                original way that this model was coded in fortran. Refer to the README in
                that folder for what the individual files do. This was abandoned because of
                the difficulties of translating from Fortran to Python, namely, the
                translation of indices. This version doesn't fully work properly.

7CountryAlphaV2-- Contains the first attempt at coding the model from the ground up in Python
                this was also abadoned because of difficulties around iterating around 
                consumption rather than iterating around savings, which is much simpler.
                It does successfully calculate the steady state, however, it's missing most
                of the complex features in the fortran code.

7CountryCurrent-- Contains the most current version of the model. Refer to the README in the
                folder for more details
