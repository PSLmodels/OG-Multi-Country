Version 2 was when we started to code up the model from the ground up. This
was we started to use a python class, so one will have to keep close track
of all the variables being passed through each of the functions. Version 2
was abandoned to change what variable we iterate around. This folder is
not currently actively edited.

Data_Files: This folder contains the data files from Kotlikoff's research.
            They include the fertility and mortality rates for each country 
            over time. Net migration and population are included as well.

Main.py- The central file where the user can change their inputs. By running
         this script, the whole model will run.

RobustChecker.py- MPI script that tests the model. It goes through different
                  combinations of numbers of countries, number of cohorts
                  and sigma levels. It outputs .CSV files that will indicate
                  if that particular combination converged.

StepbyStepv1.py- This file is where the calculations occur. It's broken
                 down into steady state and TPI sections. Examine the
                 file for more details. Note thate demographics are NOT
                 incorporated, but they are coded up.
