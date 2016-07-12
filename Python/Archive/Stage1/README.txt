Stage 1 of the project is the set up of a very simple I-Country, S-Periods
lived agents. There are no complicated demographics. Each country produces
the same single good which is moble across country borders. In each country
there is a representative firm. Examine the step by step document in the
notes section for more details.

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
