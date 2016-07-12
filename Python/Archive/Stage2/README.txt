Stage 2 incorporates demographics into the calculations. This is split
into two parts, one where the demographics have been tied into the calculation
and where it hasn't been done yet. It's split because we wanted to make sure
that the demographics where coded up correctly before putting incorporating
them into the calculations.

---------------------------------------------------------------------------
Unincorporated Demographics Folder
---------------------------------------------------------------------------

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

---------------------------------------------------------------------------
Incorporated Demographics Folder
---------------------------------------------------------------------------

Data_Files: This folder contains the data files from Kotlikoff's research.
            They include the fertility and mortality rates for each country 
            over time. Net migration and population are included as well.

Stage2Main.py- The central file where the user can change their inputs. 
               By running this script, the whole model will run.

RobustChecker.py- MPI script that tests the model. It goes through different
                  combinations of numbers of countries, number of cohorts
                  and sigma levels. It outputs .CSV files that will indicate
                  if that particular combination converged.

Stage2Functions.py- This file is where the calculations occur. It's broken
                 down into steady state and TPI sections. Examine the
                 file for more details. Demographics are fully incorporated
                 and are functioning.
