import time
import numpy as np
import AuxiliaryClass as AUX

def Multi_Country(S,I,sigma):

    #NOTE:To run the model, simply run the Multi_Country function with your chosen levels
    #of the number of cohorts (S), the number of countries (I) and slope parameter (sigma)

    #THIS SETS ALL OF THE USER PARAMETERS

    #Country Rosters
    I_dict = {"usa":0,"eu":1,"japan":2,"china":3,"india":4,"russia":5,"korea":6}
    I_touse = ["eu","russia","usa","japan","korea","china","india"] 

    #Parameters Zone
    g_A = 0.015 #Technical growth rate
    beta_ann=.95 #Annual discount rate
    delta_ann=.08 #Annual depreciation rate    
    alpha = .3 #Capital Share of production
    chi = 1.5 #Preference for lesiure
    rho = 1.3 #Other New Parameter

    #Convergence Tolerances
    tpi_tol = 1e-8 #Convergence Tolerance
    demog_ss_tol = 1e-8 #Used in getting ss for population share
    xi = .9999 #Parameter used to take the convex conjugate of paths
    MaxIters = 10000 #Maximum number of iterations on TPI.

    #Program Levers
    CalcTPI = False #Activates the calculation of Time Path Iteration

    PrintAges = False #Prints the different key ages in the demographics
    PrintLoc = False #Displays the current locations of the program inside key TPI functions
    PrintEulErrors = False #Prints the euler errors in each attempt of calculating the steady state
    PrintSS = True #Prints the result of the Steady State functions
    Print_cabqTimepaths = False #Prints the consumption, assets, and bequests timepath as it gets filled in for each iteration of TPI

    CheckerMode = False #Reduces the number of prints when checking for robustness, use in conjunction with RobustChecker.py

    DemogGraphs = False #Activates graphing graphs with demographic data and population shares
    TPIGraphs = False #Activates graphing the graphs.

    UseStaggeredAges = True #Activates using staggered ages
    UseDiffDemog = True #Turns on different demographics for each country
    UseSSDemog = False #Activates using only steady state demographics for TPI calculation
    UseDiffProductivities = True #Activates having e vary across cohorts
    UseTape = True #Activates setting any value of kd<0 to 0.001 in TPI calculation
    SAVE = False #Saves the graphs
    SHOW = True #Shows the graphs
    ADJUSTKOREAIMMIGRATION = True #Adjusts demograhpics to correct for oddities in Korea's data.

    #Adjusts the country list if we are using less than 7 Countries
    if len(I_touse) < I:
        print "WARNING: We are changing I from", I, "to", len(I_touse), "to fit the length of I_touse. So the countries we are using now are", I_touse
        I = len(I_touse)
        time.sleep(2)
    elif len(I_touse) > I:
        print "WARNING: We are changing I_touse from", I_touse, "to", I_touse[:I], "so there are", I, "regions"
        I_touse = I_touse[:I]
        time.sleep(2)



    ##INPUTS INTO THE CLASS###

    Country_Roster= (I_dict, I_touse)

    HH_params = (S,I,beta_ann,sigma)

    Firm_Params = (alpha, delta_ann, chi, rho, g_A)

    Tolerances = (tpi_tol, demog_ss_tol)

    Levers = (CalcTPI, PrintAges,PrintLoc,PrintEulErrors,PrintSS,Print_cabqTimepaths,CheckerMode,DemogGraphs,TPIGraphs,\
            UseStaggeredAges,UseDiffDemog,UseSSDemog,UseDiffProductivities,UseTape,SAVE,SHOW,ADJUSTKOREAIMMIGRATION)

    TPI_Params = (xi,MaxIters)



    ##WHERE THE MAGIC HAPPENS ##

    Model = AUX.OLG(Country_Roster,HH_params,Firm_Params,Levers, Tolerances, TPI_Params)

    #Demographics
    Model.Import_Data()
    Model.Demographics()

    #Steady State

    #STEADY STATE INITIAL GUESSES

    r_ss_guess = .2
    bq_ss_guess = np.ones(I)*.2
    Model.SteadyState(r_ss_guess, bq_ss_guess)


    #Timepath Iteration


    return None


#Input parameters for S, I and sigma here then execute this file to
#run the model.
Multi_Country(20,2,2)

