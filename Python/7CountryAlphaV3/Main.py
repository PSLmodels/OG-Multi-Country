from __future__ import division
import time
import numpy as np
import AuxiliaryClass as AUX
np.set_printoptions(threshold = 3000, linewidth=2000, suppress=True)

TimeModel=True #Activates timing the model

def Multi_Country(S,I,J,sigma):

    #NOTE:To run the model, simply run the Multi_Country function with your chosen levels
    #of the number of cohorts (S), the number of countries (I), the number of skill classes (J) and slope parameter (sigma)

    #THIS SETS ALL OF THE USER PARAMETERS

    #Country Rosters
    I_dict = {"usa":0,"eu":1,"japan":2,"china":3,"india":4,"russia":5,"korea":6} #DONT CHANGE
    #TO JEFF: I this the percentages below are for the lower skill class just FYI. 
    #         So all we'll need to do is change the variable name to I_LowSkill or something
    I_HighSkill = np.array([.7,.7,.7,.75,.75,.7,.7]) #CAN CHANGE
    I_touse = ["eu","russia","usa","japan","korea","china","india"] #CAN CHANGE

    #NOTE: I_HighSkill sets the portion of each country's population that's deemed 
    #      "high skill". The remaining portion will be divided equally among the remaining classes.

    #Parameters Zone
    g_A = 0.015 #Technical growth rate
    beta_ann=.95 #Annual discount rate
    delta_ann=.08 #Annual depreciation rate
    alpha = .3 #Capital Share of production
    alphaj = np.array([.4,.3]) #Share of production for each labor class
    chi = 1.5 #Preference for lesiure
    rho = .4 #Intratemporal elasticity of substitution


    #Convergence Tolerances
    demog_ss_tol = 1e-8 #Used in getting ss for population share

    #PROGRAM LEVERS:
    #For terminal output
    PrintAges = False #Displays the current locations of the program inside key TPI functions
    PrintSSEulErrors = False #Prints the euler errors in each attempt of calculating the steady state
    PrintSS = False #Prints the result of the Steady State functions
    Print_caTimepaths = False #Prints the consumption, assets, and bequests timepath as it gets filled in for each iteration of TPI
    Print_HH_Eulers = False #Prints whether the equations for the household decisions are satisfied (Equations 3.22, 3.19, and sum(assets) = 0)
    Print_Fill_Matricies_Time = False #Activiates Printing the total time it takes to fill the upper and lower diagonal matricies
    CheckerMode = False #Activates not printing much of anything, used in conjunction with RobustChecker.py
    Iterate = True #Shows the current iteration number and the associated Eulers
    ShaveTime = False #Shaves off a little more time for TPI.

    #For plots to display or save
    DemogGraphs = False #Activates graphing graphs with demographic data and population shares
    ShowSSGraphs = False #Activates graphs for steady-state solutions for consumption, assets, and bequests
    iterations_to_plot = set([]) #Which iterations of the timepath fsolve you want to plot
    SaveFinalTPIPlot = True #Saves the final (and hopefully converged) time path plot as a .png file

    #For using differing ways to solve the model
    UseDiffDemog = True #Turns on different demographics for each country
    UseSSDemog = False #Activates using only steady state demographics for TPI calculation
    UseDiffProductivities = True #Activates having e vary across cohorts
    UseSamePopRates = True #Activates using the same Demographics across labor classes. Leave as true.

    #Adjusts the country list if we are using less than 7 Countries
    if CheckerMode==False:
        if len(I_touse) < I:
            print "WARNING: We are changing I from", I, "to", len(I_touse), "to fit the length of I_touse. So the countries we are using now are", I_touse
            I = len(I_touse)
            I_HighSkill = I_HighSkill[:I]
            time.sleep(2)

        elif len(I_touse) > I:
            print "WARNING: We are changing I_touse from", I_touse, "to", I_touse[:I], "so there are", I, "regions"
            I_touse = I_touse[:I]
            I_HighSkill = I_HighSkill[:I]
            time.sleep(2)

    #Does a quick check on labor classes and their share of production, to make sure the conditions are met
    if len(alphaj)==J:
        if np.sum(alphaj)+alpha==1:
            print "Shares confirmed!"
        else:
            raise ValueError("Production shares MUST sum to 1!")
    else:
        raise ValueError("The number of production shares (alphaj's length) MUST equal the number of classes (J)")


    ##INPUTS INTO THE CLASS###
    Country_Roster = (I_dict, I_touse,I_HighSkill)

    HH_params = (S,I,J,beta_ann,sigma)

    Firm_Params = (alpha, delta_ann, chi, rho, g_A,alphaj)

    Levers = (PrintAges,CheckerMode,Iterate,UseDiffDemog,UseDiffProductivities,Print_Fill_Matricies_Time,ShaveTime,UseSamePopRates)

    #Initialize the class instance
    Model = AUX.OLG(Country_Roster,HH_params,Firm_Params,Levers)

    #Demographics
    Model.Demographics(demog_ss_tol, UseSSDemog)
    if DemogGraphs: Model.plotDemographics(T_touse="default", compare_across="T", data_year=0)
    #Model.immigrationplot()

    #STEADY STATE INITIAL GUESSES
    r_ss_guess = .25
    bq_ss_guess = np.ones((I,J))*.2

    #Steady State
    Model.SteadyState(r_ss_guess, bq_ss_guess, PrintSSEulErrors)
    if PrintSS: Model.PrintSSResults()
    if ShowSSGraphs: Model.plotSSResults()

    #Timepath Iteration
    
    r_init = Model.r_ss*1.05
    bq_init = Model.bqindiv_ss*.95
    a_init = Model.avec_ss*.7
    
    Model.set_initial_values(r_init, bq_init, a_init)


    Model.Timepath_optimize(Print_HH_Eulers, Print_caTimepaths, iterations_to_plot)
    if SaveFinalTPIPlot: Model.plot_timepaths(SAVE=True)



#Input parameters for S, I and sigma here then execute this file to
#run the model.

start = time.time()
# S-Number of Cohorts, I-Number of Countries, J-Number of Skill classes
# S, I, J and sigma. S and I are integers. Sigma may not be.
Multi_Country(80,7,2,4)
tottime=time.time()-start

if TimeModel==True:
    minutes=int(tottime/60)
    hours=int(minutes/60)
    seconds=tottime-minutes*60
    minutes=minutes-hours*60
    print "The code took:", hours, "hours,", minutes, "minutes and", seconds, "seconds to complete"


