from __future__ import division
import time
import numpy as np
import AuxiliaryClass as AUX
np.set_printoptions(threshold = 3000, linewidth=2000, suppress=True)

TimeModel=True #Activates timing the model

def Multi_Country(S,I,J,sigma):

    #NOTE:To run the model, simply run the Multi_Country function with your chosen levels
    #of the number of cohorts (S), the number of countries (I), the number of 
    #skill classes (J) and slope parameter (sigma)

    #THIS SETS ALL OF THE USER PARAMETERS

    #Country Rosters
    I_dict = {"usa":0,"eu":1,"japan":2,"china":3,"india":4,"russia":5,"korea":6} #DONT CHANGE
    I_HighSkill = np.array([.3,.3,.3,.25,.25,.3,.3]) #CAN CHANGE
    I_touse = ["eu","russia","usa","japan","korea","china","india"] #CAN CHANGE

    #NOTE: I_HighSkill sets the portion of each country's population that's deemed 
    #"high skill". The remaining portion will be divided equally among the remaining classes.

    #Parameters Zone
    g_A = 0.015 #Technical growth rate
    beta_ann = .95 #Annual discount rate
    delta_ann = .08 #Annual depreciation rate
    alpha = .35 #Capital Share of production
    alphaj = np.array([.25,.4]) #Share of production for each labor class
    chil = .52 #Preference for adult's lesiure 
    chik = 1.0 #Preference for Kids' lesiure
    mu = 2.29 #Unknown Parameter

    #Convergence Tolerances
    demog_ss_tol = 1e-8 #Used in getting ss for population share

    #PROGRAM LEVERS:
    #For terminal output
    PrintAges = False #Displays the current locations of the program inside key TPI functions

    PrintSSEulErrors = True#Prints the euler errors in each attempt of calculating the
                            #steady state
    PrintSS = True  #Prints the result of the Steady State functions
    Print_caTimepaths = False #Prints the consumption, assets, and bequests 
                              #timepath as it gets filled in for each iteration of TPI
    Print_HH_Eulers = True #Prints whether the equations for the household decisions 
                            #are satisfied (Equations 3.22, 3.19, and sum(assets) = 0)
    Print_Fill_Matricies_Time = False #Activiates Printing the total time it takes to 
                                      #fill the upper and lower diagonal matricies
    CheckerMode = False #Activates not printing much of anything, used in conjunction 
                        #with RobustChecker.py
    VerifyDemog = False #Verifies that all of the popluations sum to 1 and that the
                       #Fertility,Mortality and migrant rates copied correctly

    Iterate = True #Shows the current iteration number and the associated Eulers
    ShaveTime = False #Shaves off a little more time for TPI.

    #For plots to display or save
    DemogGraphs = True  #Activates graphing graphs with demographic data 
                        #and population shares
    ShowSSGraphs = True #Activates graphs for steady-state solutions for 
                        #consumption, assets, and bequests
    ShowSSSkill = True
    iterations_to_plot = set([]) #Which iterations of the timepath fsolve you want to plot
    SaveFinalTPIPlot = True #Saves the final (and hopefully converged) time 
                            #path plot as a .png file

    #For using differing ways to solve the model
    UseDiffDemog = True #Turns on different demographics for each country
    UseSSDemog = False #Activates using only steady state demographics for TPI calculation
    UseDiffProductivities = True  #Activates having e vary across cohorts

    #Adjusts the country list if we are using less than 7 Countries
    if len(I_touse) < I:
        print "WARNING: We are changing I from", I, "to", len(I_touse),\
                "to fit the length of I_touse. So the countries we are using now are", I_touse
        I = len(I_touse)
        I_HighSkill = I_HighSkill[:I]
        time.sleep(2)

    elif len(I_touse) > I:
        print "WARNING: We are changing I_touse from", I_touse,\
                "to", I_touse[:I], "so there are", I, "regions"
        I_touse = I_touse[:I]
        I_HighSkill = I_HighSkill[:I]
        time.sleep(2)

    #Does a quick check on labor classes and their share of production, 
    #to make sure the conditions are met
    if len(alphaj)==J:
        if np.sum(alphaj)+alpha==1:
            print "Shares confirmed!"
        else:
            raise ValueError("Production shares MUST sum to 1!")
    else:
        raise ValueError("The number of production shares (alphaj's length)\
                MUST equal the number of classes (J)")


    ##INPUTS INTO THE CLASS###
    Country_Roster = (I_dict, I_touse,I_HighSkill)

    HH_params = (S,I,J,beta_ann,sigma)

    Firm_Params = (alpha, delta_ann, chil, chik, mu, g_A,alphaj)


    Levers = (PrintAges,CheckerMode,Iterate,UseDiffDemog,UseDiffProductivities,\
            Print_Fill_Matricies_Time,ShaveTime)

    #Initialize the class instance
    Model = AUX.OLG(Country_Roster,HH_params,Firm_Params,Levers)

    #Demographics
    Model.Demographics(demog_ss_tol, UseSSDemog,VerifyDemog)
    if DemogGraphs: Model.plotDemographics(T_touse="default", compare_across="T", data_year=0)
    #Model.immigrationplot()

    
    #STEADY STATE OUTER FSOLVE GUESS
    k_ss_guess = np.ones((I))*.25
    kf_ss_guess = np.ones((I-1))*0
    n_ss_guess = np.ones((I,J))*.10
    bq_ss_guess = np.ones((I))*.200 #.195 is the boundary

    #STEADY STATE INNER FSOLVE GUESS
    c_innerfsolve_guess = np.ones((I,J))*.150

    #Steady State
    Model.SteadyState(k_ss_guess,kf_ss_guess,n_ss_guess, bq_ss_guess,c_innerfsolve_guess\
            ,PrintSSEulErrors)

    if PrintSS: Model.PrintSSResults()
    if ShowSSGraphs: Model.plotSSResults(ShowSSSkill)
    
    #Timepath Iteration
    '''
    r_init = Model.r_ss*1.05
    bq_init = Model.bqindiv_ss*.95
    a_init = Model.avec_ss*.7
    
    Model.set_initial_values(r_init, bq_init, a_init)


    Model.Timepath_optimize(Print_HH_Eulers, Print_caTimepaths, iterations_to_plot)
    if SaveFinalTPIPlot: Model.plot_timepaths(SAVE=True)
    '''



#Input parameters for S, I and sigma here then execute this file to
#run the model.

start = time.time()
# S-Number of Cohorts, I-Number of Countries, J-Number of Skill classes
# S, I, J and sigma. S and I are integers. Sigma may not be.
Multi_Country(30,7,2,4)
tottime=time.time()-start

if TimeModel==True:
    minutes=int(tottime/60)
    hours=int(minutes/60)
    days=int(hours/24)
    seconds=tottime-minutes*60
    minutes=minutes-hours*60
    hours=hours-days*24
    print "The code took:",days,"days,", hours, "hours,", minutes, "minutes and", seconds,\
            "seconds to complete"
