from __future__ import division
import time
import numpy as np  
import AuxiliaryClass as AUX
np.set_printoptions(threshold = np.nan, linewidth=2000, suppress=True)

TimeModel=True #Activates timing the model

def Multi_Country(S,I,sigma):

    #NOTE:To run the model, simply run the Multi_Country function with your chosen levels
    #of the number of cohorts (S), the number of countries (I) and slope parameter (sigma)

    #THIS SETS ALL OF THE USER PARAMETERS

    #Country Rosters
    I_dict = {"usa":0,"eu":1,"japan":2,"china":3,"india":4,"russia":5,"korea":6} #DONT CHANGE
    I_touse = ["eu","russia","usa","japan","korea","china","india"] #CAN CHANGE

    #Parameters Zone
    g_A = 0.015 #Technical growth rate
    beta_ann=.95 #Annual discount rate
    delta_ann=.08 #Annual depreciation rate
    alpha = .25 #Capital Share of production
    chil = .52
    chik = 1.0
    mu = 2.29

    Dem_Degree = int(4) #Must be an integer, 4 is the max, otherwise there will be collinearity

    #Convergence Tolerances
    demog_ss_tol = 1e-8 #Used in getting ss for population share

    #PROGRAM LEVERS:
    #For terminal output
    PrintAges = False #Displays the current locations of the program inside key TPI functions
    PrintSSEulErrors = True #Prints the euler errors in each attempt of 
                             #calculating the steady state
    PrintSS = True #Prints the result of the Steady State functions
    Print_caTimepaths = False #Prints the consumption, assets, and bequests 
                              #timepath as it gets filled in for each iteration of TPI
    Print_HH_Eulers = False #Prints whether the equations for the household decisions 
                            #are satisfied (Equations 3.22, 3.19, and sum(assets) = 0)
    Print_Fill_Matricies_Time = False #Activiates Printing the total time it takes to 
                                      #fill the upper and lower diagonal matricies
    Print_Reg_Coefficients = False #Prints the coefficients from the regression done to
                                   #fit a polynomial to our regression data.
    CheckerMode = False #Activates not printing much of anything, used 
                        #in conjunction with RobustChecker.py
    Iterate = True #Shows the current iteration number and the associated Eulers
    ShaveTime = False #Shaves off a little more time for TPI.
    SaveToSTATA = True #Saves the Fertility and Mortality rate matricies to a 
                        #form STATA can read. With this, we can verify Python's results.

    #For plots to display or save
    DemogGraphs = False #Activates graphing graphs with demographic data and population shares
    ShowSSGraphs = True #Activates graphs for steady-state solutions for 
                        #consumption, assets, and bequests
    ShowCompGraphs = True #Shows the graph which compares calculated 
                           #Mortality/Fertility rates against Their non-calculated counterparts
    iterations_to_plot = set([]) #Which iterations of the timepath fsolve you want to plot
    SaveFinalTPIPlot = True #Saves the final (and hopefully converged) time 
                            #path plot as a .png file

    #For using differing ways to solve the model
    UseDiffDemog = True #Turns on different demographics for each country
    UseSSDemog = False #Activates using only steady state demographics for TPI calculation
    UseCalcDemog = True #Uses the calculated mortality and fertilty rates based on regression
    UseDiffProductivities = False #Activates having e vary across cohorts

    #Adjusts the country list if we are using less than 7 Countries
    if CheckerMode==False:
        if len(I_touse) < I:
            print "WARNING: We are changing I from", I, "to", len(I_touse),\
                    "to fit the length of I_touse. So the countries we are using now are",\
                    I_touse
            I = len(I_touse)
            time.sleep(2)

        elif len(I_touse) > I:
            print "WARNING: We are changing I_touse from", I_touse, "to", I_touse[:I],\
                    "so there are", I, "regions"
            I_touse = I_touse[:I]
            time.sleep(2)

    ##INPUTS INTO THE CLASS###
    Country_Roster = (I_dict, I_touse)

    HH_params = (S,I,beta_ann,sigma,Dem_Degree)

    Firm_Params = (alpha, delta_ann, chil, chik, mu, g_A)

    Levers = (PrintAges,CheckerMode,Iterate,UseDiffDemog,UseDiffProductivities,\
            Print_Fill_Matricies_Time,ShaveTime,UseCalcDemog,ShowCompGraphs,Print_Reg_Coefficients,SaveToSTATA)

   #Initialize the class instance
    Model = AUX.OLG(Country_Roster,HH_params,Firm_Params,Levers)
    #Demographics
    Model.Demographics(demog_ss_tol, UseSSDemog)
    if DemogGraphs: Model.plotDemographics(T_touse="default", compare_across="T", data_year=0)
    #Model.immigrationplot()

    #STEADY STATE OUTER FSOLVE INITIAL GUESSES
    k_ss_guess = np.ones(I)*.55
    kf_ss_guess = np.ones(I-1)*.05
    n_ss_guess = np.ones(I)*1.25
    bq_ss_guess = np.ones(I)*.85

    #STEADY STATE INNER FSVOLE INITIAL GUESSES
    c_innerfsolve_guess = np.ones(I)*.05

    #Steady State
    Model.SteadyState(k_ss_guess,kf_ss_guess,n_ss_guess, bq_ss_guess,\
            c_innerfsolve_guess, PrintSSEulErrors)
    if PrintSS: Model.PrintSSResults()
    if ShowSSGraphs: Model.plotSSResults()
    '''

    #Timepath Iteration
    r_init = Model.r_ss*1.05
    bq_init = Model.bqindiv_ss*.95
    a_init = Model.avec_ss*.7

    
    Model.set_initial_values(r_init, bq_init, a_init)

    Model.Timepath_optimize(Print_HH_Eulers, Print_caTimepaths, iterations_to_plot)
    if SaveFinalTPIPlot: Model.plot_timepaths(SAVE=False)

    '''
#Input parameters for S, I and sigma here then execute this file to
#run the model.

start = time.time()
# S-Number of Cohorts, I-Number of Countries
# S, I, and sigma. S and I are integers. Sigma may not be.
Multi_Country(20,7,4)
tottime=time.time()-start

if TimeModel==True:
    minutes=int(tottime/60)
    hours=int(minutes/60)
    seconds=tottime-minutes*60
    minutes=minutes-hours*60
    print "The code took:", hours, "hours,", minutes, "minutes and", seconds,\
            "seconds to complete"


