import numpy as np
import Stage2Functions as Stepfuncs
import time as time

np.set_printoptions(threshold = 3000, linewidth=2000, suppress=True)

def Multi_Country(S,I,sigma):

    #Parameters Zone
    I_all = ["usa","eu","japan","china","india","russia","korea"]
    #I = 3 #Number of countries
    #S = 80 #Upper bound of age for agents
    T = int(round(2.5*S)) #Number of time periods to convergence, based on Rick Evans' function.
    I_touse = ["usa","eu","japan","russia","korea","china","india"]

    T_1 = S #This is like TransYear in the FORTRAN I think
    if S > 50:
        T_1 = 50

    g_A = 0.015 #Technical growth rate

    beta_ann=.95 #Starting future consumption discount rate
    delta_ann=.08 #Starting depreciation rate
    beta = beta_ann**(70/S) #Future consumption discount rate
    delta = 1-(1-delta_ann)**(70/S) #Depreciation Rate
    alpha = .3 #Capital Share of production

    tpi_tol = 1e-8 #Convergence Tolerance
    demog_ss_tol = 1e-8 #Used in getting ss for population share
    xi = .8 #Parameter used to take the convex conjugate of paths
    MaxIters = 50000000 #Maximum number of iterations on TPI.

    #Program Levers
    CalcTPI = True #Activates the calculation of Time Path Iteration

    PrintAges = False #Prints the different key ages in the demographics
    PrintLoc = False #Displays the current locations of the program inside key TPI functions
    PrintEulErrors = False #Prints the euler errors in each attempt of calculating the steady state
    PrintSS = False #Prints the result of the Steady State functions
    Print_cabqTimepaths = False #Prints the consumption, assets, and bequests timepath as it gets filled in for each iteration of TPI

    DemogGraphs = False #Activates graphing graphs with demographic data and population shares
    TPIGraphs = True #Activates graphing the graphs.

    UseStaggeredAges = True #Activates using staggered ages
    UseDiffDemog = True #Turns on different demographics for each country
    UseSSDemog = False #Activates using only steady state demographics for TPI calculation
    UseDiffProductivities = False #Activates having e vary across cohorts
    UseTape = True #Activates setting any value of kd<0 to 0.001 in TPI calculation
    SAVE = True #Saves the graphs
    SHOW = False #Shows the graphs

    LeaveHouseAge, FirstFertilityAge, LastFertilityAge, MaxImmigrantAge, FirstDyingAge, agestopull = Stepfuncs.getkeyages(S, PrintAges, UseStaggeredAges)

    if len(I_touse) < I:
        print "WARNING: We are changing I from", I, "to", len(I_touse), "to fit the length of I_touse"
        I = len(I_touse)	
        time.sleep(2)

    if UseDiffDemog:
        A = np.ones(I)+np.cumsum(np.ones(I)*.08)-.08 #Techonological Change, used for when countries are different
    else:
        A = np.ones(I) #Techonological Change, used for idential countries

    if UseDiffProductivities:
        e = np.ones((I, S, T+S))
        e[:,FirstDyingAge:,:] = 0.01
        e[:,:LeaveHouseAge,:] = 0.01
    else:
        e = np.ones((I, S, T+S)) #Labor productivities

    #MAIN CODE

    #Gets demographic data
    demog_params = (I, S, T, T_1, LeaveHouseAge, FirstFertilityAge, LastFertilityAge, FirstDyingAge, MaxImmigrantAge, agestopull, g_A, demog_ss_tol)
    demog_levers = PrintLoc, UseStaggeredAges, UseDiffDemog, DemogGraphs
    MortalityRates, Nhat_matrix, Nhat_ss = Stepfuncs.getDemographics(demog_params, demog_levers, I_all, I_touse)

    #Initalizes initial guesses
    assets_guess = np.ones((I, S-1))*.1
    kf_guess = np.zeros((I))

    #Gets the steady state variables
    params_ss = (I, S, beta, sigma, delta, alpha, e[:,:,-1], A, FirstFertilityAge, FirstDyingAge, Nhat_ss, MortalityRates[:,:,-1], g_A, PrintEulErrors)
    assets_ss, kf_ss, kd_ss, n_ss, y_ss, r_ss, w_ss, c_vec_ss = Stepfuncs.getSteadyState(params_ss, assets_guess, kf_guess)

    if PrintSS==True: #Prints the results of the steady state, line 23 activates this
        print "assets steady state", assets_ss
        print "kf steady state", kf_ss
        print "kd steady state", kd_ss
        print "n steady state",n_ss
        print "y steady state", y_ss
        print "r steady state",r_ss
        print "w steady state", w_ss
        print "c_vec_ss steady state",c_vec_ss

    if UseSSDemog == True:
        print "NOTE: USING SS DEMOGRAPHICS FOR TIMEPATH\n"
        Nhat_matrix = np.einsum("is,t->ist", Nhat_matrix[:,:,-1],np.ones(T+S))
        MortalityRates = np.einsum("is,t->ist", MortalityRates[:,:,-1],np.ones(T+S))
        lbar[:] = 1
        time.sleep(2)

    if CalcTPI==True: #Time Path Iteration, activated by line 24
        print "Beginning TPI..."
        #Gets initial guesses for TPI
        initialguess_params = (I, S, T, delta, alpha, e[:,:,0], A, FirstFertilityAge, FirstDyingAge, Nhat_matrix[:,:,0], MortalityRates[:,:,0], g_A)
        assets_init, wpath_initguess, rpath_initguess = \
            Stepfuncs.get_initialguesses(initialguess_params, assets_ss, kf_ss, w_ss, r_ss, PrintLoc)

        #Gets timepaths for w, r, C, K, and Y
        tp_params = (I, S, T, T_1, beta, sigma, delta, alpha, e, A, FirstFertilityAge, FirstDyingAge, Nhat_matrix, MortalityRates, g_A, tpi_tol, xi, MaxIters)
        wpath, rpath, Cpath, Kpath, Ypath = Stepfuncs.get_Timepath(tp_params, wpath_initguess, rpath_initguess, assets_init, kd_ss, kf_ss, PrintLoc, Print_cabqTimepaths, UseTape)

        print rpath
        if TPIGraphs==True:
            Stepfuncs.plotTimepaths(I, S, T, sigma, wpath, rpath, Cpath, Kpath, Ypath, I_touse, SAVE, SHOW)


Multi_Country(80,7,2)