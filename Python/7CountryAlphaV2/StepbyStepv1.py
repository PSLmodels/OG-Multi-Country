from __future__ import division
import sys
import csv
import numpy as np
import scipy as sp
import scipy.optimize as opt
import time as time
from matplotlib import pyplot as plt

#DEMOGRAPHICS FUNCTIONS
def getkeyages(S, PrintAges, UseStaggeredAges):
    """
    Description:
        -Gets the key ages for calculating demographics based on S

    Inputs:
        S                = Scalar in [10,80], Number of cohorts
        PrintAges        = Bool, Prints key ages if true
        UseStaggeredAges = Bool, Gets ages to use from a linspace rather than a sequence of integers

    Functions called:
        -None

    Objects in Function:
        agestopull        = (S) List, Contains the ages from the data files to use for the demographic data in the program
        LeaveHouseAge     = Scalar, Age where children become adults and are no longer dependent on parents for consumption
        FirstFertilityAge = Scalar in (0,S), First age when agents bear children
        LastFertilityAge  = Scalar in (0,S), Last age when agents bear children
        FirstDyingAge     = Scalar in (0,S), First age when agents die
        MaxImmigrantAge   = Scalar in (0,S), Age of the oldest immigrants
       
    Returns: LeaveHouseAge, FirstFertilityAge, LastFertilityAge, MaxImmigrantAge, FirstDyingAge, agestopull
    """
    if UseStaggeredAges:
        agestopull = np.arange(80)[np.where(np.in1d(np.arange(80), np.round(np.linspace(0, 80, num=S, endpoint=False))))]
        LeaveHouseAge = np.abs(agestopull - 21).argmin()#The age when agents become independent households and don't rely on parents consumption
        FirstFertilityAge = np.abs(agestopull - 23).argmin()#The age when agents have their first children
        LastFertilityAge = np.abs(agestopull - 45).argmin()#The age when agents have their last children
        MaxImmigrantAge = np.abs(agestopull - 65).argmin()#All immigrants are between ages 0 and MaxImmigrantAge
        FirstDyingAge = np.abs(agestopull - 68).argmin()#The first age agents can begin to die
    else:
        agestopull = range(S)
        LeaveHouseAge = int(np.round(S/80.*21))#The age when agents become independent households and don't rely on parents consumption
        FirstFertilityAge = int(np.round(S/80.*23))#The age when agents have their first children
        LastFertilityAge = int(np.round(S/80.*45))#The age when agents have their last children
        MaxImmigrantAge = int(np.round(S/80.*65))#All immigrants are between ages 0 and MaxImmigrantAge
        FirstDyingAge = int(np.round(S/80.*68))#The first age agents can begin to die

    #Insures that parents don't die with kids in the home
    if FirstDyingAge - (LeaveHouseAge + LastFertilityAge) <= 0:
        FirstDyingAge+=1

    #Makes sure we aren't pulling anything out of bounds from the data files
    if agestopull[FirstDyingAge] < 68:
        agestopull[FirstDyingAge] = 45
    if agestopull[FirstFertilityAge] < 23:
        agestopull[FirstFertilityAge] = 23
    if agestopull[LastFertilityAge] > 45:
        agestopull[LastFertilityAge] = 45

    if PrintAges:
        print "All ages to pull:", agestopull
        print "Adjusted to S =", S, ":    LeaveHouseAge =", LeaveHouseAge, "FirstFertilityAge =", FirstFertilityAge, \
        "LastFertilityAge =", LastFertilityAge, "FirstDyingAge =", FirstDyingAge, "MaxImmigrantAge =", MaxImmigrantAge
        if UseStaggeredAges:
            print "To pull from the data:  LeaveHouseAge =", agestopull[LeaveHouseAge], "FirstFertilityAge =", agestopull[FirstFertilityAge], \
            "LastFertilityAge =", agestopull[LastFertilityAge], "FirstDyingAge =", agestopull[FirstDyingAge], "MaxImmigrantAge =", agestopull[MaxImmigrantAge],"\n"

    return LeaveHouseAge, FirstFertilityAge, LastFertilityAge, MaxImmigrantAge, FirstDyingAge, agestopull

def getDemographics(params, UseStaggeredAges, DiffDemog, Graphs, I_all, I_touse):
    """
    Description:
        -Imports and stores data from csv files for initial populations, fertility rates, mortality rates, and net migrants. 
        -Calculates population distribuitons through year T and its steady-state
        -Calculates number of kids living at home for each country, generation, and year and its steady-state
        -Calculates labor endowment for each year

    Inputs:
        params            = Tuple, Contains parameters I, S, T, T_1, LeaveHouseAge, FirstFertilityAge, LastFertilityAge, 
                            FirstDyingAge, MaxImmigrantAge, g_A, and tol
        I                 = Scalar in [1,7], Number of countries
        S                 = Scalar in [10,80], Number of cohorts
        T                 = Scalar >0, Number of years away from time t=0 until we reach the steady-state
        T_1               = Scalar >0, Number of years away from time t=0 until the demographics are stationarized
        LeaveHouseAge     = Scalar, Age where children become adults and are no longer dependent on parents for consumption
        FirstFertilityAge = Scalar in (0,S), First age when agents bear children
        LastFertilityAge  = Scalar in (0,S), Last age when agents bear children
        FirstDyingAge     = Scalar in (0,S), First age when agents die
        MaxImmigrantAge   = Scalar in (0,S), Age of the oldest immigrants
        g_A               = Scalar >0, Technical growth rate
        tol               = Scalar >0, Tolerance level below which the change population share is considered to be in the steady-state
        UseStaggeredAges  = Boolean, Uses a linspace of ages from csv files if true, uses range(S) for ages if false
        DiffDemog         = Boolean, Causes data to load from different csv files for each country if==True
        Graphs            = Boolean, Plots ss population, and fertility, mortality, and immigration rates at time t=0
        countrynames      = Tuple of length (I), contains the names of each of the I countries, used for pulling data from csv files

    Functions called:
        -None

    Objects in Function:
        f_range            = Scalar >0, Number of years the agents are fertile
        N_matrix           = [I, S+1, T+S+1] Matrix, Population for each country, generation, and year
        Nhat               = [I, S+1, T+S+1] Matrix, Population share for each country, generation, and year
        all_FertilityRates = [I, S+1, f_range+T+S+1] Matrix, Fertility rates for each country, generation, and year. 
                             The first (f_range) years in the third dimension are before year t=0
                             and are used to calculate KIDs of generations alive at time t=0.
        FertilityRates     = [I, S+1, T+S+1] Matrix, Fertility rates for each country, generation, and year. 
                                             Same as all_FertilityRates, but without the first f_range years
        MortalityRates     = [I, S+1, T+S+1] Matrix, Mortality rates for each country, generation, and year.
        ImmigrationRates   = [I, S+1, T+S+1] Matrix, Immigration rates for each country, generation, and year.
        Migrants           = [I, S+1, T+S+1] Matrix, Net migration for each country, generation, and year.
        g_N                = [T+S+1,] Vector, Population growth rate for each year
        N_temp             = [I, S+1] Matrix, QUESTION WHAT IS THIS???
        index              = Scalar in [0,I], Used in for loop to gather csv demographic data. Always 0 if DiffDemog==False
        f_bar              = [S+1,] Vector, Stationarized fertility rate for each country and year beyond T_1
        rho_bar            = [S+1,] Vector, Stationarized mortality rate for each country and year beyond T_1
        Nhat_ss            = [I, S+1] Matrix, Steady-state for population share
        pop_old            = [I, S+1] Matrix, Compared to pop_new for converging to the steady-state
        pop_new            = [I, S+1] Matrix, Compared to pop_old for converging to the steady-state
        iter               = Scalar >=0, Counts the number of years beyond T+S+1 until Nhat converges to the steady-state
        lbar               = [T+S+1,] Vector, Labor endowment for each year

    Returns: MortalityRates, Nhat, KIDs, Nhatss_new, KIDs_ss, lbar
    """
    #Unpacks parameters
    I, S, T, T_1, LeaveHouseAge, FirstFertilityAge, LastFertilityAge, FirstDyingAge, MaxImmigrantAge, agestopull, g_A, tol = params

    #Parameter that represents total number of fertile years. Used for indexing purposes and to make the code easier to read
    f_range = LastFertilityAge+1-FirstFertilityAge

    #Initializes demographics matrices
    N_matrix = np.zeros((I, S, T+S))
    Nhat = np.zeros((I, S, T+S))
    all_FertilityRates = np.zeros((I, S, f_range+T+S))
    FertilityRates = np.zeros((I, S, T+S))
    MortalityRates = np.zeros((I, S, T+S))
    ImmigrationRates = np.zeros((I, S, T+S))
    Migrants = np.zeros((I, S, T+S))
    g_N = np.zeros(T+S)
    lbar = np.zeros(T+S)
    #QUESTION HOW TO INTERPRET N_temp???
    N_temp = np.ones((I, S))/(I*S)
	
    #Gathers demographic data from the csv files for each I countries in countrynames
    for i in range(I):

        #If we want the countries to have unique demographic data, have the same index
        if DiffDemog:
            index = I_all.index(I_touse[i])
            if I > len(I_all):
                sys.exit("ERROR!!! We can't have more than", len(I_all), "countries without unique data. Change either parameter I so it is less than", len(I_all), " or change DiffDemog to False")
        #If we want the countries to have identical demographic data, have the same index
        else:
            index = 0

        if UseStaggeredAges:
            N_matrix[i,:,0] = np.loadtxt(("Data_Files/population.csv"),delimiter=',',skiprows=1, usecols=[index+1])[agestopull]*1000
            all_FertilityRates[i,FirstFertilityAge:LastFertilityAge+1,:f_range+T_1] = \
                        np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_fertility.csv"),delimiter=',',skiprows=1\
                        , usecols=(agestopull[FirstFertilityAge:LastFertilityAge+1]-22))[48-f_range:48+T_1,:])
            MortalityRates[i,FirstDyingAge:,:T_1] = np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_mortality.csv"),delimiter=',',skiprows=1, usecols=(agestopull[FirstDyingAge:]-67))[:T_1,:])
            Migrants[i,:MaxImmigrantAge,:T_1] = np.einsum("s,t->st",np.loadtxt(("Data_Files/net_migration.csv"),delimiter=',',skiprows=1, usecols=[index+1])[agestopull[:MaxImmigrantAge]]*100, np.ones(T_1))
        else:
            N_matrix[i,:,0] = np.loadtxt(("Data_Files/population.csv"),delimiter=',',skiprows=1, usecols=[index+1])[:S]*1000
            all_FertilityRates[i,FirstFertilityAge:LastFertilityAge+1,:f_range+T_1] = np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_fertility.csv"),delimiter=',',skiprows=1, usecols=range(1,f_range+1))[48-f_range:48+T_1,:])
            MortalityRates[i,FirstDyingAge:-1,:T_1] = np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_mortality.csv"),delimiter=',',skiprows=1, usecols=range(1,S-FirstDyingAge))[:T_1,:])
            Migrants[i,:MaxImmigrantAge,:T_1] = np.einsum("s,t->st",np.loadtxt(("Data_Files/net_migration.csv"),delimiter=',',skiprows=1, usecols=[index+1])[:MaxImmigrantAge]*100, np.ones(T_1))
		
        print "Got demographics for", I_all[index]

    #Gets initial population share
    Nhat[:,:,0] = N_matrix[:,:,0]/np.sum(N_matrix[:,:,0])

    #The last generation dies with probability 1
    MortalityRates[:,-1,:] = np.ones((I, T+S))

    #Gets steady-state values for all countries by taking the mean at year T_1-1 across countries
    f_bar = np.mean(all_FertilityRates[:,:,f_range+T_1-1], axis=0)
    rho_bar = np.mean(MortalityRates[:,:,T_1-1], axis=0)

    #Set to the steady state for every year beyond year T_1
    all_FertilityRates[:,:,f_range+T_1:] = np.tile(np.expand_dims(f_bar, axis=2), (I,1,T-T_1+S))
    MortalityRates[:,:,T_1:] = np.tile(np.expand_dims(rho_bar, axis=2), (I,1,T-T_1+S))

    #FertilityRates is exactly like all_FertilityRates except it begins at time t=0 rather than time t=-f_range
    FertilityRates[:,FirstFertilityAge:LastFertilityAge+1,:] = all_FertilityRates[:,FirstFertilityAge:LastFertilityAge+1,f_range:]

    #Gets initial world population growth rate
    g_N[0] = 0.

    #Calculates population numbers for each country
    for t in range(1,T+S):
        #Gets the total number of children and and percentage of children and stores them in generation 0 of their respective matrices
        #See equation 2.1
        N_matrix[:,0,t] = np.sum((N_matrix[:,:,t-1]*FertilityRates[:,:,t-1]), axis=1)
        N_temp[:,0] = np.sum((Nhat[:,:,t-1]*FertilityRates[:,:,t-1]), axis=1)

        #Finds the immigration rate for each year
        if t <= T_1:
            ImmigrationRates[:,:,t-1] = Migrants[:,:,t-1]/N_matrix[:,:,t-1]

        else:
            ImmigrationRates[:,:,t-1] = np.mean(ImmigrationRates[:,:,T_1-1], axis=0)

        #Gets the population distribution for the next year, taking into account immigration and mortality
        #See equation 2.2
        N_matrix[:,1:,t] = N_matrix[:,:-1,t-1]*(1+ImmigrationRates[:,:-1,t-1]-MortalityRates[:,:-1,t-1])
        N_temp[:,1:] = Nhat[:,:-1,t-1]*(1+ImmigrationRates[:,:-1,t-1]-MortalityRates[:,:-1,t-1])

        #Gets the population share
        Nhat[:,:,t] = N_matrix[:,:,t]/np.sum(N_matrix[:,:,t])

        #Gets the growth rate for the next year
        g_N[t] = np.sum(N_temp[:,:])-1

    ImmigrationRates[:,:,t] = Migrants[:,:,t]/N_matrix[:,:,t]

    pop_old = N_matrix[:,:,-1]
    pop_new = N_matrix[:,:,-1]

    iter = 0

    while np.max(np.abs(Nhat[:,:,-1] - Nhat[:,:,-2])) > tol:
        pop_new[:,0] = np.sum((pop_old[:,:]*FertilityRates[:,:,-1]),axis=1)
        pop_new[:,1:] = pop_old[:,:-1]*(1+ImmigrationRates[:,:-1,-1]-MortalityRates[:,:-1,-1])
        Nhat = np.dstack((Nhat,pop_new/np.sum(pop_new)))
        iter+=1

    print "The SS Population Share converged in", iter, "years beyond year T"

    Nhat_ss = Nhat[:,:,-1]

    if Graphs:
        for i in range(I):
            plt.plot(range(T+S), np.sum(Nhat[i,:,:T+S], axis=0))
        plt.legend(I_touse)
        plt.title("World Super Steady State Total Population")
        plt.show()
        plt.clf()

        for i in range(I):
            plt.plot(range(MaxImmigrantAge), ImmigrationRates[i,:MaxImmigrantAge,0])
        plt.legend(I_touse)
        plt.title("ImmigrationRates")
        plt.show()
        plt.clf()

        for i in range(I):
            plt.plot(range(FirstFertilityAge,LastFertilityAge+1), FertilityRates[i,FirstFertilityAge:LastFertilityAge+1,0])
        plt.legend(I_touse)
        plt.title("FertilityRates")
        plt.show()
        plt.clf()

    	
        for i in range(I):
            plt.plot(range(FirstDyingAge, S), MortalityRates[i,FirstDyingAge:,0])
        plt.legend(I_touse)
        plt.title("MortalityRates")
        plt.show()
        plt.clf()

        plotDemographics((S,T), range(I), [T+S], Nhat, I_touse)
    
    #Gets labor endowment per household. For now it grows at a constant rate g_A
    lbar[:T] = np.cumsum(np.ones(T)*g_A)
    lbar[T:] = np.ones(S)

    return MortalityRates, Nhat[:,:,:T+S], Nhat_ss

def plotDemographics(params, indexes, years, Nhat, countrynames):
    """
    Description:
        -Plots the population distribution of a given country for any number of specified years

    Inputs:
        params       = Tuple, Contains parameters S, T
        S            = Scalar in [10,80], Number of cohorts
        T            = Scalar >0, Number of years away from time t=0 until we reach the steady-state
        indexes      = List, Contains all indexes of countries to plot. Typically either range(I) or a single-item list (like [0])
        years        = List, Contains all indexes of years to plot. Typically either range(T+S) or a list like [0,20,100]
        Nhat         = [I,S,T+S+?] Matrix, Contains the population share, 
                                           including the years until it converges to the steady-state
        countrynames = List, Contains the names of each country to be plotted. Used only for the plot legend

    Functions called:
        -None

    Objects in Function:
        -None

    Returns: None
    """

    S, T = params

    #If we want to plot total population across years...
    if len(years) == 0:
        for i in indexes:
            plt.plot(range(S), np.sum(Nhat[i,:,:], axis=0))
        plt.legend(countrynames)
    #If we want to compare years rather than countries...
    else:
        for yeartograph in years:
            for i in indexes:
            #Checks to make sure we haven't requested to plot a year past the max year
                if yeartograph <= Nhat.shape[2]:
                    plt.plot(range(S), Nhat[i,:,yeartograph])
                else:
                    print "\nERROR: WE HAVE ONLY SIMULATED UP TO THE YEAR", T
                    time.sleep(15)
        if len(years) == 1:
            plt.legend(countrynames)
        else:
            plt.legend(years)

    plt.title("Population Distribution")
    plt.show()
    plt.clf()

#STEADY STATE FUNCTIONS

def get_kd(assets, kf, Nhat):
    """
    Description: 
        -Calculates the amount of domestic capital that remains in a given country
        -Corresponds to equation 2.12

    Inputs:
        assets = [I,S+1] Matrix, Assets for given time period
        kf     = [I,] Vector, Domestic capital held by foreigners
        Nhat   = [I,S] Matrix, Population share for a given year

    Functions called:
        -None

    Objects in Function:
        kd = [I,] Vector, Capital that is used in the domestic country

    Returns kd:
    """

    kd = np.sum(assets[:,:-1]*Nhat, axis=1) - kf

    return kd

def get_n(params):
    """
    Description: 
        -Calculates the total labor supply for each country
        -Corresponds equation 2.13

    Inputs:
        e    = [I,S] or [I,S,T+S] Matrix, Labor productivities
        Nhat = [I,S] or [I,S,T+S] Matrix, Population share for a given year

    Functions called:
        -None

    Objects in Function:
        n = [I,] or [I,T+S] Matrix, Total labor supply
	
    Returns: n
    """

    e, Nhat = params
    n = np.sum(e*Nhat, axis=1)

    return n

def get_Y(params, kd, n):
    """
    Description: 
        -Calculates the output timepath
        -Corresponds to equation 2.14

    Inputs:
        params = Tuple: Contains the parameters alpha and A
        alpha  = Scalar in (0,1), Production share of capital
        A      = [I,] Vector, Technology level for each country
        kd     = [I,] Vector or [I,T+S] Matrix, Domestic held capital stock
        n      = [I,] Vector or [I,T+S] Matrix, Aggregate labor supplies

	Functions called:
        -None

    Objects in Function:
        Y = [I,] Vector or [I,S+T+1] Matrix, Aggregate output

    Returns: Y

    """
    alpha, A = params

    if kd.ndim == 1:
        Y = (kd**alpha) * ((A*n)**(1-alpha))
    elif kd.ndim == 2:
        Y = (kd**alpha) * (np.einsum("i,is->is", A, n)**(1-alpha))

    return Y

def get_r(alpha, Y, kd):
    """
    Description:
        -Calculates the steady state rental rates or rental rates timepath, depending on where this function is called
        -Corresponds to equation 2.15

    Inputs:
        alpha = Scalar in (0,1), Production share of capital
        Y     = [I,] or [T+S] Vector, Either steady state output or output timepath
        kd    = [I,] or [T+S] Vector, Either steady state domestic-owned capital or domestic-owned capital timepath

    Functions called:
        -None

    Objects in Function:
        r = [I,] or [T+S] Vector, Either steady state rental rate or rental rate timepath

    Returns: r
    """

    r = alpha * Y / kd

    return r

def get_w(alpha, Y, n):
    """
    Description:
        -Calculates the wage rate
        -Corresponds to equation 2.16

    Inputs:
        alpha = Scalar in (0,1), Production share of capital
        Y     = [I,] Vector or [I,S+T] Matrix, Aggregate output
        n     = [I,] Vector or [I,S+T] Matrix, Aggregate labor supplies

    Functions called:
        -None

    Objects in Function:
        w = [I,] Vector or [I,S+T] Matrix, Wage rates

    Returns: w
    """

    w = (1-alpha) * Y / n
    return w

def getBequests(params, assets):
    """
    Description:
        -Gets the value of the bequests given to each generation
        -Corresponds to equations 2.17-2.18

    Inputs:
        params            = Tuple that contains the parameters I, S, FirstFertilityAge, FirstDyingAge, Nhat_current, Mortality_current
        I                 = Scalar in [1,7], Number of countries
        S                 = Scalar in [10,80], Number of cohorts
        FirstFertilityAge = Scalar in (0,S), First age when agents bear children
        FirstDyingAge     = Scalar in (0,S), First age when agents die
        Nhat_current      = [I, S] Matrix, Population share for a given year
        Mortality_current = [I, S] Matrix, Mortality rates for a given year
        assets            = [I, S] Matrix, Number of assets per agent in each cohort

    Functions called:
        -None

    Objects in Function:
        bq                    = [I, S] Matrix, Number of bequests each agent receives
        BQ                    = [I,] Vector, Total assets of the people who died this year
        num_bequest_receivers = [I,] Vector, Number of people who are eligible to receive bequests this year
        bq_per_agent          = [I,] Vector, Number of bequests each agent receives

    Returns: bq
    """

    I, S, FirstFertilityAge, FirstDyingAge, Nhat_current, Mortality_current = params

    #Initializes bequests
    bq = np.zeros((I, S))

    #Gets the total assets of the people who died this year
    BQ = np.sum(assets[:,FirstDyingAge:-1]*Mortality_current[:,FirstDyingAge:]*Nhat_current[:,FirstDyingAge:], axis=1)

    #Distributes the total assets equally among the eligible population for each country
    #NOTE: This will likely change as we get a more complex function for distributing the bequests
    num_bequest_receivers = np.sum(Nhat_current[:,FirstFertilityAge:FirstDyingAge], axis=1)
    bq_per_agent = BQ/num_bequest_receivers
    bq[:,FirstFertilityAge:FirstDyingAge] = np.einsum("i,s->is", bq_per_agent, np.ones(FirstDyingAge-FirstFertilityAge))

    return bq

def get_cvecss(params, w, r, assets):
    """
    Description:
        -Gets the consumption vector in the steady state while disregarding time differences

    Inputs:
        params = Tuple that contains the parameters e, delta, bq, g_A
        e      = [I,S] Matrix, Labor productivities
        delta  = Scalar in (0,1), Depreciation rate
        bq     = [I,S] Matrix, Bequests given to each cohort
        g_A    = 
        w      =
        r      =
        assets = [I,S+1] Matrix, Assets for each cohort, along with remaining assets after period S

    Functions called:
        -None

    Objects in Function:
        c_vec = 

    Returns: c_vec
    """
    e, delta, bq, g_A = params

    c_vec = np.einsum("i, is -> is", w, e)\
          + np.einsum("i, is -> is",(1 + r - delta) , assets[:,:-1])\
          + bq - assets[:,1:]*np.exp(g_A)

    return c_vec

def check_feasible(kd, Y, w, r, c):
    """
    Description:
        -Checks to see if the matrices have infeasible values

    Inputs:
             STEADY STATE:                  or   TPI:
        kd = [I,] Vector,                   or   [I,T+S] Matrix,     Domestic-owned capital
        Y  = [I,] Vector,                   or   [I,T+S] Matrix,     Output
        w  = [I,] Vector,                   or   [I,T+S] Matrix,     Wage rates
        r  = [I,] Vector,                   or   [T+S,] Vector,      
             Rental rates for each country  or   Global interest rate timepath
        c  = [I,S] Matrix,                  or   [I,S,T+S] Matrix,   Consumption

    Functions called:
        -None

    Objects in Function:
        Feasible = Boolean, True as long as all the matrices contain all feasible values

    Returns: Feasible
    """
    Feasible = True

    if np.any(kd<0) or np.any(np.isnan(kd)):
        Feasible=False
        print "WARNING! INFEASABLE VALUE ENCOUNTERED IN kd!"
        print "The following coordinates have values less than 0:"
        print np.argwhere(kd<0)
        print "The following coordinates have nan values"
        print np.argwhere(np.isnan(kd))

    if np.any(Y<0) or np.any(np.isnan(Y)):
        Feasible=False
        print "WARNING! INFEASABLE VALUE ENCOUNTERED IN Y!"
        print "The following coordinates have values less than 0:"
        print np.argwhere(Y<0)
        print "The following coordinates have nan values"
        print np.argwhere(np.isnan(Y))

    if np.any(r<0) or np.any(np.isnan(r)):
        Feasible=False
        print "WARNING! INFEASABLE VALUE ENCOUNTERED IN r!"
        print "The following coordinates have values less than 0:"
        print np.argwhere(r<0)
        print "The following coordinates have nan values"
        print np.argwhere(np.isnan(r))

    if np.any(w<0) or np.any(np.isnan(w)):
        Feasible=False
        print "WARNING! INFEASABLE VALUE ENCOUNTERED IN w!"
        print "The following coordinates have values less than 0:"
        print np.argwhere(w<0)
        print "The following coordinates have nan values"
        print np.argwhere(np.isnan(w))

    if np.any(c<0) or np.any(np.isnan(c)):
        Feasible=False
        print "WARNING! INFEASABLE VALUE ENCOUNTERED IN c_vec!"
        print "The following coordinates have values less than 0:"
        print np.argwhere(c<0)
        print "The following coordinates have nan values"
        print np.argwhere(np.isnan(c))

    return Feasible

def SteadyStateSolution(guess, I, S, beta, sigma, delta, alpha, e_ss, A, FirstFertilityAge, FirstDyingAge, Nhat_ss, Mortality_ss, g_A, PrintEulErrors):
    """
    Description: 
        -This is the function that will be optimized by fsolve to find the steady state

    Inputs:
        guess             = [I*S,] Matrix, Contains guesses for assets and foreign capital held in the steady state
        I                 = Scalar in [1,7], Number of countries
        S                 = Scalar in [10,80], Number of cohorts
        beta              = Scalar in (0,1), Time preference
        sigma             = Scalar in (0,1), Intratemporal elasticity of substitution
        delta             = Scalar in (0,1), Depreciation rate
        alpha             = Scalar in (0,1), Capital share of production
        e_ss              = [I,S] Matrix, Labor productivities by country and cohort
        A                 = [I] Vector, Scales production
        FirstFertilityAge = Scalar in (0,S), First age when agents bear children
        FirstDyingAge     = Scalar in (0,S), First age when agents die         
        Nhat_ss           = [I,S] Matrix, Steady state population shares
        Mortality_ss      = [I,S] Matrix, Steady state mortality rates
        g_A               = Scalar >0, Technical growth rate
        PrintEulErrors    = Boolean, Prints euler errors if set to True

    Objects in Function:
        assets    = [I,S] Matrix, Asset path for each country in the steady state
        kf        = [I,] Vector, Foreign capital held by foreigners in each country in the steady state
        kd        = [I,] Vector, Capital for each country in the steady state
        n         = [I,] Vector, Labor supply for each country in the steady state
        Y         = [I,] Vector, Output for each country in the steady state
        r         = [I,] Vector, Rental rate for each country in the steady state
        w         = [I,] Vector, Wage for each country in the steady state
        bq        = [I,S] Matrix, Bequests given to each country and cohort in the steady state
        c_vec     = [I,S] Matrix, Consumption by cohort in each country in the steady state
        Euler_c   = [I,S-1] Matrix, Corresponds to (1.16)
        Euler_r   = [I,] Vector, Corresponds to (1.17)
        Euler_kf  = Scalar, Corresponds to (1.18)
        all_Euler = [I*S,] Matrix, Flattened array of all the euler errors

    Functions called:
        get_kd: Gets domestic-owned capital
        get_n: Gets labor supply
        get_Y: Gets output
        get_r: Gets rental rate
        get_w: Gets wages
        get_Bequests: Gets bequests
        get_cvecss: Gets consumption by cohort
        check_feasible: Checks for any infeasible values in the different variables

    Returns: all_Euler
    """

    #Takes a 1D guess of length I*S and reshapes it to match what the original input into the fsolve looked like since fsolve flattens numpy arrays
    guess = np.reshape(guess[:,np.newaxis], (I, S))

    #Appends a I-length vector of zeros on ends of assets to represent no assets when born and no assets when dead
    assets = np.column_stack((np.zeros(I), guess[:,:-1], np.zeros(I)))

    #Sets kf as the last element of the guess vector for each country
    kf = guess[:,-1]

    #Getting the other variables
    kd = get_kd(assets, kf, Nhat_ss)
    nparams = (e_ss, Nhat_ss)
    n = get_n(nparams)
    Yparams = (alpha, A)
    Y = get_Y(Yparams, kd, n)
    r = get_r(alpha, Y, kd)
    w = get_w(alpha, Y, n)
    bqparams = (I, S, FirstFertilityAge, FirstDyingAge, Nhat_ss, Mortality_ss)
    bq = getBequests(bqparams, assets)
    cparams = (e_ss, delta, bq, g_A)
    c_vec = get_cvecss(cparams, w, r, assets)

    Feasible = check_feasible(kd, Y, w, r, c_vec)

    if Feasible == False: #Punishes the the poor choice of negative values in the fsolve
        all_Euler=np.ones((I*S))*999.
        print "Punishing fsolve"
    else:
        #Gets Euler equations
        Euler_c = c_vec[:,:-1] ** (-sigma) - beta * (1-Mortality_ss[:,:-1])*(c_vec[:,1:]*np.exp(g_A)) ** (-sigma) * (1 + r[0] - delta)
        Euler_r = r[1:] - r[0]
        Euler_kf = np.sum(kf*np.sum(Nhat_ss, axis=1))

        #Makes a new 1D vector of length I*S that contains all the Euler equations
        all_Euler = np.append(np.append(np.ravel(Euler_c), np.ravel(Euler_r)), Euler_kf)

        if PrintEulErrors:
            print "Max Euler SS Errors:"
            print "Euler_c:", np.max(np.absolute(Euler_c)),
            print "      Euler_r:", np.max(np.absolute(Euler_r)),
            print "      Euler_kf:", np.max(np.absolute(Euler_kf))

    return all_Euler

def getSteadyState(params, assets_init, kf_init):
    """
    Description:
        -Gets the steady state values of assets and foreign-held capital
        -Uses steady state assets and foreign-held capital to solve for the rest of the system

    Inputs:
        params            = Tuple, Contains the parameters I, S, beta, sigma, delta, alpha, e_ss, A,
                                                           FirstFertilityAge, FirstDyingAge, Nhat_ss, 
                                                           Mortality_ss, g_A, PrintEulErrors
        I                 = Scalar in [1,7], Number of countries
        S                 = Scalar in [10,80], Number of cohorts
        beta              = Scalar in (0,1), Time preference
        sigma             = Scalar in (0,1), Intratemporal elasticity of substitution
        delta             = Scalar in (0,1), Depreciation rate
        alpha             = Scalar in (0,1), Capital share of production
        e_ss              = [I,S] Matrix, Labor productivities by country and cohort
        A                 = [I] Vector, Scales production
        FirstFertilityAge = Scalar in (0,S), First age when agents bear children
        FirstDyingAge     = Scalar in (0,S), First age when agents die         
        Nhat_ss           = [I,S] Matrix, Steady state population shares
        Mortality_ss      = [I,S] Matrix, Steady state mortality rates
        g_A               = Scalar >0, Technical growth rate
        PrintEulErrors    = Boolean, Prints euler errors if set to True
        assets_init       = [I,S-1] Matrix, Contains guess for stedy state assets held beginning from age 1 (since assets at age 0 are 0)
        kf_init           = [I,] Vector, Contains guess for steady state foreign-held capital

    Objects in Function:
        guess     = [I,S] Matrix, Simply assets_init and kf_init stacked together before sending it to the fsolve
        ss        = [I*S,] Vector, Flattened array for the steady-state solution
        assets_ss = [I,S] Matrix, Asset path for each country in the steady state
        kf_ss     = [I,] Vector, Foreign capital held by foreigners in each country in the steady state
        kd_ss     = [I,] Vector, Capital for each country in the steady state
        n_ss      = [I,] Vector, Labor supply for each country in the steady state
        Y_ss      = [I,] Vector, Output for each country in the steady state
        r_ss      = Scalar >0, Global rental rate in the steady state
        w_ss      = [I,] Vector, Wage for each country in the steady state
        bq_ss     = [I,S] Matrix, Bequests given to each country and cohort in the steady state
        c_vec_ss  = [I,S] Matrix, Consumption by cohort in each country in the steady state

    Functions called:
        get_kd: Gets domestic-owned capital
        get_n: Gets labor supply
        get_Y: Gets output
        get_r: Gets rental rate
        get_w: Gets wages
        get_Bequests: Gets bequests
        get_cvecss: Gets consumption by cohort

    Returns: assets_ss, kf_ss, kd_ss, n_ss, Y_ss, r_ss, w_ss, c_vec_ss
    """
    I, S, beta, sigma, delta, alpha, e_ss, A, FirstFertilityAge, FirstDyingAge, Nhat_ss, Mortality_ss, g_A, PrintEulErrors = params

    #Merges the assets and kf together into one matrix that can be inputted into the fsolve function
    guess = np.column_stack((assets_init, kf_init))

    #Solves for the steady state
    ss = opt.fsolve(SteadyStateSolution, guess, args=params)

    #Reshapes the ss code
    ss = np.array(np.split(ss, I))

    #Breaks down the steady state matrix into the two separate assets and kf matrices.
    assets_ss = np.column_stack((np.zeros(I), ss[:,:-1], np.zeros(I)))
    kf_ss = ss[:,-1]

    #Gets the other steady-state values using assets and kf
    kd_ss = get_kd(assets_ss, kf_ss, Nhat_ss)
    nparams = (e_ss, Nhat_ss)
    n_ss = get_n(nparams)
    Yparams = (alpha, A)
    Y_ss = get_Y(Yparams, kd_ss, n_ss)
    #Because of the euler conditions in the fsolve, r_ss will be the same regardless of which country we use to calculate it
    r_ss = get_r(alpha, Y_ss[0], kd_ss[0])
    w_ss = get_w(alpha, Y_ss, n_ss)
    bqparams = (I, S, FirstFertilityAge, FirstDyingAge, Nhat_ss, Mortality_ss)
    bq_ss = getBequests(bqparams, assets_ss)
    c_vec_ss = np.einsum("i, is -> is", w_ss, e_ss)\
              + (1 + r_ss - delta)*assets_ss[:,:-1]\
              + bq_ss - assets_ss[:,1:]*np.exp(g_A)

    print "\nSteady State Found!\n"

    return assets_ss, kf_ss, kd_ss, n_ss, Y_ss, r_ss, w_ss, c_vec_ss

#TIMEPATH FUNCTIONS

def get_initialguesses(params, assets_ss, kf_ss, w_ss, r_ss, PrintLoc):
    """
    Description:
        With the parameters and steady state values, this function creates
        initial guesses in a linear path.

    Inputs:
        -Params (Tuple): Tuple of parameters from Main.py
        -Assets_ss[I,S+1,T+S]: Steady state assets value
        -kf_ss[I,]: Steady State value of foreign owned domestic capital
        -w_ss[I,]: Steady state value of wages
        -r_ss[I,]: Steady state value of rental rate

    Objects in Function:
        -othervariable_params (Tuple): A tuple specifically made for GetOtherVariables


    Outputs:
        -assets_init[I,]: Initial Asset path
        -kf_init[I,]: New initial foreign held capital
        -w_initguess[I,T+S]: Initial guess wage timepath
        -r_initguess[I,T+S]: Initial guess rental rate timepath
        -k_init[I,]: total capital stock initial guess
        -n_init[I,]: total labor initial guess
        -y_init[I,]: output labor initial guess
        -c_init[I,]: consumption initial guess
    """

    if PrintLoc: print "Getting initial guesses"
    #Unpacks parameters
    I, S, T, delta, alpha, e_init, A, FirstFertilityAge, FirstDyingAge, Nhat_init, Mortality_init, g_A = params

    #Sets initial assets and kf, start with something close to the steady state
    assets_init = assets_ss*.9
    kf_init = kf_ss*0

    wpath_guess = np.zeros((I, T+S))
    rpath_guess = np.zeros((T+S))

    #Gets initial kd, n, y, r, w, and K
    kd_init = get_kd(assets_init, kf_init, Nhat_init)
    nparams = e_init, Nhat_init
    n_init = get_n(nparams)
    Yparams = (alpha, A)
    Y_init = get_Y(Yparams, kd_init, n_init)
    r_init = get_r(alpha, Y_init[0], kd_init[0])
    w_init = get_w(alpha, Y_init, n_init)

    #Gets initial guess for rental rate path. This is set up to be parabolic.
    bb = -2 * (r_init-r_ss)/(T-1)
    cc = r_init
    aa = -bb / (2*(T-1))
    rpath_guess[:T] = aa * np.arange(0,T)**2 + bb*np.arange(0,T) + cc
    rpath_guess[T:] = r_ss

    #Gets initial guess for wage rate path. This is set up to be parabolic.
    bb = -2 * (w_init-w_ss)/(T-1)
    cc = w_init
    aa = -bb / (2*(T-1))
    wpath_guess[:,:T] = np.einsum("i,it->it", aa, np.tile(np.arange(0,T), (I,1))**2)\
                      + np.einsum("i,it->it", bb, np.tile(np.arange(0,T), (I,1)))\
                      + np.einsum("i,it->it", cc, np.ones((I,T)))
    wpath_guess[:,T:] = np.einsum("i,it->it", w_ss, np.ones((I,S)))

    return assets_init, kf_init, wpath_guess, rpath_guess

def get_foreignK_path(params, Kpath, rpath, kf_ss, PrintLoc):
    """
    Description:
       This calculates the timepath of the foreign capital stock. This is based on equation (1.12 and 1.13).
    Inputs:
        apath: Asset path, from our calculations
        rpath: Rental Rate path, also from our calculation
        
    Objects in Function:
        kdpath[I,S+T+1]: Path of domestic owned capital
        n[I,S+T+1]: Path of total labor
        kf_ss[I,]: Calculated from the steady state. 
        A[I,]: Parameters from above

    Outputs:
        kfPath[I,S+T+1]: Path of domestic capital held by foreigners.
    """
    if PrintLoc: print "Entering get_foreignK_path"

    I, S, T, alpha, e, A, Nhat = params

    #Sums the labor productivities across cohorts
    n = get_n((e, Nhat))

    #Declares the array that will later be used.
    kfPath = np.zeros((I,S+T))
    kdPath = np.zeros((I,S+T))

    #Gets the domestic-owned capital stock for each country except for the first country
    kdPath[1:,:] = (rpath/alpha)**(1/(alpha-1))*np.einsum("i,is->is", A[1:], n[1:,:])

    #This is using equation 1.13 solved for the foreign capital stock to caluclate the foreign capital stock
    #For everyone except the first country
    kfPath[1:,:] = Kpath[1:,:]-kdPath[1:,:]

    #To satisfy 1.18, the first country's assets is the negative of the sum of all the other countries' assets
    kfPath[0,:] = -np.sum(kfPath[1:,:],axis=0)

    #Making every year beyond t equal to the steady-state
    kfPath[:,T:] = np.einsum("i,s->is", kf_ss, np.ones(S))
        
    if PrintLoc: print "Leaving get_foreignK_path"
    return kfPath

def get_lifetime_decisions(params, c_1, wpath_chunk, rpath_chunk, e_chunk, mortality_chunk, starting_assets, bq, current_age):
    """
    Description:
        This solves for equations 1.15 and 1.16 in the StepbyStep pdf for a certain generation
    Inputs:
        -c_1: Initial consumption (not necessarily for the year they were born)
        -wpath_chunk: Wages of an agents lifetime, a section of the timepath
        -rpath_chunk: Rental rate of an agents lifetime, a section of the timepath
        -e_chunk: Worker productivities of an agents lifetime, a section of the global matrix
        -starting_assets: Initial assets of the agent. Will be 0s if we are beginning in the year the agent was born
        -current_age: Current age of the agent

        Objects in Function:
            -NONE

    Outputs:
        -c_path[I, S]: Path of consumption until the agent dies
        -asset_path[I, S+1]: Path of assets until the agent dies
    """

    I, S, beta, sigma, delta, g_A = params

    num_decisions = S-current_age-1# -1 Because we already have (or have guessed) our starting consumption

    #Initializes the cpath and asset path vectors
    c_path = np.zeros((I, num_decisions+1))
    asset_path = np.zeros((I, num_decisions+2))

    #For each country, the cpath and asset path vectors' are the initial values provided.
    c_path[:,0] = c_1
    asset_path[:,0] = starting_assets

    #Based on the individual chunks, these are the households choices
    for p in range(1,num_decisions+1):
        c_path[:,p] = ((beta * (1-mortality_chunk[:,p-1]) * (1 + rpath_chunk[p] - delta))**(1/sigma) * c_path[:,p-1])*np.exp(-g_A)
        asset_path[:,p] = (wpath_chunk[:,p-1]*e_chunk[:,p-1] + (1 + rpath_chunk[p-1] - delta)*asset_path[:,p-1] + bq[:,p-1] - c_path[:,p-1])*np.exp(-g_A)
	
    asset_path[:,p+1] = (wpath_chunk[:,p]*e_chunk[:,p] + (1 + rpath_chunk[p] - delta)*asset_path[:,p] - c_path[:,p])*np.exp(-g_A)

    return c_path, asset_path

def find_optimal_starting_consumptions(c_1, wpath_chunk, rpath_chunk, e_chunk, mortality_chunk, starting_assets, bq, current_age, params):
    """
    Description:
       Takes the assets path from the get_householdchoices_path function and creates Euluer errors

    Inputs:
    Dimension varies
        -c_1: Initial consumption (not necessarily for the year they were born)
        -wpath_chunk: Wages of an agents lifetime, a part of the timepath
        -rpath_chunk: Rental rate of an agents lifetime, another part of the timepath.
        -epath_chunk: Worker productivities of an agents lifetime, another part.
        -starting_assets: Initial assets of the agent. It's 0 at the beginning of life.
        -current_age: Current age of the agent


    Objects in Function:
        -cpath: Path of consumption based on chunk given.
        -assets_path: Path of assets based on the chunks given

    Outputs:
        -Euler:A flattened version of the assets_path matrix
    """

    #Executes the get_household_choices_path function. Sees above.
    c_path, assets_path = get_lifetime_decisions(params, c_1, wpath_chunk, rpath_chunk, e_chunk, mortality_chunk, starting_assets, bq, current_age)
    Euler = np.ravel(assets_path[:,-1])

    if np.any(c_path<0):
        #print "WARNING! The fsolve for initial optimal kids consumption guessed a negative number"
        Euler=np.ones(Euler.shape[0])*9999.

    #print "Max Euler Error in find_optimal_starting_consumptionsK", np.max(np.absolute(Euler))

    return Euler

def get_cons_assets_matrix(params, wpath, rpath, starting_assets, PrintLoc, Print_cabqTimepaths):
    if PrintLoc: print "Entering get_cons_assets_matrix"

    I, S, T, T_1, beta, sigma, delta, e, FirstFertilityAge, FirstDyingAge, Nhat, MortalityRates, g_A = params

    #Initializes timepath variables
    c_timepath = np.zeros((I,S,S+T))
    a_timepath = np.zeros((I, S+1, S+T))
    a_timepath[:,:,0]=starting_assets
    bq_timepath = np.zeros((I, S, S+T))

    c_timepath[:,S-1,0] = wpath[:,0]*e[:,S-1,0] + (1 + rpath[0] - delta)*a_timepath[:,S-1,0]

    household_params = I, S, beta, sigma, delta, g_A

    if Print_cabqTimepaths:
        print "Initial matrices"
        print "Consumption"
        print np.round(np.transpose(c_timepath[0,:,:2]), decimals=3)
        print "Assets"
        print np.round(np.transpose(a_timepath[0,:,:2]), decimals=3)
        print "Bequests"
        print np.round(np.transpose(bq_timepath[0,:,:2]), decimals=3)

    #Fills the upper triangle (including the main diagonal) by iterating by number of periods until death
    if PrintLoc: print "Entering upper triangle loop"
    for p in range(1,S):
        #We are only doing this for all generations alive in time t=0
        t = 0

        #Getting the current age of the agent
        current_age = S-p-1

        #Uses the previous generation's consumption at age s to get the value for our guess
        c_guess = (c_timepath[:,current_age+1,t]/((beta*(1+rpath[t]-delta))**(1/sigma)))/np.exp(g_A)

        agent_assets = starting_assets[:,current_age]

        #Gets the bequests this agent will recieve in his remaining lifetime
        agent_bq = np.diagonal(bq_timepath[:,current_age:,t:t+p+1], axis1=1, axis2=2)

        #Gets labor productivities this agent will recieve in his remaining lifetime
        agent_e = np.diagonal(e[:,current_age:,t:t+p+1], axis1=1, axis2=2)

        agent_mortality = np.diagonal(MortalityRates[:,current_age:,t:t+p+1], axis1=1, axis2=2)

        #Gets optimal initial consumption beginning in the current age of the agent using chunks of w and r that span the lifetime of the given generation
        opt_consump = opt.fsolve(find_optimal_starting_consumptions, c_guess, args = \
            (wpath[:,t:t+p+1], rpath[t:t+p+2], agent_e, agent_mortality, agent_assets, agent_bq, current_age, household_params))

        #Gets optimal timepaths beginning initial consumption and starting assets
        cpath_indiv, apath_indiv = get_lifetime_decisions\
            (household_params, opt_consump, wpath[:,t:t+p+1], rpath[t:t+p+2], agent_e, agent_mortality, agent_assets, agent_bq, current_age)

        for i in xrange(I):
            np.fill_diagonal(c_timepath[i,current_age:,:], cpath_indiv[i,:])
            np.fill_diagonal(a_timepath[i,current_age:,:], apath_indiv[i,:])

        bq_params = (I, S, FirstFertilityAge, FirstDyingAge, Nhat[:,:,p-1], MortalityRates[:,:,p-1])
        bq_timepath[:,:,p-1] = getBequests(bq_params, a_timepath[:,:,p-1])

        if Print_cabqTimepaths:
            print "p =", p, "current_age =", current_age
            print "Consumption"
            print np.round(np.transpose(c_timepath[0,:,:p+2]), decimals=4)
            print "c_guess", np.round(c_guess, decimals=4)
            print "Assets"
            print np.round(np.transpose(a_timepath[0,:,:p+2]), decimals=4)
            print "Bequests"
            print np.round(np.transpose(bq_timepath[0,:,:p+2]), decimals=4)
            print "agent_bq", np.round(agent_bq[0,:], decimals=4)

    if PrintLoc: print "Entering non-upper triangle loop"
    #Fills everything except for the upper triangle (excluding the main diagonal)
    for t in xrange(1,T):
        current_age = 0
        p = S-current_age-1

        agent_assets = np.zeros((I))

        #Uses the previous generation's consumption at age s to get the value for our guess
        c_guess = c_timepath[:,current_age,t-1]

        #Gets the bequests this agent will recieve in his remaining lifetime
        agent_bq = np.diagonal(bq_timepath[:,current_age:,t:t+p+1], axis1=1, axis2=2)

        #Gets labor productivities this agent will recieve in his remaining lifetime
        agent_e = np.diagonal(e[:,current_age:,t:t+p+1], axis1=1, axis2=2)

        agent_mortality = np.diagonal(MortalityRates[:,current_age:,t:t+p+1], axis1=1, axis2=2)

        opt_consump = opt.fsolve(find_optimal_starting_consumptions, c_guess, args = \
            (wpath[:,t:t+p+1], rpath[t:t+p+2], agent_e, agent_mortality, agent_assets, agent_bq, current_age, household_params))

        #Gets optimal timepaths beginning initial consumption and starting assets
        cpath_indiv, apath_indiv = get_lifetime_decisions\
            (household_params, opt_consump, wpath[:,t:t+p+1], rpath[t:t+p+2], agent_e, agent_mortality, agent_assets, agent_bq, current_age)

        for i in range(I):
            np.fill_diagonal(c_timepath[i,:,t:], cpath_indiv[i,:])
            np.fill_diagonal(a_timepath[i,:,t:], apath_indiv[i,:])

        if t >= T_1:
            temp_t = T_1
        else:
            temp_t = t

        bq_params = (I, S, FirstFertilityAge, FirstDyingAge, Nhat[:,:,temp_t+S-2], MortalityRates[:,:,temp_t+S-2])
        bq_timepath[:,:,t+S-2] = getBequests(bq_params, a_timepath[:,:,temp_t+S-2])

        if Print_cabqTimepaths:
            print "t = ", t, "current_age =", current_age
            print "Consumption"
            print np.round(np.transpose(c_timepath[0,:,:p+2]), decimals=4)
            print "c_guess", np.round(c_guess, decimals=4)
            print "Assets"
            print np.round(np.transpose(a_timepath[0,:,:p+2]), decimals=4)
            print "Bequests"
            print np.round(np.transpose(bq_timepath[0,:,:p+2]), decimals=4)
            print "agent_bq", np.round(agent_bq[0,:], decimals=4)


    if PrintLoc: print "Leaving get_cons_assets_matrix"
    return c_timepath, a_timepath

def get_wpathnew_rpathnew(params, wpath, rpath, starting_assets, kd_ss, kf_ss, w_ss, r_ss, PrintLoc, Print_cabqTimepaths, UseTape):
    """
    Description:
        Takes initial paths of wages and rental rates, gives the consumption path and the the wage and rental paths that are implied by that consumption path.

    Inputs:
        -w_path0[I, S+T+1]: initial w path
        -r_path0[I, S+T+1]: initial r path

    Objects in Function:
    Note that these vary in dimension depending on the loop.
        -current_age: The age of the cohort at time 0
        -opt_consump: Solved for consumption
        -starting_assets: Initial assets for the cohorts. 
        -cpath_indiv: The small chunk of cpath.
        -assetpath_indiv: The small chunk of assetpath_indiv
        -optimalconsumption: Solved from the chunks
        -c_timepath: Overall consumption path
        -a_timepath: Overall assets timepath
        -kfpath: Foreign held domestic capital
        -agent assets: Assets held by individuals.

    Outputs:
        -w_path1[I,S+T+1]: calculated w path
        -r_path1[I,S+T+1]: calculated r path
        -CPath[I,S+T+1]: Calculated aggregate consumption path for each country
        -Kpath[I,S+T+1]: Calculated capital stock path.
        -Ypath1[I, S+T+1]: timepath of assets implied from initial guess

    """
    if PrintLoc: print "Entering get_wpathnew_rpathnew"

    I, S, T, T_1, beta, sigma, delta, alpha, e, A, FirstFertilityAge, FirstDyingAge, Nhat, MortalityRates, g_A = params

    ca_params = (I, S, T, T_1, beta, sigma, delta, e, FirstFertilityAge, FirstDyingAge, Nhat, MortalityRates, g_A)
    c_timepath, a_timepath = get_cons_assets_matrix(ca_params, wpath, rpath, starting_assets, PrintLoc, Print_cabqTimepaths)

    #Calculates the total amount of capital in each country
    Kpath=np.sum(a_timepath[:,:-1,:]*Nhat, axis=1)

    #Calculates Aggregate Consumption
    Cpath=np.sum(c_timepath, axis=1)

    #After time period T, the total capital stock and total consumption is forced to be the steady state
    Kpath[:,T:] = np.einsum("i,t->it", kd_ss+kf_ss, np.ones(S))
    Cpath[:,T:] = np.einsum("i,t->it", Cpath[:,T-1], np.ones(S))

    #Gets the foriegned owned capital
    kf_params = (I, S, T, alpha, e, A, Nhat)
    kfpath = get_foreignK_path(kf_params, Kpath, rpath, kf_ss, PrintLoc)

    #Based on the overall capital path and the foreign owned capital path, we get new w and r paths.
    kdpath = Kpath - kfpath

    if UseTape:
        kdpath_with_tape = np.clip(kdpath, 0.001, np.max(kdpath))
    else:
        kdpath_with_tape = kdpath

    nparams = (e, Nhat)
    npath = get_n(nparams)
    Yparams = (alpha, A)
    Ypath = get_Y(Yparams, kdpath_with_tape, npath)
    rpath_new = get_r(alpha, Ypath[0], kdpath_with_tape[0])
    wpath_new = get_w(alpha, Ypath, npath)

    #Checks to see if any of the timepaths have negative values or nans
    Feasible = check_feasible(kdpath, Ypath, wpath, rpath, c_timepath)

    if PrintLoc: print "Leaving get_wpathnew_rpathnew"
    return wpath_new, rpath_new, Cpath, Kpath, Ypath

def get_Timepath(params, wstart, rstart, assets_init, kd_ss, kf_ss, w_ss, r_ss, PrintLoc, Print_cabqTimepaths, UseTape):

    I, S, T, T_1, beta, sigma, delta, alpha, e, A, FirstFertilityAge, FirstDyingAge, Nhat, MortalityRates, g_A, distance, diff, xi, MaxIters = params

    Iter=1 #Serves as the iteration counter
    wr_params = (I, S, T, T_1, beta, sigma, delta, alpha, e, A, FirstFertilityAge, FirstDyingAge, Nhat, MortalityRates, g_A)

    while distance>diff and Iter<MaxIters: #The timepath iteration runs until the distance gets below a threshold or the iterations hit the maximum
            wpath_new, rpath_new, Cpath_new, Kpath_new, Ypath_new = \
            get_wpathnew_rpathnew(wr_params, wstart, rstart, assets_init, kd_ss, kf_ss, w_ss, r_ss, PrintLoc, Print_cabqTimepaths, UseTape)
            try:
                dist1=sp.linalg.norm(wstart-wpath_new,2) #Norm of the wage path
                dist2=sp.linalg.norm(rstart-rpath_new,2) #Norm of the intrest rate path
                distance=max([dist1,dist2]) #We take the maximum of the two norms to get the distance
                print "Iteration:",Iter,", Norm Distance: ", distance#, "Euler Error, ", EError
                #print "Ypath"
                #print Ypath_new[0,:]
                #print "rpath"
                #print rpath_new[0,:]
            except:
                distance = diff+333
                print "Iteration:",Iter,", Error in calculating the distance"
                sys.exit("\nERROR!!: We are getting nowhere. Take a look at those nans. Quitting the program\n")

            Iter+=1 #Updates the iteration counter
            if distance<diff or Iter==MaxIters: #When the distance gets below the tolerance or the maximum of iterations is hit, then the TPI finishes.
                wend=wpath_new
                rend=rpath_new
                Cend=Cpath_new
                Kend=Kpath_new
                Yend=Ypath_new
            if Iter==MaxIters: #In case it never gets below the tolerance, it will throw this warning and give the last timepath.
                print "Doesn't converge within the maximum number of iterations"
                print "Providing the last iteration"

            wstart=wstart*xi+(1-xi)*wpath_new #Convex conjugate of the wage path
            rstart=rstart*xi+(1-xi)*rpath_new #Convex conjugate of the intrest rate path

    return wend, rend, Cend, Kend, Yend

def plotTimepaths(I, S, T, sig, wpath, rpath, cpath, kpath, Ypath, I_touse,save,show):

    if save==True and show==True:
        print "Cannot save and show graphs at the same time, showing only!"
        save=False

    for i in xrange(I): #Wages
        plt.plot(np.arange(0,T),wpath[i,:T], label=I_touse[i])
    plt.title("Time path for Wages")
    plt.ylabel("Wages")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    if show: plt.show()
    if save: 
        name= "wages_"+str(I)+"_"+str(S)+"_"+str(sig)+".png"
        plt.savefig(name)
        plt.cla()

    #Rental Rates  
    plt.plot(np.arange(0,T),rpath[:T], label='Global Interest Rate')
    plt.title("Time path for Rental Rates")
    plt.ylabel("Rental Rates")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    if show: plt.show()
    if save: 
        name= "rentalrate_"+str(I)+"_"+str(S)+"_"+str(sig)+".png"
        plt.savefig(name)
        plt.cla()


    for i in xrange(I): #Aggregate Consumption
        plt.plot(np.arange(0,S+T),cpath[i,:],label=I_touse[i])
    plt.title("Time Path for Aggregate Consumption")
    plt.ylabel("Consumption Level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    if show: plt.show()
    if save: 
        name= "aconsump_"+str(I)+"_"+str(S)+"_"+str(sig)+".png"
        plt.savefig(name)
        plt.cla()


    for i in xrange(I): #Aggregate Capital Stock
        plt.plot(np.arange(0,T),kpath[i,:T],label=I_touse[i])
    plt.title("Time path for Capital Path")
    plt.ylabel("Capital Stock level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    if show: plt.show()
    if save: 
        name= "acapitalstock_"+str(I)+"_"+str(S)+"_"+str(sig)+".png"
        plt.savefig(name)
        plt.cla()


    for i in xrange(I):
        plt.plot(np.arange(0,T),Ypath[i,:T],label=I_touse[i])
    plt.title("Time path for Output")
    plt.ylabel("Output Stock level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    if show: plt.show()
    if save: 
        name= "aoutput_"+str(I)+"_"+str(S)+"_"+str(sig)+".png"
        plt.savefig(name)
        plt.cla()

