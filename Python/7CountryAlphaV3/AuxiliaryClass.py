from __future__ import division
import csv
import time
import numpy as np
import scipy as sp
import scipy.optimize as opt
from matplotlib import pyplot as plt

import AuxiliaryDemographics as demog

class OLG(object):
    """
    TODO: SOME SORT OF DESCRIPTION OF WHAT THE OLG CLASS DOES

    For each function there are the following categories:
        Description:                    Brief description of what the function does
        Inputs:                         Lists the inputs that the function uses
        Variables Called From Object:   Lists the variables that the function calls from storage
        Variables Stored in Object:     Lists the variables that are put into storage
        Other Functions Called:         Lists the other non-library functions needed to complete the process of the current function
        Objects in Function:            Lists the variables that are exclusive to that function
        Outputs:                        Lists the outputs that the function puts out.

    Note that if a category isn't listed, then there aren't any variables/functions that fit that category in that
        particular function.
    """

    def __init__(self, countries, HH_Params, Firm_Params, Lever_Params):
        """
        Description: 
            -This creates the object and stores all of the parameters into the object.
            -The initialization is the starting point for model.

        Inputs:
            -self: "self" stores all of the components of the model. To access any part,
                 simply type "self.variable_name" while in the object and "objectname.variable_name"
                 outside the object. Every other object function will just take this as given, so 
                 future mentions of self won't be rewritten.

            -countries              = tuple: contains a dictionary and tuple for countries and their associated number
            -Firm_Params            = tuple: contains alpha, annualized delta, chi, rho and g_A
            -HH_Params              = tuple: contains S, I, annualized Beta and sigma.
            -Lever_Params           = tuple: contains the following boolean levers indicated by the users:
                                             PrintAges,self.PrintLoc,self.Print_cabqTimepaths,self.Iterate,self.UseDiffDemog,
                                             self.UseDiffProductivities,self.ADJUSTKOREAIMMIGRATION,self.VectorizeHouseholdSolver

        Variables Stored in Object:

            - self.A                = Array: [I], Technology level for each country
            - self.agestopull       = Array: [S], Contains which ages to be used from the data when S<80
            - self.e                = Array: [I,S,T], Labor Productivities
            - self.e_ss             = Array: [I,S], Labor produtivities for the Steady State
            - self.lbar             = Array: [T+S], Time endowment in each year

            - self.CheckerMode      = Boolean: Used in conjunction with Checker.py, a MPI code that checks the
                                               robustness of the code. With this activated, the code only prints
                                               the statements that are necessary. This speeds up the robust check
                                               process.
            - self.Iterate          = Boolean: Activates printing the iteration number and euler errors at each
                                                   step of the TPI process.
            - self.PrintAges        = Boolean: Prints the ages calculated in the demographics
            - self.PrintLoc         = Boolean: Prints the location of the code, used for debugging purposes
            - self.UseDiffDemog     = Boolean: Allows each country to have different demographics.
                
            - self.I_dict           = Dictionary: [I], Associates a country with a number
            - self.I_touse          = List: [I], Roster of countries that are being used

            - self.alpha            = Scalar: Capital share of production
            - self.beta             = Scalar: Calculated overall future discount rate
            - self.chi              = Scalar: Leisure preference Parameter
            - self.delta            = Scalar: Calulated overall depreciation rate
            - self.g_A              = Scalar: Growth rate of technology
            - self.rho              = Scalar: The intratemporal elasticity of substitution between consumption and leisure
            - self.sigma            = Scalar: Rate of Time Preference
            - self.FirstDyingAge    = Int: First age where mortality rates effect agents
            - self.FirstFertilityAge= Int: First age where agents give birth
            - self.I                = Int: Number of Countries
            - self.LastFertilityAge = Int: Last age where agents give birth
            - self.LeaveHouseAge    = Int: First age where agents don't count as children in utility function
            - self.MaxImmigrantAge  = Int: No immigration takes place for cohorts older than this age
            - self.S                = Int: Number of Cohorts
            - self.T                = Int: Number of time periods
            - self.T_1              = Int: Transition year for the demographics
            - self.Timepath_counter = Int: Counter that keeps track of the number of iterations in solving for the time paths
            - self.IterationsToShow = Set: A set of user inputs of iteration of TPI graphs to show

        Other Functions Called:
            - getkeyages = Gets the important ages for calculating demographic dynamics like FirstFertilityAge, etc. 
            - Importdata = Imports the demographic data from CSV files

        Objects in Function:
            - beta_annual           = Scalar: Original value for beta. Adjusted by S and stored as self.beta
            - delta_annual          = Scalar: Original value for delta. Adjusted by S and stored as self.delta
        """

        #PARAMETER SET UP

        #HH Parameters
        (self.S, self.I, beta_annual, self.sigma) = HH_Params
        
        self.beta=beta_annual**(70/self.S)

        self.T = int(round(4*self.S))

        self.T_1 = self.S

        if self.S > 50:
            self.T_1 = 50

        #Demographics Parameters
        self.I_dict, self.I_touse = countries

        #Firm Parameters
        (self.alpha,delta_annual,self.chi,self.rho, self.g_A)= Firm_Params
        self.delta=1-(1-delta_annual)**(70/self.S)

        #Lever Parameters
        (PrintAges,self.PrintLoc,self.CheckerMode,self.Iterate,self.UseDiffDemog,\
         self.UseDiffProductivities,self.ADJUSTKOREAIMMIGRATION,self.VectorizeHouseholdSolver) = Lever_Params

        #Getting key ages for calculating demographic dynamics
        self.LeaveHouseAge, self.FirstFertilityAge, self.LastFertilityAge,\
        self.MaxImmigrantAge, self.FirstDyingAge, self.agestopull = demog.getkeyages(self.S,PrintAges)

        if self.UseDiffDemog:
            self.A = np.ones(self.I)+np.cumsum(np.ones(self.I)*.05)-.05 #Techonological Change, used for when countries are different

        else:
            self.A = np.ones(self.I) #Techonological Change, used for idential countries

        #Initialize Labor Productivities
        if self.UseDiffProductivities:
            self.e = np.ones((self.I, self.S, self.T+self.S))
            self.e[:,self.FirstDyingAge:,:] = 0.3
            self.e[:,:self.LeaveHouseAge,:] = 0.3
        else:
            self.e = np.ones((self.I, self.S, self.T+self.S)) #Labor productivities

        self.e_ss=self.e[:,:,-1]

        #Initilize Time Endowment
        self.lbar = np.cumsum(np.ones(self.T+self.S)*self.g_A)
        self.lbar[self.T:] = np.ones(self.S)
        self.lbar[:self.T] = np.ones(self.T)
        self.lbar_ss=self.lbar[-1]

        #Imports all of the data from .CSV files needed for the model
        self.Import_Data()

        #Initialize counter that will keep track of the number of iterations the time path solver takes
        self.Timepath_counter = 1

    #DEMOGRAPHICS SET-UP

    def Import_Data(self):
        """
        Description:
            - This function activates importing the .CSV files that contain our demographics data

        Variables Called from Object:
            - self.agestopull             = Array: [S], Contains which ages to be used from the data when S<80
            - self.I                      = Int: Number of Countries
            - self.S                      = Int: Number of Cohorts
            - self.T                      = Int: Number of Time Periods
            - self.FirstFertilityAge      = Int: First age where agents give birth
            - self.LastFertilityAge       = Int: Last age where agents give birth
            - self.UseDiffDemog           = Boolean: True activates using unique country demographic data
            - self.PrintLoc               = Boolean: True prints the location of the code, used for debugging purposes
            - self.ADJUSTKOREAIMMIGRATION = Boolean: True will correctly adjust Korea's immigration, which is off by a factor of 100

        Variables Stored in Object:
            - self.N                      = Array: [I,S,T], Population of each country for each age cohort and year
            - self.Nhat                   = Array: [I,S,T], World opulation share of each country for each age cohort and year
            - self.all_FertilityAges      = Array: [I,S,f_range+T], Fertility rates from a f_range years ago to year T
            - self.FertilityRates         = Array: [I,S,T], Fertility rates from the present time to year T
            - self.MortalityRates         = Array: [I,S,T], Mortality rates of each country for each age cohort and year
            - self.Migrants               = Array: [I,S,T], Number of immigrants
            - self.g_N                    = Array: [T], Population growth rate for each year

        Other Functions Called:
            - None

        Objects in Function:
            - f_range                     = Int: Number of fertile years, will be used to correctly store the fertilty data
            - index                       = Int: Unique index for a given country that corresponds to the I_dict
            - f_bar                       = Array: [I,S], Average fertility rate across all countries and cohorts in year T_1, 
                                            used to get the SS demographics
            - rho_bar                     = Array: [I,S], Average mortality rate across all countries and cohorts in year T_1, 
                                            used to get the SS demographics

        Outputs:
            - None

        """

        f_range=self.LastFertilityAge+1-self.FirstFertilityAge


        self.N=np.zeros((self.I,self.S,self.T))
        self.Nhat=np.zeros((self.I,self.S,self.T))
        self.all_FertilityRates = np.zeros((self.I, self.S, f_range+self.T))
        self.FertilityRates = np.zeros((self.I, self.S, self.T))
        self.MortalityRates = np.zeros((self.I, self.S, self.T))
        self.Migrants = np.zeros((self.I, self.S, self.T))
        self.g_N = np.zeros(self.T)

        I_all = list(sorted(self.I_dict, key=self.I_dict.get))

        #We loop over each country to import its demographic data
        for i in xrange(self.I):

            #If the bool UseDiffDemog == True, we get the unique country index number for importing from the .CSVs
            if self.UseDiffDemog:
                index = self.I_dict[self.I_touse[i]]

            #Otherwise we just only use the data for one specific country
            else:
                index = 0

            #Importing the data and correctly storing it in our demographics matrices
            self.N[i,:,0] = np.loadtxt(("Data_Files/population.csv"),delimiter=',',\
                    skiprows=1, usecols=[index+1])[self.agestopull]*1000

            self.all_FertilityRates[i,self.FirstFertilityAge:self.LastFertilityAge+1,\
                    :f_range+self.T_1] =  np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_fertility.csv"),delimiter=',',skiprows=1\
                    ,usecols=(self.agestopull[self.FirstFertilityAge:self.LastFertilityAge+1]-22))[48-f_range:48+self.T_1,:])

            self.MortalityRates[i,self.FirstDyingAge:,:self.T_1] = np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_mortality.csv")\
                    ,delimiter=',',skiprows=1, usecols=(self.agestopull[self.FirstDyingAge:]-67))[:self.T_1,:])

            self.Migrants[i,:self.MaxImmigrantAge,:self.T_1] = np.einsum("s,t->st",np.loadtxt(("Data_Files/net_migration.csv"),delimiter=','\
                    ,skiprows=1, usecols=[index+1])[self.agestopull[:self.MaxImmigrantAge]]*100, np.ones(self.T_1))

            if self.PrintLoc: print "Got Demographics for", I_all[index]

            if self.ADJUSTKOREAIMMIGRATION and I_all[index] == "korea":
                self.Migrants[i,:]/=100

        #Gets initial population share
        self.Nhat[:,:,0] = self.N[:,:,0]/np.sum(self.N[:,:,0])

        #The last generation dies with probability 1
        self.MortalityRates[:,-1,:] = np.ones((self.I, self.T))

        #Gets steady-state values for all countries by taking the mean at year T_1-1 across countries
        f_bar = np.mean(self.all_FertilityRates[:,:,f_range+self.T_1-1], axis=0)
        rho_bar = np.mean(self.MortalityRates[:,:,self.T_1-1], axis=0)

        #Set to the steady state for every year beyond year T_1
        self.all_FertilityRates[:,:,f_range+self.T_1:] = np.tile(np.expand_dims(f_bar, axis=2), (self.I,1,self.T-self.T_1))
        self.MortalityRates[:,:,self.T_1:] = np.tile(np.expand_dims(rho_bar, axis=2), (self.I,1,self.T-self.T_1))

        #FertilityRates is exactly like all_FertilityRates except it begins at time t=0 rather than time t=-f_range
        self.FertilityRates[:,self.FirstFertilityAge:self.LastFertilityAge+1,:] = self.all_FertilityRates[:,self.FirstFertilityAge:self.LastFertilityAge+1,f_range:]

        #Gets initial world population growth rate
        self.g_N[0] = 0.

    def Demographics(self, demog_ss_tol, UseSSDemog=False):
        """
        Description:
            - This function calculates the population dynamics and steady state from the imported data by doing the following:
                1. For each year from now until year T, uses equations 3.11 and 3.12 to find the net population in a new year.
                    (Note however that after year T_1 the fertility, mortality, and immigration rates are set to be the same across countries)
                2. Divides the new population by the world population to get the population share for each country and cohort
                3. While doing steps 1-2, finds the immigration rate since the data only gives us net migration
                4. After getting the population dynamics until year T, we continue to get population shares of future years beyond time T 
                    as explained in steps 1-2 until it converges to a steady state
                5. Stores the new steady state and non-steady state variables of population shares and mortality in the OLG object

        Inputs:
            - demog_ss_tol              = Scalar: The tolerance for the greatest absolute difference between 2 years' population shares 
                                                  before it is considered to be the steady state
            - UseSSDemog                = Boolean: True uses the steady state demographics in calculating the transition path. Mostly used for debugging purposes

        Variables Called from Object:
            - self.I                    = Int: Number of Countries
            - self.S                    = Int: Number of Cohorts
            - self.T                    = Int: Number of Time Periods
            - self.T_1                  = Int: Transition year for the demographics

            - self.N                    = Array: [I,S,T], Population of each country for each age cohort and year
            - self.Nhat                 = Array: [I,S,T], World opulation share of each country for each age cohort and year
            - self.FertilityRates       = Array: [I,S,T], Fertility rates from the present time to year T
            - self.Migrants             = Array: [I,S,T], Number of immigrants
            - self.MortalityRates       = Array: [I,S,T], Mortality rates of each country for each age cohort and year

            - self.PrintLoc             = Boolean: True prints the location of the code, used for debugging purposes

        Variables Stored in Object:
            - self.ImmigrationRates     = Array: [I,S,T], Immigration rates of each country for each age cohort and year
            - self.N                    = Array: [I,S,T], UPDATED population of each country for each age cohort and year
            - self.Nhat                 = Array: [I,S,T+S], UPDATED world population share of each country for each age cohort and year
            - self.g_N                  = Array: [T], Population growth rate each year
            - self.Nhat_ss              = Array: [I,S], Population of each country for each age cohort in the steady state
            - self.Mortality_ss         = Array: [I,S], Mortality rates of each country for each age cohort in the steady state
            - self.MortalityRates       = Array: [I,S,T+S], UPDATED mortality rates of each country for each age cohort and year

        Other Functions Called:
            - None

        Objects in Function:
            - N_temp                    = Array: [I,S,T], Matrix created to help calculate the population dynamics TODO: CHANGE StepbyStep to have this way of using the equations
            - pop_old                   = Array: [I,S,T], Population shares in a given year beyond T
                                                          that is compared with pop_new to determine the steady state
            - pop_new                   = Array: [I,S,T], Population shares in a given year beyond T
                                                          that is compared with pop_old to determine the steady state
            - future_year_iter          = Int: Counter that keeps track of how many years beyond T it takes 
                                               for the population shares to converge to the steady state

        Outputs:
            - None

        """

        #Initializes immigration rates
        self.ImmigrationRates = np.zeros((self.I,self.S,self.T))

        #Initialize helper matrix in calculating population dynamics
        N_temp = np.ones((self.I,self.S))/(self.I*self.S)

        #Getting the population and population shares from the present to year T
        for t in xrange(1,self.T):

            #Gets new babies born this year (Equation 3.11)
            self.N[:,0,t] = np.sum((self.N[:,:,t-1]*self.FertilityRates[:,:,t-1]), axis=1)
            N_temp[:,0] = np.sum((self.Nhat[:,:,t-1]*self.FertilityRates[:,:,t-1]), axis=1)

            #Get the immigration RATES for the past year
            #If before the transition year T_1, just divide total migrants by population
            if t <= self.T_1:
                self.ImmigrationRates[:,:,t-1] = self.Migrants[:,:,t-1]/self.N[:,:,t-1]

            #If beyond the transition year T_1, average the immigration rates in year T_1 itself
            else:
                self.ImmigrationRates[:,:,t-1] = np.mean(self.ImmigrationRates[:,:,self.T_1-1],\
                        axis=0)

            #Gets the non-newborn population for the next year (Equation 3.12)
            self.N[:,1:,t] = self.N[:,:-1,t-1]*(1+self.ImmigrationRates[:,:-1,t-1]-self.MortalityRates[:,:-1,t-1])
            N_temp[:,1:] = self.Nhat[:,:-1,t-1]*(1+self.ImmigrationRates[:,:-1,t-1]-self.MortalityRates[:,:-1,t-1])
            
            #Gets the population share by taking a fraction of the total world population this year
            self.Nhat[:,:,t] = self.N[:,:,t]/np.sum(self.N[:,:,t])

            #Getting the growth rate
            self.g_N[t] = np.sum(N_temp)-1

        #Gets Immigration rates for the final year
        self.ImmigrationRates[:,:,t] = self.Migrants[:,:,t]/self.N[:,:,t]

        #Initialize iterating variables to find the steady state population shares
        pop_old = self.N[:,:,-1]
        pop_new = self.N[:,:,-1]
        future_year_iter = 0

        #Calculates new years of population shares until the greatest absolute difference between 2 consecutive years is less than demog_ss_tol
        while np.max(np.abs(self.Nhat[:,:,-1] - self.Nhat[:,:,-2])) > demog_ss_tol:
            pop_new[:,0] = np.sum((pop_old[:,:]*self.FertilityRates[:,:,-1]),axis=1)
            pop_new[:,1:] = pop_old[:,:-1]*(1+self.ImmigrationRates[:,:-1,-1]\
                    -self.MortalityRates[:,:-1,-1])
            self.Nhat = np.dstack((self.Nhat,pop_new/np.sum(pop_new)))
            future_year_iter += 1

        if self.PrintLoc: print "The SS Population Share converged in", future_year_iter, "years beyond year T"

        #Stores the steady state year in a seperate matrix
        self.Nhat_ss = self.Nhat[:,:,-1]
        self.Mortality_ss=self.MortalityRates[:,:,-1]
        
        #Deletes all the years between t=T and the steady state calculated in the while loop
        self.Nhat = self.Nhat[:,:,:self.T]

        #Imposing the ss for years after self.T
        self.Nhat = np.dstack((  self.Nhat[:,:,:self.T], np.einsum("is,t->ist",self.Nhat_ss,np.ones(self.S))  ))

        #Imposing the ss for years after self.T
        self.MortalityRates = np.dstack((  self.MortalityRates[:,:,:self.T], np.einsum("is,t->ist",self.Mortality_ss, np.ones(self.S))  ))        

        #Overwrites all the years in the transition path with the steady state if UseSSDemog == True
        if UseSSDemog == True:
            self.Nhat = np.einsum("is,t->ist",self.Nhat_ss,np.ones(self.T+self.S))
            self.MortalityRates = np.einsum("is,t->ist",self.Mortality_ss,np.ones(self.T+self.S))

    def plotDemographics(self, T_touse="default", compare_across="T", data_year=0):
        """
        Description: This just calls the plotDemographics function from the AuxiliaryDemographics.py file. See it for details
        """

        ages = self.FirstFertilityAge, self.LastFertilityAge, self.FirstDyingAge, self.MaxImmigrantAge
        datasets = self.FertilityRates, self.MortalityRates, self.ImmigrationRates, self.Nhat

        #Calls the Auxiliary Demographics file for this function
        demog.plotDemographics(ages, datasets, self.I, self.S, self.T, self.I_touse, T_touse, compare_across, data_year)

    #STEADY STATE

    def get_Psi(self, w, e):
        """
        Description:
            - Calculates the variable Psi using equation 3.21 for the steady state and for transition path functions

        Inputs:
            - w          = Array: [I,T] or [I], Wage rate for each country and year if called from transition path.
                                  Otherwise it is the steady state wage rate
            - e          = Array: [I,S,T] or [I,S], Labor Productivities for each country, cohort and year if called from transition path. 
                                  Otherwise it is the steady state labor productivities

        Variables Called from Object:
            - self.chi   = Scalar: Leisure Preference Parameter
            - self.rho   = Scalar: The intratemporal elasticity of substitution between consumption and leisure
            - self.sigma = Scalar: Rate of Time Preference

        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Objects in Function:
            - we         = Array: [I,S,T] or [I,S], Matrix product of w and e
            - psi        = Array: [I,S,T] or [I,S], Variable made just to simplify calculation of household decision equations

        Outputs:
            - psi
        """
        #If getting the SS
        if e.ndim == 2:
            we =  np.einsum("i,is->is",w,e)

        #If getting transition path
        elif e.ndim == 3:
            we = np.einsum("it, ist -> ist", w, e)

        psi = ( 1 + self.chi*( (self.chi/we)**(self.rho-1) ) )**( (1-self.rho*self.sigma)/(self.rho-1) )

        return psi

    def get_lhat(self,c,w,e):
        """
        Description:
            - Gets household leisure based on equation 3.20

        Inputs:
            - c             = Array: [I,S,T] or [I,S], Consumption for either the transition path or the steady steady-state
            - w             = Array: [I,T] or [I], Wage rate for either the transition path or the steady steady-state
            - e             = Array: [I,S,T] or [I,S], Labor productivities for either the transition path or the steady steady-state

        Variables Called from Object:
            - self.chi      = Scalar: Leisure preference parameter
            - self.rho      = Scalar: The intratemporal elasticity of substitution between consumption and leisure

        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Objects in Function:
            - None

        Outputs:
            - lhat          = Array: [I,S,T] or [I,S], Leisure for either the transition path or the steady steady-state
        """

        if e.ndim == 2:
            lhat=c*(self.chi/np.einsum("i,is->is",w,e))**self.rho
        elif e.ndim == 3:
            lhat=c*(self.chi/np.einsum("it,ist->ist",w,e))**self.rho

        return lhat

    def get_n(self, lhat):
        """
        Description:
            -Calculates the aggregate labor productivity based on equation (3.14)

        Inputs:
            - lhat          = Array: [I,S,T] or [I,S], Leisure for either the transition path or the steady steady-state

        Variables Called from Object:
            - self.e        = Array: [I,S,T], Labor productivities for the transition path  
            - self.e_ss     = Array: [I,S], Labor produtivities for the Steady State
            - self.lbar     = Array: [T+S], Time endowment in each year
            - self.Nhat     = Array: [I,S,T+S], World population share of each country for each age cohort and year
            - self.Nhat_ss  = Array: [I,S], Population of each country for each age cohort in the steady state
            - self.lbar_ss  = Int: Steady state time endowment. Normalized to 1.0
            - self.T        = Int: Number of Time Periods

        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Objects in Function:
            - None

        Outputs:
            - n          = Array: [I,S,T] or [I,S], Aggregate labor productivity for either the transition path or the steady steady-state
        """

        if lhat.ndim == 2:
            n = np.sum(self.e_ss*(self.lbar_ss-lhat)*self.Nhat_ss,axis=1)
        elif lhat.ndim == 3:
            n = np.sum(self.e[:,:,:self.T]*(self.lbar[:self.T]-lhat)*self.Nhat[:,:,:self.T],axis=1)

        return n

    def get_Y(self, kd, n):
        """
        Description:
            -Calculates the aggregate output based on equation (3.15)

        Inputs:
            - kd         = Array: [I,S,T] or [I,S], Domestic owned capital path for either the transition path or steady-state.
            - n          = Array: [I,S,T] or [I,S], Aggregate labor productivity for either the transition path or the steady steady-state

        Variables Called from Object:
            - self.A     = Array: [I,1], Technology level for each country
            - self.alpha = Scalar: Capital share of production

        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Objects in Function:
            - None

        Outputs:
            - Y          = Array: [I,S,T] or [I,S], Total output from firms for either the transition path or the steady steady-state
        """

        if kd.ndim ==1:
            Y = (kd**self.alpha) * ((self.A*n)**(1-self.alpha))
        elif kd.ndim== 2:
            Y = (kd**self.alpha) * (np.einsum("i,is->is",self.A,n)**(1-self.alpha))

        return Y

    def GetSSComponents(self, bq_ss, r_ss):
        """
        Description:
            - Solves for all the other variables in the model using bq_ss and r_ss

        Inputs:
            - bq_ss                     = Array: [I,S], 
            - r_ss                      = Scalar: Steady-state intrest rate

        Variables Called from Object:
            - self.alpha             = Scalar: Capital share of production
            - self.A                    = Array: [I], Technology level for each country
            - self.e_ss                 = Array: [I,S], Labor produtivities for the Steady State
            - self.I                    = Int: Number of Countries
            - self.Nhat_ss              = Array: [I,S,T+S], World population share of each country for each age cohort and year

        Variables Stored in Object:
            - None

        Other Functions Called:
            - householdEuler_SS = System of Euler equations to solve the household problem. Used by opt.fsolve
            - get_Psi = Solves for the Psi variable as in Equation 3.21
            - get_lhat = Solves for leisure as in Equation 3.20
            - get_n = Solves for labor supply as in Equation 3.14
            - get_Y = Solves for output as in Equation 3.15

        Objects in Function:
            - w_ss                      = Array: [I], Steady state wage rate
            - psi_ss                    = Array: [I,S], Steady state Psi variable (see equation 3.21)
            - c1_guess                  = Array: [I,S], Initial guess for consumption of the youngest cohort 
            - opt_c1                    = Array: [I,S], Optimal consumption of the youngest cohort 
            - cvec_ss                   = Array: [I,S], Steady state consumption for each country and cohort
            - avec_ss                   = Array: [I,S], Steady state assets holdings for each country and cohort
            - lhat_ss                   = Array: [I,S], Steady state leisure decision for each country and cohort.
            - n_ss                      = Array: [I], Steady state labor supply
            - kd_ss                     = Array: [I], Steady state total capital holdings for each country
            - y_ss                      = Array: [I], Steady state output of each country
            - kf_ss                     = Array: [I], Steady state foreign capital in each country

        Outputs:
            - w_ss, cvec_ss, avec_ss, kd_ss, kf_ss, n_ss, y_ss, and lhat_ss
        """

        def get_lifetime_decisionsSS(c_1, w_ss, r_ss, psi_ss):
            """
            Description:
                -Description of the Function

            Inputs:
                - c_1
                - w_ss
                - r_ss
                - psi_ss

            Variables Called from Object:
                - self.I                 = Int: Number of Countries
                - self.S                     = Int: Number of Cohorts
                - self.beta              = Scalar: Calculated overall future discount rate
                - self.Mortality_ss
                - self.delta             = Scalar: Calulated overall depreciation rate
                - self.sigma             = Scalar: Rate of Time Preference
                - self.g_A              = Scalar: Growth rate of technology
                - self.e_ss
                - self.chi       = Scalar: Leisure preference parameter

            Variables Stored in Object:
                - None

            Other Functions Called:
                - None

            Objects in Function:
                - None

            Outputs:
                - cvec_ss
                - avec_ss
            """


            cvec_ss = np.zeros((self.I,self.S))
            avec_ss = np.zeros((self.I,self.S+1))
            cvec_ss[:,0] = c_1

            for s in xrange(self.S-1):
                cvec_ss[:,s+1] = (self.beta * (1-self.Mortality_ss[:,s]) * (1 + r_ss - self.delta)\
                        *psi_ss[:,s+1]/psi_ss[:,s])**(1/self.sigma) * cvec_ss[:,s]*np.exp(-self.g_A)

                avec_ss[:,s+1] = (w_ss*self.e_ss[:,s] + (1 + r_ss - self.delta)*avec_ss[:,s] + \
                        bq_ss[:,s] - cvec_ss[:,s]*(1+w_ss*self.e_ss[:,s]*\
                        (self.chi/(w_ss*self.e_ss[:,s]))**self.rho))*np.exp(-self.g_A)

            avec_ss[:,s+2] = (w_ss*self.e_ss[:,s+1] + (1 + r_ss - self.delta)*avec_ss[:,s+1] \
                    - cvec_ss[:,s+1]*(1+w_ss*self.e_ss[:,s+1]*(self.chi/(w_ss*self.e_ss[:,s+1]))\
                    **self.rho))*np.exp(-self.g_A)

            return cvec_ss, avec_ss

        def householdEuler_SS(c_1, w_ss, r_ss, psi_ss):
            """
            Description:
                -Description of the Function

            Inputs:
                - c_1
                - w_ss
                - r_ss
                - psi_ss

            Variables Called from Object:
                - None

            Variables Stored in Object:
                - None

            Other Functions Called:
                - get_lifetimedecisionsSS

            Objects in Function:
                - None

            Outputs:
                - Euler
            """


            cpath, assets_path = get_lifetime_decisionsSS(c_1, w_ss, r_ss, psi_ss)

            Euler = np.ravel(assets_path[:,-1])

            if np.any(cpath<0):
                print "WARNING! The fsolve for initial optimal consumption guessed a negative number"
                Euler = np.ones(Euler.shape[0])*9999.

            return Euler

        #Equation 3.25
        w_ss = (self.alpha*self.A/r_ss)**(self.alpha/(1-self.alpha))*(1-self.alpha)*self.A

        #Equation 3.21
        psi_ss = self.get_Psi(w_ss,self.e_ss)

        #Initial guess for the first cohort's consumption
        c1_guess = np.ones(self.I)*.02

        #Finds the optimal consumption for the first cohort
        opt_c1 = opt.fsolve(householdEuler_SS, c1_guess, args = (w_ss, r_ss, psi_ss))

        #Gets the optimal paths for consumption and assets as a function of the first cohort's consumption
        cvec_ss, avec_ss = get_lifetime_decisionsSS(opt_c1,w_ss,r_ss,psi_ss)

        #Snips off the final entry of assets since it is just 0 if the equations solved correctly
        avec_ss = avec_ss[:,:-1]

        #Equation 3.20
        lhat_ss = self.get_lhat(cvec_ss, w_ss, self.e_ss)

        #Equation 3.14
        n_ss = self.get_n(lhat_ss)

        #Equation 3.26
        kd_ss = np.sum(avec_ss*self.Nhat_ss,axis=1)

        #Equation 3.15
        y_ss = self.get_Y(kd_ss,n_ss)

        #Equation 3.27
        kf_ss = (self.alpha*self.A/r_ss)**(1/(1-self.alpha)) * n_ss-kd_ss

        return w_ss, cvec_ss, avec_ss, kd_ss, kf_ss, n_ss, y_ss, lhat_ss

    def EulerSystemSS(self, guess, PrintSSEulErrors=False):
        """
        Description:
            - System of Euler equations that must be satisfied (or = 0) for the ss to solve. 

        Inputs:
            - guess                     = Array: [I+1], Contains guesses for individual bequests in each country 
                                                        and the guess for the world intrest rate
             - PrintSSEulErrors          = Boolean, True prints the Euler Errors in each iteration of calculating the steady state

        Variables Called from Object:
            - self.I                    = Int: Number of Countries
            - self.S                    = Int: Number of Cohorts
            - self.FirstFertilityAge    = Int: First age where agents give birth
            - self.FirstDyingAge        = Int: First age where mortality rates effect agents
            - self.Mortality_ss         = Array: [I,S], Mortality rates of each country for each age cohort in the steady state
            - self.Nhat_ss              = Array: [I,S,T+S], World population share of each country for each age cohort and year

        Variables Stored in Object:
            - None

        Other Functions Called:
            - GetSSComponents = System of equations that solves for wages, consumption, assets, 
                                capital stocks, labor input, domestic output, and leisure in terms 
                                of the world intrest rate and bequests

        Objects in Function:
            - bqindiv_ss                = Array: [I], Current guess for the amount of bequests each eligible-aged 
                                                      individual will receive in each country
            - r_ss                      = Scalar: Current guess for the steady-state intrest rate
            - bq_ss                     = Array: [I,S], Vector of bequests received for each cohort and country.
                                                        Basically bqindiv_ss copied for each eligible-aged individual.
            - w_ss                      = Array: [I], Current guess for each countries ss wage rate as a function of r_ss and bqvec_ss
            - cvec_ss                   = Array: [I,S], Current guess for ss consumption for each country and cohort
            - avec_ss                   = Array: [I,S], Current guess for the ss assets holdings for each country and cohort
            - kd_ss                     = Array: [I], Current guess for ss total capital holdings for each country
            - kf_ss                     = Array: [I], Current guess for ss foreign capital in each country
            - n_ss                      = Array: [I], Current guess for ss labor supply
            - y_ss                      = Array: [I], Current guess for ss output of each country
            - lhat_ss                   = Array: [I,S], Current guess for ss leisure decision for each country and cohort.
            - alldeadagent_assets       = Array: [I], Sum of assets of all the individuals who die in the steady state. 
                                                      Evenly distributed to eligible-aged cohorts.
            - Euler_bq                  = Array: [I], Distance between bqindiv_ss and the actual bqindiv_ss calculated in the system. 
                                                      Must = 0 for the ss to correctly solve.
            - Euler_kf                  = Scalar: Sum of the foreign capital stocks. Must = 0 for the ss to correctly solve

        Outputs:
            - Euler_all                 = Array: [I+1], Euler_bq and Euler_kf stacked together. Must = 0 for the ss to correctly solve
        """
        #Breaking up the input into its 2 components
        bqindiv_ss = guess[:-1]
        r_ss = guess[-1]

        #Initializes a vector of bequests received for each individial. Will be = 0 for a block of young and a block of old cohorts
        bq_ss = np.zeros((self.I,self.S))
        bq_ss[:,self.FirstFertilityAge:self.FirstDyingAge] = \
                np.einsum("i,s->is", bqindiv_ss, np.ones(self.FirstDyingAge-self.FirstFertilityAge))

        #Calls self.GetSSComponents, which solves for all the other ss variables in terms of bequests and intrest rate
        w_ss, cvec_ss, avec_ss, kd_ss, kf_ss, n_ss, y_ss, lhat_ss = self.GetSSComponents(bq_ss, r_ss)

        #Sum 
        alldeadagent_assets = np.sum(avec_ss[:,self.FirstDyingAge:]*\
                self.Mortality_ss[:,self.FirstDyingAge:]*self.Nhat_ss[:,self.FirstDyingAge:], axis=1)

        Euler_bq = bqindiv_ss - alldeadagent_assets/np.sum(self.Nhat_ss[:,self.FirstFertilityAge:self.FirstDyingAge],\
                axis=1)

        Euler_kf = np.sum(kf_ss)

        Euler_all = np.append(Euler_bq, Euler_kf)

        if PrintSSEulErrors: print "Euler Errors:", Euler_all
        
        return Euler_all

    def SteadyState(self, rss_guess, bqss_guess, PrintSSEulErrors=False):
        """
        Description:
            -Description of the Function

        Inputs:
            -

        Variables Called from Object:
            -

        Variables Stored in Object:
            -

        Other Functions Called:
            -

        Objects in Function:
            -

        Outputs:
            -
        """

        guess = np.append(bqss_guess, rss_guess)

        ss = opt.fsolve(self.EulerSystemSS, guess, args=PrintSSEulErrors)

        self.bq_ss = ss[:-1]
        self.r_ss = ss[-1]

        self.bqvec_ss = np.zeros((self.I,self.S))
        self.bqvec_ss[:,self.FirstFertilityAge:self.FirstDyingAge] = np.einsum("i,s->is",self.bq_ss,\
                np.ones(self.FirstDyingAge-self.FirstFertilityAge))

        self.w_ss, self.cvec_ss, self.avec_ss, self.kd_ss, self.kf_ss, self.n_ss, self.y_ss, self.lhat_ss \
                = self.GetSSComponents(self.bqvec_ss,self.r_ss)


        alldeadagent_assets = np.sum(self.avec_ss[:,self.FirstDyingAge:]*self.Mortality_ss[:,self.FirstDyingAge:]*\
                self.Nhat_ss[:,self.FirstDyingAge:], axis=1)

        Euler_bq = self.bq_ss - alldeadagent_assets/np.sum(self.Nhat_ss[:,self.FirstFertilityAge:self.FirstDyingAge],\
                axis=1)
        Euler_kf = np.sum(self.kf_ss)

        print "\n\nSTEADY STATE FOUND!"
        if self.CheckerMode==False:
            print "-Euler for bq satisfied:", np.isclose(np.max(np.absolute(Euler_bq)), 0)
            print "-Euler for r satisfied:", np.isclose(Euler_kf, 0), "\n\n"

    def checkSSEulers(self):
        """
        Description:
            -Description of the Function

        Inputs:
            - None

        Variables Called from Object:
            - self.w_ss
            - self.e_ss
            - self.psi_ss
            - self.cvec_ss
            - self.sigma             = Scalar: Rate of Time Preference
            - self.beta              = Scalar: Calculated overall future discount rate
            - self.Mortality_ss
            - self.psi_ss
            - self.g_A              = Scalar: Growth rate of technology
            - self.r_ss
            - self.delta             = Scalar: Calulated overall depreciation rate
            - self.avec_ss
            - self.bqvec_ss

        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Objects in Function:
            - we

        Outputs:
            - None
        """

        we = np.einsum("i,is->is",self.w_ss,self.e_ss[:,:-1])

        print self.psi_ss[:,:-1]*self.cvec_ss[:,:-1]**(-self.sigma) - self.beta*(1-self.Mortality_ss[:,:-1])*self.psi_ss[:,1:]*(self.cvec_ss[:,1:]*np.exp(self.g_A))**(-self.sigma)*(1+self.r_ss-self.delta)
        
        print self.cvec_ss[:,:-1] - \
        (we + (1+self.r_ss-self.delta)*self.avec_ss[:,:-1] + self.bqvec_ss[:,:-1] - self.avec_ss[:,1:]*np.exp(self.g_A)) / \
        (1 + we*(self.chi/we)**self.rho)
  
    def PrintSSResults(self):
        """
        Description:
            -Prints the final result of steady state calculations

        Inputs:
            - None

        Variables Called from Object:
            - self.avec_ss
            - self.kf_ss
            - self.kd_ss
            - self.n_ss
            - self.y_ss
            - self.r_ss
            - self.w_ss
            - self.cvec_ss

        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Objects in Function:
            - None

        Outputs:
            - None
        """
        print "assets steady state:", self.avec_ss
        print "kf steady state", self.kf_ss
        print "kd steady state", self.kd_ss
        print "n steady state", self.n_ss
        print "y steady state", self.y_ss
        print "r steady state", self.r_ss
        print "w steady state", self.w_ss
        print "c_vec_ss steady state", self.cvec_ss

    def plotSSResults(self):
        """
        Description:
            - Plots the final calculations of the Steady State

        Inputs:
            - None

        Variables Called from Object:
            - self.I                 = Int: Number of Countries
            - self.S                     = Int: Number of Cohorts
            - self.cvec_ss
            - self.avec_ss
            - self.bqvec_ss

        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Objects in Function:
            - None

        Outputs:
            - None
        """
        for i in range(self.I):
            plt.plot(range(self.S),self.cvec_ss[i,:])
        plt.title("Consumption")
        plt.legend(self.I_touse[:self.I])
        plt.show()
        for i in range(self.I):
            plt.plot(range(self.S),self.avec_ss[i,:])
        plt.title("Assets")
        plt.legend(self.I_touse[:self.I])
        plt.show()
        for i in range(self.I):
            plt.plot(range(self.S),self.bqvec_ss[i,:])
        plt.title("Bequests")
        plt.legend(self.I_touse[:self.I])
        plt.show()

    #TIMEPATH-ITERATION

    def set_initial_values(self, r_init, bq_init, a_init):
        """
        Description:
            - Saves the initial guesses of r, bq and a given by the user into the object

        Inputs:
            - r_init        = Scalar: Initial interest rate given by User
            - bq_init       = Array: [I], Initial bequests given by User
            - a_init        = Array: [I,S], Initial asset distribution given by User

        Variables Called from Object:
            - None


        Variables Stored in Object:
            - self.r_init   = Scalar: Initial interest rate given by User
            - self.bq_init  = Array: [I], Initial bequests given by User
            - self.a_init   = Array: [I,S], Initial asset distribution given by Users

        Other Functions Called:
            - None

        Objects in Function:
            - None

        Outputs:
            - None

        """

        self.r_init = r_init
        self.bq_init = bq_init
        self.a_init = a_init

    def get_initialguesses(self):
        """
        Description:
            - Generates an initial guess path used for beginning TPI calculation. It follows the form
              of a quadratic function:

              y = aa x^2 + bb x + cc

        Inputs:
            - None

        Variables Called from Object:
            - self.bq_init  = Array: [I], Initial bequests given by User
            - self.I        = Int: Number of Countries
            - self.T        = Int: Number of Time Periods
            - self.r_init   = Scalar: Initial interest rate given by User
            - self.r_ss     = Scalar: Steady state interest rate
      
        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Objects in Function:
            - aa            = coefficient for x^2 term
            - bb            = coefficient for x term
            - cc            = coefficient for constant term

        Outputs:
            - rpath_guess   = Array: [T], Initial path of interest rates in quadratic form
            - bqpath_guess  = Array: [I,T], Initial path of bequests in quadratic form

        """


        rpath_guess = np.zeros(self.T)
        bqpath_guess = np.zeros((self.I,self.T))

        cc = self.r_init
        bb = -2 * (self.r_init-self.r_ss)/(self.T-1)
        aa = -bb / (2*(self.T-1))
        rpath_guess[:self.T] = aa * np.arange(0,self.T)**2 + bb*np.arange(0,self.T) + cc

        for i in range(self.I):
            bqpath_guess[i,:self.T] = np.linspace(self.bq_init[i], self.bq_ss[i], self.T)

        return rpath_guess, bqpath_guess

    def GetTPIComponents(self, bqvec_path, r_path, Print_HH_Eulers, Print_cabqTimepaths):
        """
        Description:
            -Description of the Function

        Inputs:
            - bqvec_path
            - r_path
            - Print_HH_Eulers
            - Print_cabqTimepaths

        Variables Called from Object:
            - None

        Variables Stored in Object:
            - None

        Other Functions Called:
            - self.get_Psi
            - get_c_a_matrices
            - self.get_lhat
            - self.get_n
            - self.get_Y

        Objects in Function:
            - psi

        Outputs:
            - w_path
            - c_matrix
            - a_matrix
            - kd_path
            - kf_path
            - n_path
            - y_path
            - lhat_path
        """

        #Functions that solve lower-diagonal household decisions in vectors
        def get_lifetime_decisions_LOWERTRIANGLETEST(c0_guess, c_uppermat, a_uppermat, w_path, r_path, psi, bqvec_path):
            """
            Description:
                -Description of the Function

            Inputs:
                -c0_guess
                -c_uppermat
                -a_uppermat
                -w_path
                -r_path
                -psi
                -bqvec_path

            Variables Called from Object:
                - self.I                 = Int: Number of Countries
                - self.T                 = Int: Number of time periods
                - self.e
                - self.S                     = Int: Number of Cohorts
                - self.beta              = Scalar: Calculated overall future discount rate
                - self.MortalityRates
                - self.delta             = Scalar: Calulated overall depreciation rate
                - self.sigma             = Scalar: Rate of Time Preference
                - self.g_A              = Scalar: Growth rate of technology
                - self.chi       = Scalar: Leisure preference parameter
                - self.rho       = Scalar: The intratemporal elasticity of substitution between consumption and leisure

            Variables Stored in Object:
                - None

            Other Functions Called:
                - None

            Objects in Function:
                - None

            Outputs:
                - c_matrix
                - a_matrix
            """


            #Initializes consumption and assets with all of the upper triangle already filled in
            c_matrix = c_uppermat
            a_matrix = a_uppermat
            c_matrix[:,0,:self.T] = c0_guess.reshape(self.I,self.T)

            #Gets we ahead of time for easier calculation
            we = np.einsum("it,ist->ist",w_path,self.e)

            #Loops through each year (across S) and gets decisions for every agent in the next year
            for s in range(self.S-1):

                #Gets consumption for every agents' next year using Equation 3.22
                c_matrix[:,s+1,s+1:self.T+s+1] = ((self.beta * (1-self.MortalityRates[:,s,s:self.T+s]) * (1 + r_path[s+1:self.T+s+1] - self.delta)\
                                                 * psi[:,s+1,s+1:self.T+s+1])/psi[:,s,s:self.T+s])**(1/self.sigma) * c_matrix[:,s,s:self.T+s]*np.exp(-self.g_A)
                #Gets assets for every agents' next year using Equation 3.19
                a_matrix[:,s+1,s+1:self.T+s+1] = (  (we[:,s,s:self.T+s] + (1 + r_path[s:self.T+s] - self.delta)*a_matrix[:,s,s:self.T+s] + bqvec_path[:,s,s:self.T+s])\
                                                 -c_matrix[:,s,s:self.T+s]*(1+we[:,s,s:self.T+s]*(self.chi/we[:,s,s:self.T+s])**self.rho)  )*np.exp(-self.g_A)

            #Gets assets in the final period of every agents' lifetime
            a_matrix[:,-1,s+2:self.T+s+2] = (  (we[:,-1,s+1:self.T+s+1] + (1 + r_path[s+1:self.T+s+1] - self.delta)*a_matrix[:,-2,s+1:self.T+s+1])\
                                            -c_matrix[:,-1,s+1:self.T+s+1]*(1+we[:,-1,s+1:self.T+s+1]*(self.chi/we[:,-1,s+1:self.T+s+1])**self.rho)  )*np.exp(-self.g_A)


            return c_matrix, a_matrix

        def get_lower_triangle_Euler_TEST(c0_guess, c_uppermat, a_uppermat, w_path, r_path, psi, bqvec_path):
            """
            Description:
                -Description of the Function

            Inputs:
                - c0_guess
                - c_uppermat
                - a_uppermat
                - w_path
                - r_path
                - psi
                - bqvec_path

            Variables Called from Object:
                - None

            Variables Stored in Object:
                - None

            Other Functions Called:
                - get_lifetime_decisions_LOWERTRIANGLETEST

            Objects in Function:
                - None 

            Outputs:
                - Euler

            """

            
            #Gets the decisions paths for each agent
            c_matrix, a_matrix = get_lifetime_decisions_LOWERTRIANGLETEST(c0_guess, c_uppermat, a_uppermat, w_path, r_path, psi, bqvec_path)
            
            #Household Eulers are solved when the agents have no assets at the end of their life
            Euler = np.ravel(a_matrix[:,-1,self.S:])

            #print np.round(a_matrix[0,-1,self.S:], decimals=3)

            return Euler

        #Functions that solve upper-diagonal household decisions in vectors
        def get_lifetime_decisions_UPPERTRIANGLETEST(c0_guess, c_matrix, a_matrix, w_path, r_path, psi, bqvec_path):
            """
            Description:
                -Description of the Function

            Inputs:
                - c0_guess
                - c_matrix
                - a_matrix
                - w_path
                - r_path
                - psi
                - bqvec_path

            Variables Called from Object:
                - self.beta              = Scalar: Calculated overall future discount rate
                - self.MortalityRates
                - self.delta             = Scalar: Calulated overall depreciation rate
                - self.sigma             = Scalar: Rate of Time Preference
                - self.g_A              = Scalar: Growth rate of technology
                - self.chi       = Scalar: Leisure preference parameter
                - self.rho       = Scalar: The intratemporal elasticity of substitution between consumption and leisure

            Variables Stored in Object:
                - None

            Other Functions Called:
                - None

            Objects in Function:
                - we

            Outputs:
                - c_matrix
                - a_matrix
            """

            
            c_matrix[:,:-1,0] = c0_guess.reshape(self.I,self.S-1)
            we = np.einsum("it,ist->ist",w_path,self.e)
            #print np.round(np.transpose(c_matrix[0,:,:self.S+3]), decimals=3)

            for s in range(self.S):
                t = s
                c_matrix[:,s+1:,t+1] = ((self.beta * (1-self.MortalityRates[:,s:-1,t]) * (1 + r_path[t+1] - self.delta)\
                                                 * psi[:,s+1:,t+1])/psi[:,s:-1,t])**(1/self.sigma) * c_matrix[:,s:-1,t]*np.exp(-self.g_A)
                #print np.round(np.transpose(c_matrix[0,:,:self.S+3]), decimals=3)
                
                a_matrix[:,s+1:,t+1] = (  (we[:,s:,t] + (1 + r_path[t] - self.delta)*a_matrix[:,s:-1,t] + bqvec_path[:,s:,t])\
                                                 -c_matrix[:,s:,t]*(1+we[:,s:,t]*(self.chi/we[:,s:,t])**self.rho)  )*np.exp(-self.g_A)

            #Gets assets in the final period of every agents' lifetime
            a_matrix[:,-1,t+2] = (  (we[:,-1,t+1] + (1 + r_path[t+1] - self.delta)*a_matrix[:,-2,t+1])\
                                            -c_matrix[:,-1,t+1]*(1+we[:,-1,t+1]*(self.chi/we[:,-1,t+1])**self.rho)  )*np.exp(-self.g_A)

            return c_matrix, a_matrix

        def get_upper_triangle_Euler_TEST(c0_guess, c_matrix, a_matrix, w_path, r_path, psi, bqvec_path):
            """
            Description:
                -Description of the Function

            Inputs:
                - c0_guess
                - c_matrix
                - a_matrix
                - w_path
                - r_path
                - psi
                - bqvec_path

            Variables Called from Object:
                - self.S                     = Int: Number of Cohorts

            Variables Stored in Object:
                - None

            Other Functions Called:
                - get_lifetime_decisions_UPPERTRIANGLETEST

            Objects in Function:
                - None

            Outputs:
                - Euler
            """

            #Gets the decisions paths for each agent
            c_matrix, a_matrix = get_lifetime_decisions_UPPERTRIANGLETEST(c0_guess, c_matrix, a_matrix, w_path, r_path, psi, bqvec_path)
            
            #Household Eulers are solved when the agents have no assets at the end of their life
            Euler = np.ravel(a_matrix[:,-1,1:self.S])

            #print np.round(a_matrix[0,-1,self.S:], decimals=3)

            return Euler

        #Functions that solve household decisions with for-loops (old)
        def get_lifetime_decisionsTPI(c_1, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age):

            """
            Description:
                -Description of the Function

            Inputs:
                - c_1
                - w_life
                - r_life
                - mort_life
                - e_life
                - psi_life
                - bq_life
                - a_current
                - age

            Variables Called from Object:
                - self.beta              = Scalar: Calculated overall future discount rate
                - self.delta             = Scalar: Calulated overall depreciation rate
                - self.g_A              = Scalar: Growth rate of technology
                - self.rho       = Scalar: The intratemporal elasticity of substitution between consumption and leisure
                - self.chi       = Scalar: Leisure preference parameter

            Variables Stored in Object:
                - None

            Other Functions Called:
                - None

            Objects in Function:
                - decisions

            Outputs:
                - cvec_path
                - avec_path

            """

            #Number of decisions the agent needs to make in its lifetime
            decisions = self.S - age -1

            #Initializes the consumption and assets vectors for this agent
            cvec_path = np.zeros((self.I,decisions+1))
            avec_path = np.zeros((self.I,decisions+2))
            cvec_path[:,0] = c_1
            avec_path[:,0] = a_current

            #Loops through each decision the agent makes
            for s in range(decisions):

                #Algebraically manipulated version of Equation 3.22
                cvec_path[:,s+1] = ((self.beta * (1-mort_life[:,s]) * (1 + r_path[s+1] - self.delta)\
                                   * psi_life[:,s+1])/psi_life[:,s])**(1/self.sigma) * cvec_path[:,s]*np.exp(-self.g_A)

                #Algebraically manipulataed version of Equation 3.19
                avec_path[:,s+1] = (w_life[:,s]*e_life[:,s] + (1 + r_life[s] - self.delta)*avec_path[:,s] + \
                        bq_life[:,s] - cvec_path[:,s]*(1+w_life[:,s]*e_life[:,s]*\
                        (self.chi/(w_life[:,s]*e_life[:,s]))**self.rho))*np.exp(-self.g_A)

            #Gets the remaining assets in the final year of the agent't lifetime
            avec_path[:,s+2] = (w_life[:,s+1]*e_life[:,s+1] + (1 + r_life[s+1] - self.delta)*avec_path[:,s+1] \
                    - cvec_path[:,s+1]*(1+w_life[:,s+1]*e_life[:,s+1]*(self.chi/(w_path[:,s+1]*e_life[:,s+1]))\
                    **self.rho))*np.exp(-self.g_A)


            return cvec_path, avec_path

        def optc1_Euler_TPI(c1_guess, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age):

            """
            Description:
                -Description of the Function

            Inputs:
                - c1_guess
                - w_life
                - r_life
                - mort_life
                - e_life
                - psi_life
                - bq_life
                - a_current
                - age

            Variables Called from Object:
                - None

            Variables Stored in Object:
                - None

            Other Functions Called:
                - get_lifetime_decisionsTPI

            Objects in Function:
                - cpath_indiv
                - apath_indiv

            Outputs:
                - Euler

            """


            #Gets the individual decisions paths of the agent to check if finals assets are 0
            cpath_indiv, apath_indiv = get_lifetime_decisionsTPI(c1_guess, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age)
            
            #Household Eulers are solved when the agent doesn't have any assets at the end of its life
            Euler = np.ravel(apath_indiv[:,-1])

            if np.any(cpath_indiv<0):
                print "WARNING! The fsolve for initial optimal consumption guessed a negative number"
                Euler = np.ones(Euler.shape[0])*9999.

            return Euler

        #Checks various household condidions
        def check_household_conditions(w_path, r_path, c_matrix, a_matrix, psi, bqvec_path):  

            """
            Description:
                -Description of the Function

            Inputs:
                - w_path
                - r_path
                - c_matrix
                - a_matrix
                - psi
                - bqvec_path

            Variables Called from Object:
                - self.T                 = Int: Number of time periods
                - self.e
                - self.sigma             = Scalar: Rate of Time Preference
                - self.beta              = Scalar: Calculated overall future discount rate
                - self.g_A              = Scalar: Growth rate of technology
                - self.delta             = Scalar: Calulated overall depreciation rate
                - self.chi       = Scalar: Leisure preference parameter
                - self.rho       = Scalar: The intratemporal elasticity of substitution between consumption and leisure

            Variables Stored in Object:
                - None

            Other Functions Called:
                - None

            Objects in Function:
                - we

            Outputs:
                - Chained_C_Condition
                - Modified_Budget_Constraint
                - Modified_Budget_Constraint2
                - Household_Euler

            """

            
            #Multiplies wages and productivities ahead of time for easy calculations of the first two equations below
            we = np.einsum("it,ist->ist",w_path[:,:self.T-1],self.e[:,:-1,:self.T-1])

            #Disparity between left and right sides of Equation 3.22
            Chained_C_Condition = psi[:,:-1,:self.T-1]*c_matrix[:,:-1,:self.T-1]**(-self.sigma)\
                                  - self.beta*(1-self.MortalityRates[:,:-1,:self.T-1])*psi[:,1:,1:self.T]\
                                  *(c_matrix[:,1:,1:self.T]*np.exp(self.g_A))**(-self.sigma)*(1+r_path[1:self.T]-self.delta)
            
            #Disparity between left and right sides of Equation 3.19
            Modified_Budget_Constraint = c_matrix[:,:-1,:self.T-1]\
                                         -  (we + (1+r_path[:self.T-1]-self.delta)*a_matrix[:,:-2,:self.T-1] + bqvec_path[:,:-1,:self.T-1]\
                                         - a_matrix[:,1:-1,1:self.T]*np.exp(self.g_A))\
                                         /(1 + we*(self.chi/we)**self.rho)
            
            #Multiplies wages and productivities ahead of time for easy calculations of the two equations below
            we = np.einsum("it,ist->ist",w_path[:,:self.T],self.e[:,:-1,:self.T])
            
            #Disparity between left and right sides of Equation 3.19 algebraically manipulated so assets next period is the left-side variable
            Modified_Budget_Constraint2 = - a_matrix[:,1:-1,1:self.T+1]\
                                          + ( we + (1+r_path[:self.T] - self.delta)*a_matrix[:,:-2,:self.T] + bqvec_path[:,:-1,:self.T]\
                                          - (c_matrix[:,:-1,:self.T]*(1+we*(self.chi/we)**self.rho)) )*np.exp(-self.g_A)

            #Any remaining assets each agent has at the end of its lifetime. Should be 0 if other Eulers are solving correctly
            Household_Euler = a_matrix[:,-1,:]

            return Chained_C_Condition, Modified_Budget_Constraint, Modified_Budget_Constraint2, Household_Euler

        #Gets consumption and assets matrices using fsolve
        def get_c_a_matrices(w_path, r_path, psi, bqvec_path, Print_HH_Eulers, Print_cabqTimepaths):
            """
            Description:
                -Description of the Function

            Inputs:
                - w_path
                - r_path
                - psi
                - bqvec_path
                - Print_HH_Eulers
                - Print_cabqTimepaths

            Variables Called from Object:
                - self.I                 = Int: Number of Countries
                - self.S                     = Int: Number of Cohorts
                - self.T                 = Int: Number of time periods
                - self.e
                - self.delta             = Scalar: Calulated overall depreciation rate
                - self.a_init
                - self.chi       = Scalar: Leisure preference parameter
                - self.rho       = Scalar: The intratemporal elasticity of substitution between consumption and leisure
                - self.VectorizeHouseholdSolver
                - self.cvec_ss
                - self.MortalityRates


            Variables Stored in Object:
                - self.c0_alive
                - self.c0_future

            Other Functions Called:
                - get_upper_traingle_Euler_TEST
                - get_lower_triangle_Euler_TEST
                - optc1_Euler_TPI
                - get_lifetime_decisionsTPI

            Objects in Function:
                - c0alive_guess
                - p
                - c1_guess
                - w_life
                - r_life
                - mort_life
                - e_life
                - psi_life
                - bq_life
                - a_current
                - cpath_indiv
                - apath_indiv


            Outputs:
                - c_matrix[:,:,:self.T]
                - a_matrix[:,:-1,:self.T]

            """

            #Initializes the consumption and assets matrices
            c_matrix = np.zeros((self.I,self.S,self.T+self.S))
            a_matrix = np.zeros((self.I,self.S+1,self.T+self.S))
            a_matrix[:,:-1,0] = self.a_init

            #Equation 3.19 for the oldest agent in time t=0. Note that this agent chooses to consume everything so that it has no assets in the following period
            c_matrix[:,self.S-1,0] = (w_path[:,0]*self.e[:,self.S-1,0] + (1 + r_path[0] - self.delta)*self.a_init[:,self.S-1] + bqvec_path[:,self.S-1,0])\
            /(1+w_path[:,0]*self.e[:,self.S-1,0]*(self.chi/(w_path[:,0]*self.e[:,self.S-1,0]))**(self.rho))

            #Will solve the household matrices using vectorization if = True and by agent if = False
            if self.VectorizeHouseholdSolver:
            
                c0alive_guess = np.ones((self.I, self.S-1))*.3

                opt.fsolve(get_upper_triangle_Euler_TEST, c0alive_guess, args=(c_matrix, a_matrix, w_path, r_path, psi, bqvec_path))

                #Initializes a guess for the first vector for the fsolve to use

                c0future_guess = np.zeros((self.I,self.T))
                for i in range(self.I):
                    c0future_guess[i,:] = np.linspace(c_matrix[i,1,0], self.cvec_ss[i,-1], self.T)

                #Solves for the entire consumption and assets matrices
                opt.fsolve(get_lower_triangle_Euler_TEST, c0future_guess, args=(c_matrix, a_matrix, w_path, r_path, psi, bqvec_path))

                self.c0_alive = c_matrix[:,:-1,0]
                self.c0_future = c_matrix[:,0,:self.T]

            else:
            #Loops over each agent's lifetime decisions who is alive today (Upper triangle)
                for age in range(self.S-2,0,-1):

                    p = self.S-age #Remaining decisions

                    #Makes a guess for the fsolve for this agent's consumption that is a function of the consumption of the agent one year older
                    c1_guess = (c_matrix[:,age+1,0]*(psi[:,age,0]/psi[:,age+1,1])\
                        /((self.beta*(1+r_path[0]-self.delta))**(1/self.sigma)))/np.exp(self.g_A)

                    #All the variables this agent will face during its lifetime. (Note these are diagonal vectors)
                    w_life = w_path[:,:p]
                    r_life = r_path[:p+1]
                    mort_life = np.diagonal(self.MortalityRates[:,age:,age:], axis1=1, axis2=2)
                    e_life = np.diagonal(self.e[:,age:,:p+1], axis1=1, axis2=2)
                    psi_life = np.diagonal(psi[:,age:,:p+2], axis1=1, axis2=2)
                    bq_life = np.diagonal(bqvec_path[:,age:,:p+1], axis1=1, axis2=2)
                    a_current = self.a_init[:,age]

                    #Solves for this agent's optimal initial consumption in time t=0 using an fsolve
                    opt_c1 = opt.fsolve(optc1_Euler_TPI, c1_guess, args = (w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age))

                    #Solves for all the remaining lifetime decisions for this agent as a function of its optimal initial consumption using Equations 3.19 and 3.22
                    cpath_indiv, apath_indiv = get_lifetime_decisionsTPI(opt_c1, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age)
                    
                    #Fills the agents consumption and assets decision vectors as diagonals in the main matrix
                    for i in xrange(self.I):
                        np.fill_diagonal(c_matrix[i,age:,:], cpath_indiv[i,:])
                        np.fill_diagonal(a_matrix[i,age:,:], apath_indiv[i,:])

                    #Prints Consumption and Assets matrices if Print_cabqTimepaths = True in Main.py
                    if Print_cabqTimepaths:
                        print "Consumption for generation of age", age
                        print np.round(np.transpose(c_matrix[0,:,:self.T]), decimals=3)
                        print "Assets for generation of age", age
                        print np.round(np.transpose(a_matrix[0,:,:self.T]), decimals=3)

                #Loops through each agent yet to be born and gets that agents lifetime decisions from age 0 to death (Upper Triangle)
                for t in range(self.T):

                    age = 0

                    #Guess for initial consumption is based on agent one year older
                    c1_guess = c_matrix[:,0,t-1]

                    #Variables this agent will face over its lifetime
                    w_life = w_path[:,t:t+self.S]
                    r_life = r_path[t:t+self.S+1]
                    mort_life = np.diagonal(self.MortalityRates[:,:,t:t+self.S+1], axis1=1, axis2=2)
                    e_life = np.diagonal(self.e[:,:,t:t+self.S+1], axis1=1, axis2=2)
                    psi_life = np.diagonal(psi[:,:,t:t+self.S+1], axis1=1, axis2=2)
                    bq_life = np.diagonal(bqvec_path[:,:,t:t+self.S+1], axis1=1, axis2=2)
                    #Agents are born with no assets
                    a_current = np.zeros(self.I)

                    #Gets optimal initial consumption for this agent using an fsolve. Solved using equations 3.19 and 3.22
                    opt_c1 = opt.fsolve(optc1_Euler_TPI, c1_guess, args = (w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age))

                    #Gets optimal decisions paths for this agent
                    cpath_indiv, apath_indiv = get_lifetime_decisionsTPI(opt_c1, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age)

                    #Fills the agents consumption and assets decision vectors as diagonals in the main matrix
                    for i in xrange(self.I):
                        np.fill_diagonal(c_matrix[i,:,t:], cpath_indiv[i,:])
                        np.fill_diagonal(a_matrix[i,:,t:], apath_indiv[i,:])

                    #Prints Consumption and Assets matrices if Print_cabqTimepaths = True in Main.py                    
                    if Print_cabqTimepaths:
                        print "Consumption for year", t
                        print np.round(np.transpose(c_matrix[0,:,:self.T]), decimals=3)
                        print "Assets for year", t
                        print np.round(np.transpose(a_matrix[0,:,:self.T]), decimals=3)

            #Gets matrices for the disparities of critical household conditions and constraints
            Chained_C_Condition, Modified_Budget_Constraint, Modified_Budget_Constraint2, Household_Euler = check_household_conditions(w_path, r_path, c_matrix, a_matrix, psi, bqvec_path)
            
            #Prints if each set of conditions are satisfied or not
            if Print_HH_Eulers:
                print "\nEuler Household satisfied:", np.isclose(np.max(np.absolute(Household_Euler)), 0)
                print "Equation 3.22 satisfied:", np.isclose(np.max(np.absolute(Chained_C_Condition)), 0)
                print "Equation 3.19 satisfied:", np.isclose(np.max(np.absolute(Modified_Budget_Constraint)), 0)
                print "Equation 3.19 (other) satisfied:", np.isclose(np.max(np.absolute(Modified_Budget_Constraint2)), 0)

                #print np.round(np.transpose(100000*Chained_C_Condition[0,:,:self.T]), decimals=4)
                #print np.round(np.transpose(Modified_Budget_Constraint[0,:,:self.T]), decimals=4) #Eulers
                #print np.transpose(Chained_C_Condition[0,:,:self.T])
                #print np.round(np.transpose(self.MortalityRates[0,:,:self.T]), decimals=4)

            #Returns only up until time T and not the vector
            return c_matrix[:,:,:self.T], a_matrix[:,:-1,:self.T]

        #Equation 3.25
        w_path = np.einsum("it,i->it",np.einsum("i,t->it",self.alpha*self.A,1/r_path)**(self.alpha/(1-self.alpha)),(1-self.alpha)*self.A)

        #Equation 3.21
        psi = self.get_Psi(w_path,self.e)

        #Equations 3.19, 3.22
        c_matrix, a_matrix = get_c_a_matrices(w_path, r_path, psi, bqvec_path, Print_HH_Eulers, Print_cabqTimepaths)

        #Equation 3.20
        lhat_path = self.get_lhat(c_matrix, w_path[:,:self.T], self.e[:,:,:self.T])

        #Equation 3.14
        n_path = self.get_n(lhat_path)

        #Equation 3.26
        kd_path = np.sum(a_matrix*self.Nhat[:,:,:self.T],axis=1)

        #Equation 3.15
        y_path = self.get_Y(kd_path,n_path)

        #Equation 3.27
        kf_path = np.outer(self.alpha*self.A, 1/r_path[:self.T])**( 1/(1-self.alpha) )*n_path - kd_path

        return w_path, c_matrix, a_matrix, kd_path, kf_path, n_path, y_path, lhat_path

    def EulerSystemTPI(self, guess, Print_HH_Eulers, Print_cabqTimepaths):
        """
        Description:
            -Description of the Function

        Inputs:
            - guess
            - Print_HH_Eulers
            - Print_cabqTimepaths

        Variables Called from Object:
            - self.T                 = Int: Number of time periods
            - self.I                 = Int: Number of Countries
            - self.S                     = Int: Number of Cohorts
            - self.FirstDyingAge     = Int: First age where mortality rates effect agents
            - self.FirstFertilityAge = Int: First age where agents give birth
            - self.Nhat
            - self.MortalityRates
            - self.Timepath_counter
            - self.IterationsToShow
            -

        Variables Stored in Object:
            - None

        Other Functions Called:
            - self.GetTPIComponents
            - self.plot_timepaths

        Objects in Function:
            - r_path
            - bq_path
            - bqvec_path
            - w_path
            - c_matrix
            - a_matrix
            - kd_path
            - kf_path
            - n_path
            - y_path
            - lhat_path
            - alldeadagent_assets
            - Euler_bq
            - Euler_kf

        Outputs:
            - Euler_all
        """

        guess = np.expand_dims(guess, axis=1).reshape((self.I+1,self.T))
        r_path = guess[0,:]
        bq_path = guess[1:,:]

        r_path = np.hstack((r_path, np.ones(self.S)*self.r_ss))
        bq_path = np.column_stack((  bq_path,   np.outer(self.bq_ss,np.ones(self.S))  ))

        bqvec_path = np.zeros((self.I,self.S,self.T+self.S))
        bqvec_path[:,self.FirstFertilityAge:self.FirstDyingAge,:] = np.einsum("it,s->ist", bq_path, \
                np.ones(self.FirstDyingAge-self.FirstFertilityAge))

        w_path, c_matrix, a_matrix, kd_path, \
        kf_path, n_path, y_path, lhat_path = self.GetTPIComponents(bqvec_path, r_path,Print_HH_Eulers, Print_cabqTimepaths)

        alldeadagent_assets = np.sum(a_matrix[:,self.FirstDyingAge:,:]*\
                self.MortalityRates[:,self.FirstDyingAge:,:self.T]*self.Nhat[:,self.FirstDyingAge:,:self.T], axis=1)

        Euler_bq = bq_path[:,:self.T] - alldeadagent_assets/np.sum(self.Nhat[:,self.FirstFertilityAge:self.FirstDyingAge,:self.T],\
                axis=1)

        Euler_kf = np.sum(kf_path,axis=0)

        Euler_all = np.append(Euler_bq, Euler_kf)

        if self.Iterate: 
            print "Iteration:", self.Timepath_counter, "Min Euler:", np.min(np.absolute(Euler_all)), "Mean Euler:", np.mean(np.absolute(Euler_all)), "Max Euler_bq:", np.max(np.absolute(Euler_bq)), "Max Euler_kf", np.max(np.absolute(Euler_kf))

        if self.Timepath_counter in self.IterationsToShow:
            self.plot_timepaths(SAVE=False, Paths = (r_path, bq_path, w_path, c_matrix, lhat_path, n_path, kd_path, kf_path))

        self.Timepath_counter += 1
        
        return Euler_all

    def Timepath_fsolve(self, Print_HH_Eulers, Print_cabqTimepaths, to_plot = set([])):
        """
        Description:
            -Description of the Function

        Inputs:
            - Print_HH_Eulers
            - Print_cabqTimepaths
            - to_plot

        Variables Called from Object:
            - self.S                     = Int: Number of Cohorts
            - self.r_ss
            - self.I                 = Int: Number of Countries
            - self.T                 = Int: Number of time periods
            - self.bq_ss
            - self.FirstFertilityAge = Int: First age where agents give birth
            - self.FirstDyingAge     = Int: First age where mortality rates effect agents
            - self.bqvec_path
            - self.r_path
            - self.IterationsToShow

        Variables Stored in Object:
            - self.r_path
            - self.bq_path
            - self.bqvec_path
            - self.w_path
            - self.c_matrix
            - self.a_matrix
            - self.kd_path
            - self.kf_path
            - self.n_path
            - self.y_path
            - self.lhat_path

        Other Functions Called:
            - self.get_initialguesses()
            - self.EulerSystemTPI
            - self.GetTPIComponents

        Objects in Function:
            - guess
            - paths
            - r_path
            - bq_path

        Outputs:
            - None
        """

        
        self.IterationsToShow = to_plot

        rpath_guess, bqpath_guess = self.get_initialguesses()

        guess = np.append(rpath_guess, bqpath_guess)

        paths = opt.fsolve(self.EulerSystemTPI, guess, args=(Print_HH_Eulers, Print_cabqTimepaths) )

        paths = np.expand_dims(paths, axis=1).reshape((self.I+1,self.T))
        r_path = paths[0,:]
        bq_path = paths[1:,:]           
        
        self.r_path = np.hstack((r_path, np.ones(self.S)*self.r_ss))

        self.bq_path = np.column_stack(( bq_path,  np.outer(self.bq_ss,np.ones(self.S)) ))
        self.bqvec_path = np.zeros((self.I,self.S,self.T+self.S))
        self.bqvec_path[:,self.FirstFertilityAge:self.FirstDyingAge,:] = np.einsum("it,s->ist", self.bq_path, \
                np.ones(self.FirstDyingAge-self.FirstFertilityAge))

        self.w_path, self.c_matrix, self.a_matrix, self.kd_path, self.kf_path, self.n_path, self.y_path, self.lhat_path = \
                self.GetTPIComponents(self.bqvec_path, self.r_path, Print_HH_Eulers, Print_cabqTimepaths)

    def plot_timepaths(self, SAVE=False, Paths = None):
        """
        Description:
        - Take the timepaths and plots them into one sheet of graphs
            TODO: Fill in the array sizes
            TODO: Update the inputs so it is one tuple

        Inputs:
            - r_path                = Array:[], Given interest rate path
            - bq_path               = Array:[], Given bequests path
            - c_matrix              = Array:[], Given consumption matrix
            - lhat_path             = Array:[], Given time endowment
            - n_path                = Array:[], Given aggregate labor productivity
            - kd_path               = Array:[], Given domestic capital path
            - kf_path               = Array:[], Given foreign capital path
            - SAVE                  = Boolean:[], Switch that determines whether we save the graphs or simply show it.

        Variables Called from Object:
            - self.S                = Int: Number of Cohorts
            - self.T                = Int: Number of time periods
            - self.I                = Int: Number of Countries
            - self.Timepath_counter = Int: Counter that keeps track of the number of iterations in solving for the time paths
            - self.I_touse          = List: [I], Roster of countries that are being used
            - self.cvec_ss          = Array:[],
            - self.lhat_ss          = Array:[],
            - self.n_ss             = Array:[],
            - self.kd_ss            = Array:[],

        Variables Stored in Object:
            - None

        Other Functions Called:
            - None

        Objects in Function:
            - title                 = String: Overall title of the sheet of graphs
            - name                  = String: Name of the .png file that will save the graphs.
            - ax                    = String: TODO

        Outputs:
            - None

        """
        if Paths is None:
            r_path, bq_path, w_path, c_matrix, lhat_path, n_path, kd_path, kf_path = \
            self.r_path, self.bq_path, self.w_path, self.c_matrix, self.lhat_path, self.n_path, self.kd_path, self.kf_path
        else:
            r_path, bq_path, w_path, c_matrix, lhat_path, n_path, kd_path, kf_path = Paths


        title = str("S = " + str(self.S) + ", T = " + str(self.T))
        plt.suptitle(title)

        ax = plt.subplot(331)
        for i in range(self.I):
            plt.plot(range(self.S+self.T), r_path)
        plt.title(str("r_path "+"iteration: "+str(self.Timepath_counter)))
        plt.legend(self.I_touse)

        plt.subplot(332)
        for i in range(self.I):
            plt.plot(range(self.S+self.T), bq_path[i,:])
        plt.title(str("bqvec_path "+"iteration: "+str(self.Timepath_counter)))

        plt.subplot(333)
        for i in range(self.I):
            plt.plot(range(self.S+self.T), w_path[i,:])
        plt.title(str("w_path "+"iteration: "+str(self.Timepath_counter)))

        plt.subplot(334)
        for i in range(self.I):
            plt.plot( range(self.S+self.T), np.hstack((np.sum(c_matrix[i,:,:],axis=0),np.ones(self.S)*np.sum(self.cvec_ss[i,:]))) )
        plt.title(str("C_path "+"iteration: "+str(self.Timepath_counter)))

        plt.subplot(335)
        for i in range(self.I):
            plt.plot( range(self.S+self.T), np.hstack((np.sum(lhat_path[i,:,:],axis=0),np.ones(self.S)*np.sum(self.lhat_ss[i,:]))) )
        plt.title(str("Lhat_path "+"iteration: "+str(self.Timepath_counter)))

        plt.subplot(336)
        for i in range(self.I):
            plt.plot(range(self.S+self.T), np.hstack((n_path[i,:],np.ones(self.S)*self.n_ss[i])))
        plt.title(str("n_path "+"iteration: "+str(self.Timepath_counter)))

        plt.subplot(337)
        for i in range(self.I):
            plt.plot(range(self.S+self.T), np.hstack((kd_path[i,:],np.ones(self.S)*self.kd_ss[i])) )
        plt.title(str("kd_path "+"iteration: "+str(self.Timepath_counter)))
        
        plt.subplot(338)
        for i in range(self.I):
            plt.plot(range(self.S+self.T), np.hstack((kf_path[i,:],np.ones(self.S)*self.kf_ss[i])))
        plt.title(str("kf_path "+"iteration: "+str(self.Timepath_counter)))

        plt.subplot(339)
        for i in range(self.I):
            plt.plot(range(self.S+self.T), np.hstack((kf_path[i,:]+kd_path[i,:],np.ones(self.S)*(self.kf_ss[i]+self.kd_ss[i]))))
        plt.title(str("K_path "+"iteration: "+str(self.Timepath_counter)))


        if SAVE:
            name= "Graphs/OLGresult_Iter"+str(self.Timepath_counter)+"_"+str(self.I)+"_"+str(self.S)+"_"+str(self.sigma)+".png"
            plt.savefig(name)
            plt.clf()

        else:
            plt.show()

