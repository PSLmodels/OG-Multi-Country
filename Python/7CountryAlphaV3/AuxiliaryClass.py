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
    (NEW!) This object centralizes all of the operations of the model. Before, we had to pass in and keep track of
    different parameters. With this, we have all important parts of the model to the OLG object, which can be easily
    accessed by any of the functions.    
    Unlike the previous versions of the code, the comments for this version will indicate inputs, outputs, and which
    objects will be stored in the objects.


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

    def __init__(self, countries, HH_Params, Firm_Params, Lever_Params, Tol_Params):
        """
        Description: 
            -This creates the object and stores all of the parameters into the object.
             The initialization is the starting point for model.

        Inputs:
            -self: "self" stores all of the components of the model. To access any part,
             simply type "self.variable_name" while in the object and "objectname.variable_name"
             outside the object. Every other object function will just take this as given, so 
             future mentions of self won't be rewritten.

            -countries              = tuple: contains a dictionary and tuple for countries and their associated number
            -Firm_Params            = tuple: contains alpha, annualized delta, chi, rho and g_A
            -HH_Params              = tuple: contains S, I, annualized Beta and sigma.
            -Lever_Params           = tuple: contains boolean levers indicated by the users such as:
                                        CalcTPI,PrintAges,PrintLoc,PrintSSEulErrors,PrintSS,ShowSSGraphs,Print_cabqTimepaths
                                        Print_HH_Eulers,CheckerMode,Iterate,DemogGraphs,TPIGraphs,UseStaggeredAges,
                                        UseDiffDemog,UseSSDemog,UseDiffProductivities,UseTape,ADJUSTKOREAIMMIGRATION,
                                        VectorizeHouseholdSolver,PinInitialValues,UsePrev_C0
            -Tol_Params             = tuple: contains the xi parameter and the maximum number of iterations.

        Variables Stored in Object:

            - self.A                = Array: [I,1], Technology level for each country
            - self.c0_alive         = Array:
            - self.c0_future        = Array:
            - self.e                = Array: [I,S,T], Labor Productivities
            - self.I_touse          = Array: [I], Roster of countries that are being used
            - self.lbar             = Array:
            - self.MortalityRates   = Array:
            - self.Mortality_ss     = Array:
            - self.Nhat             = Array:
            - self.Nhat_ss          = Array:
            - self.CalcTPI          = Boolean: Activates Calculating the TPI
            - self.CheckerMode      = Boolean: Used in conjunction with Checker.py, a MPI code that checks the
                                               robustness of the code. With this activated, the code only prints
                                               the statements that are necessary. This speeds up the robust check
                                               process.
            - self.DemogGraphs      = Boolean: Activates the showing of the deomgraphics graphs
            - self.Iterate          = Boolean: Activates printing the iteration number and euler errors at each
                                               step of the TPI process.
            - self.PinInitialValues = Boolean: REDUNDANT, WILL REMOVE SOON
            - self.PrintAges        = Boolean: Prints the ages calculated in the demographics
            - self.Print_cabqTimepaths = Boolean: Prints the assests and consumption matrices we are filling
            - self.Print_HH_Eulers  = Boolean: Prints the Euler Errors for households
            - self.PrintLoc         = Boolean: Prints the location of the code, used for debugging purposes
            - self.PrintSS          = Boolean: Prints the results of the steady state calculation
            - self.PrintSSEulErrors = Boolean: Prints the Euler Errors calculated in the steady state
            - self.ShowSSGraphs     = Boolean: Activates showing the graphs that result from the steady state
            - self.TPIGraphs        = Boolean: Activates the final TPI graph
            - self.UseDiffDemog     = Boolean: Allows each country to have different demographics.
            - self.UsePrev_c0       = Boolean:
            - self.UseStaggeredAges = Boolean:
            - self.I_dict           = Dictionary: [I], Associates a country with a number

            - self.agestopull       = Scalar:
            - self.alpha            = Scalar: Capital share of production
            - self.beta             = Scalar: calculated overall future discount rate
            - self.chi              = Scalar:
            - self.delta            = Scalar: calulated overall depreciation rate
            - self.demog_ss_tol     = Scalar:
            - self.FirstDyingAge    = Scalar: From the Auxiliary Demographics module, see that page for details
            - self.FirstFertilityAge= Scalar: From the Auxiliary Demographics module, see that page for details
            - self.g_A              = Scalar:
            - self.I                = Scalar: Number of Countries
            - self.LastFertilityAge = Scalar: From the Auxiliary Demographics module, see that page for details
            - self.LeaveHouseAge    = Scalar: From the Auxiliary Demographics module, see that page for details

            - self.MaxImmigrantAge  = Scalar: From the Auxiliary Demographics module, see that page for details
            - self.rho              = Scalar:
            - self.S                = Scalar: Number of Cohorts
            - self.T                = Scalar: of the total amount of time periods
            - self.T_1              = Scalar: Transition year for the demographics
            - self.Timepath_counter = Scalar:
            - self.IterationsToShow = Set: A set of user inputs of iteration TPI graphs to show 

        Objects in Function:
            - beta_annual           = Scalar:
            - delta_annual          = Scalar:
        """
        #PARAMETER SET UP

        #HH Parameters
        (self.S, self.I, beta_annual,self.sigma) = HH_Params
        
        self.beta=beta_annual**(70/self.S)

        self.T = int(round(4*self.S))

        self.T_1 = self.S

        if self.S > 50:
            self.T_1 = 50

        #Demographics Parameters

        self.I_dict, self.I_touse = countries
        self.Nhat = np.ones((self.I, self.S, self.T))
        self.Nhat_ss = np.ones((self.I, self.S))
        self.MortalityRates = np.zeros((self.I, self.S, self.T+self.S))
        self.Mortality_ss = np.zeros((self.I, self.S))


        #Firm Parameters
        (self.alpha,delta_annual,self.chi,self.rho, self.g_A)= Firm_Params
        self.delta=1-(1-delta_annual)**(70/self.S)

        #Lever Parameters
        (self.CalcTPI,self.PrintAges,self.PrintLoc,self.PrintSSEulErrors,self.PrintSS,self.ShowSSGraphs,self.Print_cabqTimepaths,\
         self.Print_HH_Eulers, self.CheckerMode,self.Iterate,self.DemogGraphs,self.TPIGraphs,self.UseStaggeredAges,self.UseDiffDemog,\
         self.UseSSDemog,self.UseDiffProductivities,self.UseTape,self.ADJUSTKOREAIMMIGRATION, self.VectorizeHouseholdSolver,\
         self.PinInitialValues,self.UsePrev_c0) = Lever_Params

        self.IterationsToShow = set([])

        #Tolerance Parameters

        (self.demog_ss_tol) = Tol_Params

        self.LeaveHouseAge, self.FirstFertilityAge, self.LastFertilityAge, self.MaxImmigrantAge, self.FirstDyingAge,\
                self.agestopull = demog.getkeyages(self.S,self.PrintAges,self.UseStaggeredAges)

        if self.UseDiffDemog:
            self.A = np.ones(self.I)+np.cumsum(np.ones(self.I)*.05)-.05 #Techonological Change, used for when countries are different

        else:
            self.A = np.ones(self.I) #Techonological Change, used for idential countries

        if self.UseDiffProductivities:
            self.e = np.ones((self.I, self.S, self.T+self.S))
            self.e[:,self.FirstDyingAge:,:] = 0.3
            self.e[:,:self.LeaveHouseAge,:] = 0.3
        else:
            self.e = np.ones((self.I, self.S, self.T+self.S)) #Labor productivities

        self.e_ss=self.e[:,:,-1]

        self.lbar = np.cumsum(np.ones(self.T+self.S)*self.g_A)
        self.lbar[self.T:] = np.ones(self.S)
        self.lbar[:self.T] = np.ones(self.T)
        self.lbar_ss=self.lbar[-1]

        self.Timepath_counter = 1

        self.c0_alive = np.ones((self.I, self.S-1))*.3

        self.c0_future = np.ones((self.I,self.T))*.3

    #DEMOGRAPHICS SET-UP

    def Import_Data(self):
        """
        Description:
            - This function activates importing the .CSV files that contain our demographics.

        Variables Called from Object:
            - self.S                = Scalar: Number of Cohorts
            - self.T                = Scalar: of the total amount of time periods
            - self.I                = Scalar: Number of Countries



        Variables Stored in Object:
            - self.N                = Scalar: Number of Countries
            - self.Nhat             = Scalar: Number of Countries
            - self.all_FertilityAges= Scalar: Number of Countries
            - self.FertilityRates   = Scalar: Number of Countries
            - self.MortalityRates   = Scalar: Number of Countries
            - self.Migrants         = Scalar: Number of Countries
            - self.g_N              = Scalar: Number of Countries


        Other Functions Called:
            -

        Objects in Function:
            - f_range              = Scalar: Range of fertilities, based on start age
            - I_all

        Outputs:
            -

        """

        f_range=self.LastFertilityAge+1-self.FirstFertilityAge


        self.N=np.zeros((self.I,self.S,self.T))
        self.Nhat=np.zeros((self.I,self.S,self.T))
        self.all_FertilityRates = np.zeros((self.I, self.S, f_range+self.T))
        self.FertilityRates = np.zeros((self.I, self.S, self.T))
        self.MortalityRates = np.zeros((self.I, self.S, self.T))
        self.Migrants = np.zeros((self.I, self.S, self.T))
        self.g_N = np.zeros(self.T)

        I_all = ["usa", "eu", "japan", "china", "india", "russia", "korea"]

        for i in xrange(self.I):

            if self.UseDiffDemog:
                index = I_all.index(self.I_touse[i])
                if self.I > len(I_all):
                    sys.exit("Error! There can't be more than", len(I_all),"countries without\
                            unique data. \n You must change I so it's less than", len(I_all),\
                            "or change DiffDemog to False")
            else:
                index = 0


            if self.UseStaggeredAges:
                self.N[i,:,0] = np.loadtxt(("Data_Files/population.csv"),delimiter=',',\
                        skiprows=1, usecols=[index+1])[self.agestopull]*1000

                self.all_FertilityRates[i,self.FirstFertilityAge:self.LastFertilityAge+1,\
                        :f_range+self.T_1] =  np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_fertility.csv"),delimiter=',',skiprows=1\
                        , usecols=(self.agestopull[self.FirstFertilityAge:self.LastFertilityAge+1]-22))[48-f_range:48+self.T_1,:])

                self.MortalityRates[i,self.FirstDyingAge:,:self.T_1] = np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_mortality.csv"),delimiter=','\
                        ,skiprows=1, usecols=(self.agestopull[self.FirstDyingAge:]-67))[:self.T_1,:])

                self.Migrants[i,:self.MaxImmigrantAge,:self.T_1] = np.einsum("s,t->st",np.loadtxt(("Data_Files/net_migration.csv"),delimiter=','\
                        ,skiprows=1, usecols=[index+1])[self.agestopull[:self.MaxImmigrantAge]]*100, np.ones(self.T_1))

            else:
                self.N[i,:,0] = np.loadtxt(("Data_Files/population.csv"),delimiter=',',skiprows=1, usecols=[index+1])[:self.S]*1000

                self.all_FertilityRates[i,self.FirstFertilityAge:self.LastFertilityAge+1,:f_range+self.T_1] = \
                        np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_fertility.csv"),delimiter=','\
                        ,skiprows=1, usecols=range(1,f_range+1))[48-f_range:48+self.T_1,:])

                self.MortalityRates[i,self.FirstDyingAge:-1,:self.T_1] = \
                        np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_mortality.csv"),delimiter=','\
                        ,skiprows=1, usecols=range(1,self.S-self.FirstDyingAge))[:self.T_1,:])

                self.Migrants[i,:self.MaxImmigrantAge,:self.T_1] = np.einsum("s,t->st",np.loadtxt(("Data_Files/net_migration.csv"),delimiter=','\
                        ,skiprows=1, usecols=[index+1])[:self.MaxImmigrantAge]*100, np.ones(self.T_1))

            if self.PrintLoc: print "Got Demographics for", I_all[index]

            if self.ADJUSTKOREAIMMIGRATION and I_all[index] == "korea":
                self.Migrants[i,:]/=100

        #Gets initial population share
        self.Nhat[:,:,0] = self.N[:,:,0]/np.sum(self.N[:,:,0])

        #The last generation dies with probability 1
        self.MortalityRates[:,-1,:] = np.ones((self.I, self.T))

        #Gets steady-state values for all countries by taking the mean at year T_1-1 across countries
        self.f_bar = np.mean(self.all_FertilityRates[:,:,f_range+self.T_1-1], axis=0)
        self.rho_bar = np.mean(self.MortalityRates[:,:,self.T_1-1], axis=0)

        #Set to the steady state for every year beyond year T_1
        self.all_FertilityRates[:,:,f_range+self.T_1:] = np.tile(np.expand_dims(self.f_bar, axis=2), (self.I,1,self.T-self.T_1))
        self.MortalityRates[:,:,self.T_1:] = np.tile(np.expand_dims(self.rho_bar, axis=2), (self.I,1,self.T-self.T_1))

        #FertilityRates is exactly like all_FertilityRates except it begins at time t=0 rather than time t=-f_range
        self.FertilityRates[:,self.FirstFertilityAge:self.LastFertilityAge+1,:] = self.all_FertilityRates[:,self.FirstFertilityAge:self.LastFertilityAge+1,f_range:]

        #Gets initial world population growth rate
        self.g_N[0] = 0.

    def Demographics(self):

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




        self.ImmigrationRates = np.zeros((self.I,self.S,self.T))

        N_temp = np.ones((self.I,self.S))/(self.I*self.S)

        for t in xrange(1,self.T):

            self.N[:,0,t] = np.sum((self.N[:,:,t-1]*self.FertilityRates[:,:,t-1]), axis=1)
            N_temp[:,0] = np.sum((self.Nhat[:,:,t-1]*self.FertilityRates[:,:,t-1]), axis=1)


            if t <= self.T_1:
                self.ImmigrationRates[:,:,t-1] = self.Migrants[:,:,t-1]/self.N[:,:,t-1]

            else:
                self.ImmigrationRates[:,:,t-1] = np.mean(self.ImmigrationRates[:,:,self.T_1-1],\
                        axis=0)

            self.N[:,1:,t] = self.N[:,:-1,t-1]*(1+self.ImmigrationRates[:,:-1,t-1]-\
                    self.MortalityRates[:,:-1,t-1])
            N_temp[:,1:] = self.Nhat[:,:-1,t-1]*(1+self.ImmigrationRates[:,:-1,t-1]-\
                    self.MortalityRates[:,:-1,t-1])
            
            self.Nhat[:,:,t] = self.N[:,:,t]/np.sum(self.N[:,:,t])

            self.g_N[t] = np.sum(N_temp[:,:])-1


        self.ImmigrationRates[:,:,t] = self.Migrants[:,:,t]/self.N[:,:,t]

        pop_old = self.N[:,:,-1]
        pop_new = self.N[:,:,-1]

        iteration = 0

        while np.max(np.abs(self.Nhat[:,:,-1] - self.Nhat[:,:,-2])) > self.demog_ss_tol:
            pop_new[:,0] = np.sum((pop_old[:,:]*self.FertilityRates[:,:,-1]),axis=1)
            pop_new[:,1:] = pop_old[:,:-1]*(1+self.ImmigrationRates[:,:-1,-1]\
                    -self.MortalityRates[:,:-1,-1])
            self.Nhat = np.dstack((self.Nhat,pop_new/np.sum(pop_new)))
            iteration += 1

        if self.PrintLoc: print "The SS Population Share converged in", iter, "years beyond year T"


        if self.CheckerMode==False:
            print "\nDemographics obtained!"


        self.Nhat_ss = self.Nhat[:,:,-1]
        self.Nhat = self.Nhat[:,:,:self.T]

        #Imposing the ss for years after self.T
        self.Nhat = np.dstack((  self.Nhat[:,:,:self.T], np.einsum("is,t->ist",self.Nhat_ss,np.ones(self.S))  ))

        self.Mortality_ss=self.MortalityRates[:,:,-1]
        #Imposing the ss for years after self.T
        self.MortalityRates = np.dstack((  self.MortalityRates[:,:,:self.T], np.einsum("is,t->ist",self.Mortality_ss,np.ones(self.S))  ))        

        if self.UseSSDemog == True:
            self.Nhat = np.einsum("is,t->ist",self.Nhat_ss,np.ones(self.T+self.S))
            self.MortalityRates = np.einsum("is,t->ist",self.Mortality_ss,np.ones(self.T+self.S))

        if self.DemogGraphs:
            ages = self.FirstFertilityAge, self.LastFertilityAge, self.FirstDyingAge, \
                    self.MaxImmigrantAge
            datasets = self.FertilityRates, self.MortalityRates, self.ImmigrationRates, self.Nhat
            demog.plotDemographics(ages, datasets, self.I, self.S, self.T, self.I_touse, T_touse = [0,1,2,3,20]\
                    , compare_across="T", data_year=0)

    #STEADY STATE

    def get_Psi(self, w, e):

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

        if e.ndim == 2:
            we =  np.einsum("i,is->is",w,e)

        elif e.ndim == 3:
            we = np.einsum("it, ist -> ist", w, e)

        part1 = (self.chi/we)**(self.rho-1)

        psi = (1+self.chi*part1)**( (1-self.rho*self.sigma)/(self.rho-1) )

        return psi

    def get_lhat(self,c,w,e):

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


        if e.ndim == 2:
            lhat=c*(self.chi/np.einsum("i,is->is",w,e))**self.rho
        elif e.ndim == 3:
            lhat=c*(self.chi/np.einsum("it,ist->ist",w,e))**self.rho

        return lhat

    def get_n(self, lhat):
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

        if lhat.ndim == 2:
            n = np.sum(self.e_ss*(self.lbar_ss-lhat)*self.Nhat_ss,axis=1)
        elif lhat.ndim == 3:
            n = np.sum(self.e[:,:,:self.T]*(self.lbar[:self.T]-lhat)*self.Nhat[:,:,:self.T],axis=1)

        return n

    def get_Y(self, kd, n):
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

        if kd.ndim ==1:
            Y = (kd**self.alpha) * ((self.A*n)**(1-self.alpha))
        elif kd.ndim== 2:
            Y = (kd**self.alpha) * (np.einsum("i,is->is",self.A,n)**(1-self.alpha))

        return Y

    def GetSSComponents(self, bq_ss, r_ss):
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

        def get_lifetime_decisionsSS(c_1, w_ss, r_ss):
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


            cvec_ss = np.zeros((self.I,self.S))
            avec_ss = np.zeros((self.I,self.S+1))
            cvec_ss[:,0] = c_1

            for s in xrange(self.S-1):
                cvec_ss[:,s+1] = (self.beta * (1-self.Mortality_ss[:,s]) * (1 + r_ss - self.delta)\
                        *self.psi_ss[:,s+1]/self.psi_ss[:,s])**(1/self.sigma) * cvec_ss[:,s]*np.exp(-self.g_A)

                avec_ss[:,s+1] = (w_ss*self.e_ss[:,s] + (1 + r_ss - self.delta)*avec_ss[:,s] + \
                        bq_ss[:,s] - cvec_ss[:,s]*(1+w_ss*self.e_ss[:,s]*\
                        (self.chi/(w_ss*self.e_ss[:,s]))**self.rho))*np.exp(-self.g_A)

            avec_ss[:,s+2] = (w_ss*self.e_ss[:,s+1] + (1 + r_ss - self.delta)*avec_ss[:,s+1] \
                    - cvec_ss[:,s+1]*(1+w_ss*self.e_ss[:,s+1]*(self.chi/(w_ss*self.e_ss[:,s+1]))\
                    **self.rho))*np.exp(-self.g_A)

            return cvec_ss, avec_ss


        def householdEuler_SS(c_1, w_ss, r_ss):
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


            cpath, assets_path = get_lifetime_decisionsSS(c_1, w_ss, r_ss)

            Euler = np.ravel(assets_path[:,-1])

            if np.any(cpath<0):
                print "WARNING! The fsolve for initial optimal consumption guessed a negative number"
                Euler = np.ones(Euler.shape[0])*9999.

            return Euler

        w_ss = (self.alpha*self.A/r_ss)**(self.alpha/(1-self.alpha))*(1-self.alpha)*self.A

        self.psi_ss = self.get_Psi(w_ss,self.e_ss)

        c1_guess = np.ones(self.I)*.02

        opt_c1 = opt.fsolve(householdEuler_SS, c1_guess, args = (w_ss, r_ss))

        cvec_ss, avec_ss = get_lifetime_decisionsSS(opt_c1,w_ss,r_ss)

        avec_ss = avec_ss[:,:-1]

        lhat_ss = self.get_lhat(cvec_ss, w_ss, self.e_ss)

        n_ss = self.get_n(lhat_ss)

        kd_ss = np.sum(avec_ss*self.Nhat_ss,axis=1)
        y_ss = self.get_Y(kd_ss,n_ss)

        kf_ss = (self.alpha*self.A/r_ss)**(1/(1-self.alpha)) * n_ss-kd_ss

        K_ss_with_tape = np.clip(kd_ss + kf_ss, .0001, np.max(kd_ss + kf_ss))

        return w_ss, cvec_ss, avec_ss, kd_ss, kf_ss, n_ss, y_ss, lhat_ss

    def EulerSystemSS(self,guess):
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


        bq_ss = guess[:-1]
        r_ss = guess[-1]

        bqvec_ss = np.zeros((self.I,self.S))
        bqvec_ss[:,self.FirstFertilityAge:self.FirstDyingAge] = np.einsum("i,s->is", bq_ss, \
                np.ones(self.FirstDyingAge-self.FirstFertilityAge))

        w_ss, cvec_ss, avec_ss, kd_ss, kf_ss, n_ss, y_ss, lhat_ss = self.GetSSComponents(bqvec_ss, r_ss)


        alldeadagent_assets = np.sum(avec_ss[:,self.FirstDyingAge:]*\
                self.Mortality_ss[:,self.FirstDyingAge:]*self.Nhat_ss[:,self.FirstDyingAge:], axis=1)

        Euler_bq = bq_ss - alldeadagent_assets/np.sum(self.Nhat_ss[:,self.FirstFertilityAge:self.FirstDyingAge],\
                axis=1)

        Euler_kf = np.sum(kf_ss)

        Euler_all = np.append(Euler_bq, Euler_kf)

        if self.PrintSSEulErrors: print "Euler Errors:", Euler_all

        return Euler_all

    def SteadyState(self, rss_guess, bqss_guess):
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

        ss = opt.fsolve(self.EulerSystemSS, guess)

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
        print "-Euler for bq satisfied:", np.isclose(np.max(np.absolute(Euler_bq)), 0)
        print "-Euler for r satisfied:", np.isclose(Euler_kf, 0), "\n\n"

        if self.PrintSS:
            if self.ShowSSGraphs:
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

            print "assets steady state:", self.avec_ss
            print "kf steady state", self.kf_ss
            print "kd steady state", self.kd_ss
            print "n steady state", self.n_ss
            print "y steady state", self.y_ss
            print "r steady state", self.r_ss
            print "w steady state", self.w_ss
            print "c_vec_ss steady state", self.cvec_ss

    def checkSSEulers(self):
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

        we = np.einsum("i,is->is",self.w_ss,self.e_ss[:,:-1])

        print self.psi_ss[:,:-1]*self.cvec_ss[:,:-1]**(-self.sigma) - self.beta*(1-self.Mortality_ss[:,:-1])*self.psi_ss[:,1:]*(self.cvec_ss[:,1:]*np.exp(self.g_A))**(-self.sigma)*(1+self.r_ss-self.delta)
        
        print self.cvec_ss[:,:-1] - \
        (we + (1+self.r_ss-self.delta)*self.avec_ss[:,:-1] + self.bqvec_ss[:,:-1] - self.avec_ss[:,1:]*np.exp(self.g_A)) / \
        (1 + we*(self.chi/we)**self.rho)

    #TIMEPATH-ITERATION

    def set_initial_values(self, r_init, bq_init, a_init):
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

        self.r_init = r_init
        self.bq_init = bq_init
        self.a_init = a_init

    def get_initialguesses(self):
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


        rpath_guess = np.zeros(self.T)
        bqpath_guess = np.zeros((self.I,self.T))

        cc = self.r_init
        bb = -2 * (self.r_init-self.r_ss)/(self.T-1)
        aa = -bb / (2*(self.T-1))
        rpath_guess[:self.T] = aa * np.arange(0,self.T)**2 + bb*np.arange(0,self.T) + cc

        for i in range(self.I):
            bqpath_guess[i,:self.T] = np.linspace(self.bq_init[i], self.bq_ss[i], self.T)

        return rpath_guess, bqpath_guess

    def GetTPIComponents(self, bqvec_path, r_path):

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


        #Functions that solve lower-diagonal household decisions in vectors
        def get_lifetime_decisions_LOWERTRIANGLETEST(c0_guess, c_uppermat, a_uppermat, w_path, r_path, psi, bqvec_path):
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
        def get_c_a_matrices(w_path, r_path, psi, bqvec_path):
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

            
            #Initializes the consumption and assets matrices
            c_matrix = np.zeros((self.I,self.S,self.T+self.S))
            a_matrix = np.zeros((self.I,self.S+1,self.T+self.S))
            a_matrix[:,:-1,0] = self.a_init

            #Equation 3.19 for the oldest agent in time t=0. Note that this agent chooses to consume everything so that it has no assets in the following period
            c_matrix[:,self.S-1,0] = (w_path[:,0]*self.e[:,self.S-1,0] + (1 + r_path[0] - self.delta)*self.a_init[:,self.S-1] + bqvec_path[:,self.S-1,0])\
            /(1+w_path[:,0]*self.e[:,self.S-1,0]*(self.chi/(w_path[:,0]*self.e[:,self.S-1,0]))**(self.rho))

            #Will solve the household matrices using vectorization if = True and by agent if = False
            if self.VectorizeHouseholdSolver:
            
                if self.UsePrev_c0:
                    c0alive_guess = self.c0_alive
                else:
                    c0alive_guess = np.ones((self.I, self.S-1))*.3

                opt.fsolve(get_upper_triangle_Euler_TEST, c0alive_guess, args=(c_matrix, a_matrix, w_path, r_path, psi, bqvec_path))

                #Initializes a guess for the first vector for the fsolve to use
                if self.UsePrev_c0:
                    c0future_guess = self.c0_future
                else:
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
                    if self.Print_cabqTimepaths:
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
                    if self.Print_cabqTimepaths:
                        print "Consumption for year", t
                        print np.round(np.transpose(c_matrix[0,:,:self.T]), decimals=3)
                        print "Assets for year", t
                        print np.round(np.transpose(a_matrix[0,:,:self.T]), decimals=3)

            #Gets matrices for the disparities of critical household conditions and constraints
            Chained_C_Condition, Modified_Budget_Constraint, Modified_Budget_Constraint2, Household_Euler = check_household_conditions(w_path, r_path, c_matrix, a_matrix, psi, bqvec_path)
            
            #Prints if each set of conditions are satisfied or not
            if self.Print_HH_Eulers:
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
        c_matrix, a_matrix = get_c_a_matrices(w_path, r_path, psi, bqvec_path)

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

    def EulerSystemTPI(self, guess):
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


        if self.PinInitialValues:
            guess = np.expand_dims(guess, axis=1).reshape((self.I+1,self.T-1))
            r_path = np.zeros(self.T)
            r_path[0] = self.r_init
            r_path[1:] = guess[0,:]

            bq_path = np.zeros((self.I,self.T))
            bq_path[:,0] = self.bq_init
            bq_path[:,1:] = guess[1:,:]

        else:
            guess = np.expand_dims(guess, axis=1).reshape((self.I+1,self.T))
            r_path = guess[0,:]
            bq_path = guess[1:,:]

        r_path = np.hstack((r_path, np.ones(self.S)*self.r_ss))
        bq_path = np.column_stack((  bq_path,   np.outer(self.bq_ss,np.ones(self.S))  ))

        bqvec_path = np.zeros((self.I,self.S,self.T+self.S))
        bqvec_path[:,self.FirstFertilityAge:self.FirstDyingAge,:] = np.einsum("it,s->ist", bq_path, \
                np.ones(self.FirstDyingAge-self.FirstFertilityAge))

        w_path, c_matrix, a_matrix, kd_path, kf_path, n_path, y_path, lhat_path = self.GetTPIComponents(bqvec_path, r_path)

        alldeadagent_assets = np.sum(a_matrix[:,self.FirstDyingAge:,:]*\
                self.MortalityRates[:,self.FirstDyingAge:,:self.T]*self.Nhat[:,self.FirstDyingAge:,:self.T], axis=1)

        Euler_bq = bq_path[:,:self.T] - alldeadagent_assets/np.sum(self.Nhat[:,self.FirstFertilityAge:self.FirstDyingAge,:self.T],\
                axis=1)

        Euler_kf = np.sum(kf_path,axis=0)

        if self.PinInitialValues:
            Euler_all = np.append(Euler_bq[:,1:], Euler_kf[1:])
        else:
            Euler_all = np.append(Euler_bq, Euler_kf)

        if self.Iterate: 
            print "Iteration:", self.Timepath_counter, "Min Euler:", np.min(np.absolute(Euler_all)), "Mean Euler:", np.mean(np.absolute(Euler_all)), "Max Euler_bq:", np.max(np.absolute(Euler_bq)), "Max Euler_kf", np.max(np.absolute(Euler_kf))

        if self.Timepath_counter in self.IterationsToShow:
            self.plot_timepaths(r_path, bq_path, w_path, c_matrix, lhat_path, n_path, kd_path, kf_path, SAVE=False)

        self.Timepath_counter += 1
        
        return Euler_all

    def Timepath_fsolve(self, to_plot = set([])):
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

        
        self.IterationsToShow = to_plot

        rpath_guess, bqpath_guess = self.get_initialguesses()

        if self.PinInitialValues:
            rpath_guess = rpath_guess[1:]
            bqpath_guess = bqpath_guess[:,1:]

        guess = np.append(rpath_guess, bqpath_guess)

        paths = opt.fsolve(self.EulerSystemTPI, guess)

        if self.PinInitialValues:
            paths = np.expand_dims(paths, axis=1).reshape((self.I+1,self.T-1))
            r_path = np.zeros(self.T)
            r_path[0] = self.r_init
            r_path[1:] = paths[0,:]

            bq_path = np.zeros((self.I, self.T))
            bq_path[:,0] = self.bq_init
            bq_path[:,1:] = paths[1:,:]
        else:
            paths = np.expand_dims(paths, axis=1).reshape((self.I+1,self.T))
            r_path = paths[0,:]
            bq_path = paths[1:,:]           
        
        self.r_path = np.hstack((r_path, np.ones(self.S)*self.r_ss))

        self.bq_path = np.column_stack(( bq_path,  np.outer(self.bq_ss,np.ones(self.S)) ))
        self.bqvec_path = np.zeros((self.I,self.S,self.T+self.S))
        self.bqvec_path[:,self.FirstFertilityAge:self.FirstDyingAge,:] = np.einsum("it,s->ist", self.bq_path, \
                np.ones(self.FirstDyingAge-self.FirstFertilityAge))

        self.w_path, self.c_matrix, self.a_matrix, self.kd_path, self.kf_path, self.n_path, self.y_path, self.lhat_path = self.GetTPIComponents(self.bqvec_path, self.r_path)

        self.plot_timepaths(self.r_path, self.bq_path, self.w_path, self.c_matrix, self.lhat_path, self.n_path, self.kd_path, self.kf_path, SAVE=True)

    def plot_timepaths(self, r_path, bq_path, w_path, c_matrix, lhat_path, n_path, kd_path, kf_path, SAVE=False):
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
            name= "OLGresult_Iter"+str(self.Timepath_counter)+"_"+str(self.I)+"_"+str(self.S)+"_"+str(self.sigma)+".png"
            plt.savefig(name)
            plt.clf()

        else:
            plt.show()

