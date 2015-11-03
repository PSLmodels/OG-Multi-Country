from __future__ import division
import csv
import time

import numpy as np
import scipy as sp
import scipy.optimize as opt
from matplotlib import pyplot as plt

import AuxiliaryDemographics as demog
#import AuxiliarySteadyState as SS
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#from mpl_toolkits.plot3d import Axes3D

class OLG(object):


    def __init__(self, countries, HH_Params, Firm_Params, Lever_Params, Tol_Params, SS_Params, TPI_Params):
        #PARAMETER SET UP

        #HH Parameters
        (self.S, self.I, self.beta_annual,self.sigma)=HH_Params
        
        self.beta=self.beta_annual**(70/self.S)
        self.T = int(round(4*self.S))

        self.T_1=self.S

        if self.S > 50:
            self.T_1=50

        #Demographics Parameters

        self.I_dict, self.I_touse = countries

        
        #Firm Parameters
        (self.alpha,self.delta_annual,self.chi,self.rho, self.g_A)= Firm_Params
        self.delta=1-(1-self.delta_annual)**(70/self.S)

        #Lever Parameters
        (self.CalcTPI,self.PrintAges,self.PrintLoc,self.EulErrors,self.PrintSS,self.Print_cabqTimepaths,self.CheckerMode,\
                self.DemogGraphs,self.TPIGraphs,self.UseStaggeredAges,self.UseDiffDemog,self.UseSSDemog,\
                self.UseDiffProductivities,self.UseTape,self.SAVE,self.SHOW,self.ADJUSTKOREAIMMIGRATION) = Lever_Params

        #Tolerance Parameters

        (self.tpi_tol, self.demog_ss_tol) = Tol_Params

        #Steady State Parameters

        #Timepath Iteration Parameters
        #(self.tpi_tol, self.xi, self.MaxIters) = TPI_Params

        self.LeaveHouseAge, self.FirstFertilityAge, self.LastFertilityAge, self.MaxImmigrantAge, self.FirstDyingAge,\
                self.agestopull = demog.getkeyages(self.S,self.PrintAges,self.UseStaggeredAges)

        if self.UseDiffDemog:
            self.A = np.ones(self.I)+np.cumsum(np.ones(self.I)*.05)-.05 #Techonological Change, used for when countries are different
            #A = np.ones(I)
        else:
            self.A = self.np.ones(I) #Techonological Change, used for idential countries

        if self.UseDiffProductivities:
            self.e = np.ones((self.I, self.S, self.T))
            self.e[:,self.FirstDyingAge:,:] = 0.3
            self.e[:,:self.LeaveHouseAge,:] = 0.3
        else:
            self.e = np.ones((self.I, self.S, self.T)) #Labor productivities

    #DEMOGRAPHICS SET-UP

    def Import_Data(self):

        self.f_range=self.LastFertilityAge+1-self.FirstFertilityAge


        self.N=np.zeros((self.I,self.S,self.T))
        self.Nhat=np.zeros((self.I,self.S,self.T))
        self.all_FertilityRates = np.zeros((self.I, self.S, self.f_range+self.T))
        self.FertilityRates = np.zeros((self.I, self.S, self.T))
        self.MortalityRates = np.zeros((self.I, self.S, self.T))
        self.Migrants = np.zeros((self.I, self.S, self.T))
        self.g_N = np.zeros(self.T)
        self.lbar = np.zeros(self.T)

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
                        :self.f_range+self.T_1] =  np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_fertility.csv"),delimiter=',',skiprows=1\
                        , usecols=(self.agestopull[self.FirstFertilityAge:self.LastFertilityAge+1]-22))[48-self.f_range:48+self.T_1,:])

                self.MortalityRates[i,self.FirstDyingAge:,:self.T_1] = np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_mortality.csv"),delimiter=','\
                        ,skiprows=1, usecols=(self.agestopull[self.FirstDyingAge:]-67))[:self.T_1,:])

                self.Migrants[i,:self.MaxImmigrantAge,:self.T_1] = np.einsum("s,t->st",np.loadtxt(("Data_Files/net_migration.csv"),delimiter=','\
                        ,skiprows=1, usecols=[index+1])[self.agestopull[:self.MaxImmigrantAge]]*100, np.ones(self.T_1))

            else:
                self.N[i,:,0] = np.loadtxt(("Data_Files/population.csv"),delimiter=',',skiprows=1, usecols=[index+1])[:self.S]*1000

                self.all_FertilityRates[i,self.FirstFertilityAge:self.LastFertilityAge+1,:self.f_range+self.T_1] = \
                        np.transpose(np.loadtxt(str("Data_Files/" + I_all[index] + "_fertility.csv"),delimiter=','\
                        ,skiprows=1, usecols=range(1,self.f_range+1))[48-self.f_range:48+self.T_1,:])

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
        self.f_bar = np.mean(self.all_FertilityRates[:,:,self.f_range+self.T_1-1], axis=0)
        self.rho_bar = np.mean(self.MortalityRates[:,:,self.T_1-1], axis=0)

        #Set to the steady state for every year beyond year T_1
        self.all_FertilityRates[:,:,self.f_range+self.T_1:] = np.tile(np.expand_dims(self.f_bar, axis=2), (self.I,1,self.T-self.T_1))
        self.MortalityRates[:,:,self.T_1:] = np.tile(np.expand_dims(self.rho_bar, axis=2), (self.I,1,self.T-self.T_1))

        #FertilityRates is exactly like all_FertilityRates except it begins at time t=0 rather than time t=-f_range
        self.FertilityRates[:,self.FirstFertilityAge:self.LastFertilityAge+1,:] = self.all_FertilityRates[:,self.FirstFertilityAge:self.LastFertilityAge+1,self.f_range:]

        #Gets initial world population growth rate
        self.g_N[0] = 0.

    def Demographics(self):

        self.ImmigrationRates = np.zeros((self.I,self.S,self.T))

        N_temp = np.ones((self.I,self.S))/(self.I*self.S)

        for t in xrange(1,self.T):

            self.N[:,0,t] = np.sum((self.N[:,:,t-1]*self.FertilityRates[:,:,t-1]), axis=1)
            N_temp[:,0] = np.sum((self.Nhat[:,:,t-1]*self.FertilityRates[:,:,t-1]), axis=1)


            if t <= self.T_1:
                self.ImmigrationRates[:,:,t-1] = self.Migrants[:,:,t-1]/self.N[:,:,t-1]

            else:
                self.ImmigrationRates[:,:,t-1] = np.mean(self.ImmigrationRates[:,:,self.T_1],\
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

        itera = 0

        while np.max(np.abs(self.Nhat[:,:,-1] - self.Nhat[:,:,-2])) > self.demog_ss_tol:
            pop_new[:,0] = np.sum((pop_old[:,:]*self.FertilityRates[:,:,-1]),axis=1)
            pop_new[:,1:] = pop_old[:,:-1]*(1+self.ImmigrationRates[:,:-1,-1]\
                    -self.MortalityRates[:,:-1,-1])
            self.Nhat = np.dstack((self.Nhat,pop_new/np.sum(pop_new)))
            itera += 1

        if self.PrintLoc: print "The SS Population Share converged in", iter, "years beyond year T"

        self.Nhat_ss = self.Nhat[:,:,-1] 

        if self.DemogGraphs:
            ages = self.FirstFertilityAge, self.LastFertilityAge, self.FirstDyingAge, \
                    self.MaxImmigrantAge
            datasets = self.FertilityRates, self.MortalityRates, self.ImmigrationRates, self.Nhat
            demog.plotDemographics(ages, datasets, self.I, self.S, self.T, self.I_touse, T_touse = [0,1,2,3,20]\
                    , compare_across="T", data_year=0)

        self.lbar[:self.T] = np.cumsum(np.ones(self.T)*self.g_A)

        self.lbar[self.T-self.S:] = np.ones(self.S)

        if self.CheckerMode==False:
            print "\nDemographics obtained!"

        self.Nhat=self.Nhat[:,:,self.T]

        self.e_ss=self.e[:,:,-1]

        self.Mortality_ss=self.MortalityRates[:,:,-1]


        self.lbar_ss=self.lbar[-1]






    #STEADY STATE

    def get_Psi(self, w, e):

        if e.ndim == 2:
            denom = np.einsum("i,is->is",w,e)
            denom = (self.chi/denom)**self.rho

        elif e.ndim == 3:
            denom = np.einsum("it, ist -> ist", w, e)
            denom = (self.chi/denom)**self.rho

        psi = (1+self.chi*denom)**((1-self.rho*(self.sigma))/self.rho)

        return psi

    def get_lhat(self,c,w,e):

        if e.ndim == 2:
            lhat=c*(self.chi/(np.einsum("i,is->is",w,e)))**self.rho
        elif e.ndim == 3:
            lhat=c*(self.chi/(np.einsum("it,ist->ist",w,e)))**self.rho

        return lhat

    def get_n(self, lhat):

        if self.Nhat.shape == 3:
            self.I, self.S, self.T = self.Nhat.shape
            lhat = np.einsum("t,is->ist",lhat,np.ones((self.I,self.S)))

        #print self.e_ss.shape
        #print self.lbar.shape
        #print lhat.shape
        #print self.Nhat.shape

        n = np.sum(self.e_ss*(self.lbar_ss-lhat)*self.Nhat,axis=1)

        return n

    def get_Y(self, kd, n):

        if kd.ndim ==1:
            print kd
            Y = (kd**self.alpha) * ((self.A*n)**(1-self.alpha))
        elif kd.ndim== 2:
            Y = (kd**self.alpha) * (np.einsum("i,is->is",self.A,n)**(1-self.alpha))

        return Y


    def GetSSComponents(self, bq_ss, r_ss):

        def get_lifetime_decisionsSS(c_1, w_ss, r_ss):

            #print self.I

            cvec_ss = np.zeros((self.I,self.S))
            avec_ss = np.zeros((self.I,self.S+1))
            cvec_ss[:,0] = c_1

            for s in xrange(self.S-1):
                cvec_ss[:,s+1] = (self.beta * (1-self.Mortality_ss[:,s]) * (1 + r_ss - self.delta)\
                        *self.psi_ss[:,s+1]/self.psi_ss[:,s])**(1/self.sigma) * cvec_ss[:,s]*np.exp(-self.g_A)

                avec_ss[:,s+1] = (w_ss*self.e_ss[:,s] + (1 + r_ss - self.delta)*avec_ss[:,s] + \
                        bq_ss[:,s] - cvec_ss[:,s]*(1+w_ss*self.e_ss[:,s]*\
                        (self.chi/(w_ss*self.e_ss[:,s])**self.rho)))*np.exp(-self.g_A)

            avec_ss[:,s+2] = (w_ss*self.e_ss[:,s+1] + (1 + r_ss - self.delta)*avec_ss[:,s+1] \
                    - cvec_ss[:,s+1]*(1+w_ss*self.e_ss[:,s+1]*(self.chi/(w_ss*self.e_ss[:,s+1])\
                    **self.rho)))*np.exp(-self.g_A)
            return cvec_ss, avec_ss

        def householdEuler_SS(c_1, w_ss, r_ss):

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

        return w_ss, cvec_ss, avec_ss, kd_ss, kf_ss, n_ss, y_ss


    def EulerSystemSS(self,guess):

        bq_ss = guess[:-1]
        r_ss = guess[-1]

        bqvec_ss = np.zeros((self.I,self.S))
        bqvec_ss[:,self.FirstFertilityAge:self.FirstDyingAge] = np.einsum("i,s->is", bq_ss, \
                np.ones(self.FirstDyingAge-self.FirstFertilityAge))

        w_ss, cvec_ss, avec_ss, kd_ss, kf_ss, n_ss, y_ss = self.GetSSComponents(bqvec_ss, r_ss)


        alldeadagent_assets = np.sum(avec_ss[:,self.FirstDyingAge:]*\
                self.Mortality_ss[:,self.FirstDyingAge:]*self.Nhat_ss[:,self.FirstDyingAge:], axis=1)

        Euler_bq = bq_ss - alldeadagent_assets/np.sum(self.Nhat_ss[:,self.FirstFertilityAge:self.FirstDyingAge],\
                axis=1)

        Euler_kf = np.sum(kf_ss)


        Euler_all = np.append(Euler_bq, Euler_kf)


        if self.EulErrors: print "Euler Errors:", Euler_all

        return Euler_all


    def SteadyState(self, rss_guess, bqss_guess):

        guess = np.append(bqss_guess, rss_guess)

        ss = opt.fsolve(self.EulerSystemSS, guess)

        self.bq_ss = ss[:-1]

        self.r_ss = ss[-1]

        self.bqvec_ss = np.zeros((self.I,self.S))
        self.bqvec_ss[:,self.FirstFertilityAge:self.FirstDyingAge] = np.einsum("i,s->is",self.bq_ss,\
                np.ones(self.FirstDyingAge-self.FirstFertilityAge))

        self.w_ss, self.cvec_ss, self.avec_ss, self.kd_ss, self.kf_ss, self.n_ss, self.y_ss \
                = self.GetSSComponents(self.bqvec_ss,self.r_ss)


        print "\n\nSTEADY STATE FOUND!"
        
        alldeadagent_assets = np.sum(self.avec_ss[:,self.FirstDyingAge:]*self.Mortality_ss[:,self.FirstDyingAge:]*\
                self.Nhat_ss[:,self.FirstDyingAge:], axis=1)

        print "Euler_bq", self.bq_ss - alldeadagent_assets/np.sum(self.Nhat_ss[:,self.FirstFertilityAge:self.FirstDyingAge],axis=1)

        print np.isclose(alldeadagent_assets/np.sum(self.Nhat_ss[:,self.FirstFertilityAge:self.FirstDyingAge],axis=1),0)



        if self.PrintSS:
            print "assets steady state:", self.avec_ss
            print "kf steady state", self.kf_ss
            print "kd steady state", self.kd_ss
            print "n steady state", self.n_ss
            print "y steady state", self.y_ss
            print "r steady state", self.r_ss
            print "w steady state", self.w_ss
            print "c_vec_ss steady state", self.cvec_ss








    #TIMEPATH-ITERATION




def Multi_Country(S,I,sigma):

    #THIS SETS ALL OF THE USER PARAMETERS
    I_dict = {"usa":0,"eu":1,"japan":2,"china":3,"india":4,"russia":5,"korea":6}
    #Parameters Zone
    #T = int(round(4*S)) #Number of time periods to convergence, based on Rick Evans' function.
    I_touse = ["eu","russia","usa","japan","korea","china","india"] 

    g_A = 0.015 #Technical growth rate

    beta_ann=.95 #Starting future consumption discount rate
    delta_ann=.08 #Starting depreciation rate
    alpha = .3 #Capital Share of production
    chi = 1.5 #New Parameter
    rho = 1.3 #Other New Parameter


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
    CheckerMode = False #Reduces the number of prints when checking for robustness

    DemogGraphs = False #Activates graphing graphs with demographic data and population shares
    TPIGraphs = False #Activates graphing the graphs.

    UseStaggeredAges = True #Activates using staggered ages
    UseDiffDemog = True #Turns on different demographics for each country
    UseSSDemog = False #Activates using only steady state demographics for TPI calculation
    UseDiffProductivities = True #Activates having e vary across cohorts
    UseTape = True #Activates setting any value of kd<0 to 0.001 in TPI calculation
    SAVE = False #Saves the graphs
    SHOW = True #Shows the graphs
    ADJUSTKOREAIMMIGRATION = True

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

    HH_params = (S,I,beta_ann,sigma)#,nvec)

    Firm_Params = ( alpha, delta_ann, chi, rho, g_A)

    Tolerances = (tpi_tol, demog_ss_tol)

    Levers = (CalcTPI, PrintAges,PrintLoc,PrintEulErrors,PrintSS,Print_cabqTimepaths,CheckerMode,DemogGraphs,TPIGraphs,\
            UseStaggeredAges,UseDiffDemog,UseSSDemog,UseDiffProductivities,UseTape,SAVE,SHOW,ADJUSTKOREAIMMIGRATION)

    TPI_Params = ()
    SS_Params = ()



    ##WHERE THE MAGIC HAPPENS ##

    Model = OLG(Country_Roster,HH_params,Firm_Params,Levers, Tolerances, SS_Params, TPI_Params)

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


Multi_Country(20,2,2)

