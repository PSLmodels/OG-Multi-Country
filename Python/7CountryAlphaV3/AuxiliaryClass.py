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
    accessed by any of the functions. Inputting all of the parameters starts at about line 585.
    
    Unlike the previous versions of the code, the comments for this version will indicate inputs, outputs, and which
    objects will be stored in the objects.

    KEY:

    * - Variable that isn't stored in object
    # - Variable that is stored in object 

    """

    def __init__(self, countries, HH_Params, Firm_Params, Lever_Params, Tol_Params, TPI_Params):
        """
        Description: 
            -This creates the object and stores all of the parameters into the object.
             The initialization is the starting point for model.

        Inputs:
            -self: "self" stores all of the components of the model. To access any part,
             simply type "self.variable_name" while in the object and "objectname.variable_name"
             outside the object.

            -countries = tuple: contains a dictionary and tuple for countries and their associated number

            -HH_Params = tuple: contains S, I, annualized Beta and sigma.

            -Firm_Params = tuple: contains alpha, annualized delta, chi, rho and g_A

            -Lever_Params = tuple: Very large tuple that contains all of the binary variables given by the user in the input section.

            -Tol_Params = tuple: contains the entered tolerance levels for TPI and demographics

            -TPI_Params = tuple: contains the xi parameter and the maximum number of iterations.


        
        Variables Stored in Object:
            - self.T = Scalar: of the total amount of time periods
            - self.T_1 = Scalar: Transition year for the demographics
            - self.delta = Scalar: calulated overall depreciation rate
            - self.beta = Scalar: calculated overall future discount rate
            - self.LeaveHouseAge = Scalar:  From the Auxiliary Demographics module, see that page for details
            - self.FirstFertilityAge = Scalar: From the Auxiliary Demographics module, see that page for details
            - self.LastFertilityAge = Scalar: From the Auxiliary Demographics module, see that page for details
            - self.MaxImmigrantAge = Scalar: From the Auxiliary Demographics module, see that page for details
            - self.FirstDyingAge = Scalar: From the Auxiliary Demographics module, see that page for details
            - self.agestopull = Scalar: From the Auxiliary Demographics module, see that page for details
            - self.A = Vector [I,1], Technology level for each country
            - self.e = Matrix [I,S,T], Labor Productivities
            - All of these Objects in Function, as well as the contents of the Tuples were saved in the object



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
        (self.CalcTPI,self.PrintAges,self.PrintLoc,self.EulErrors,self.PrintSS,self.Print_cabqTimepaths,self.CheckerMode,\
                self.DemogGraphs,self.TPIGraphs,self.UseStaggeredAges,self.UseDiffDemog, self.UseSSDemog,\
                self.UseDiffProductivities,self.UseTape,self.SAVE,self.SHOW,self.ADJUSTKOREAIMMIGRATION) = Lever_Params

        #Tolerance Parameters

        (self.tpi_tol, self.demog_ss_tol) = Tol_Params


        #Timepath Iteration Parameters
        (self.xi, self.MaxIters) = TPI_Params

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


        self.rpathlist = np.empty((1,self.T+self.S))


    #DEMOGRAPHICS SET-UP

    def Import_Data(self):
        """
        Description: Imports the data files (.csv). Additionally, it creates and stores the
        data in the object


        Variables used in Object: 
        -LastFertility
        -FirstFertility
        -I
        -S
        -T
        -I_all
        -
        
        Variables stored in object:
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
        Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

        Inputs:

        
        Objects in Function:


        Outputs:

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

        #print np.sum(np.sum(self.Nhat,axis=0),axis=0)
        if self.DemogGraphs:
            ages = self.FirstFertilityAge, self.LastFertilityAge, self.FirstDyingAge, \
                    self.MaxImmigrantAge
            datasets = self.FertilityRates, self.MortalityRates, self.ImmigrationRates, self.Nhat
            demog.plotDemographics(ages, datasets, self.I, self.S, self.T, self.I_touse, T_touse = [0,1,2,3,20]\
                    , compare_across="T", data_year=0)


    #STEADY STATE

    def get_Psi(self, w, e):

        """
        Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

        Inputs:

        
        Objects in Function:


        Outputs:

        """


        if e.ndim == 2:
            we =  np.einsum("i,is->is",w,e)

        elif e.ndim == 3:
            we = np.einsum("it, ist -> ist", w, e)

        part1 = (self.chi/we)**(self.rho-1)

        psi = (1+self.chi*part1)**((1-self.rho*(self.sigma))/(self.rho-1))

        return psi

    def get_lhat(self,c,w,e):

        """
        Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

        Inputs:

        
        Objects in Function:


        Outputs:

        """


        if e.ndim == 2:
            lhat=c*(self.chi/(np.einsum("i,is->is",w,e)))**self.rho
        elif e.ndim == 3:
            lhat=c*(self.chi/(np.einsum("it,ist->ist",w,e)))**self.rho

        return lhat

    def get_n(self, lhat):
        """
        Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

        Inputs:

        
        Objects in Function:


        Outputs:

        """

        if lhat.ndim == 2:
            n = np.sum(self.e_ss*(self.lbar_ss-lhat)*self.Nhat_ss,axis=1)
        elif lhat.ndim == 3:
            n = np.sum(self.e[:,:,:self.T]*(self.lbar[:self.T]-lhat)*self.Nhat[:,:,:self.T],axis=1)

        return n

    def get_Y(self, kd, n):
        """
        Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

        Inputs:

        
        Objects in Function:


        Outputs:

        """

        if kd.ndim ==1:
            Y = (kd**self.alpha) * ((self.A*n)**(1-self.alpha))
        elif kd.ndim== 2:
            Y = (kd**self.alpha) * (np.einsum("i,is->is",self.A,n)**(1-self.alpha))

        return Y

    def GetSSComponents(self, bq_ss, r_ss):
        """
        Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

        Inputs:

        
        Objects in Function:


        Outputs:

        """


        def get_lifetime_decisionsSS(c_1, w_ss, r_ss):
            """
            Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

            Inputs:

        
            Objects in Function:


            Outputs:

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
            Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

            Inputs:

        
            Objects in Function:


            Outputs:

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
        Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

        Inputs:

        
        Objects in Function:


        Outputs:

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

        if self.EulErrors: print "Euler Errors:", Euler_all

        return Euler_all

    def SteadyState(self, rss_guess, bqss_guess):
        """
        Description: Initializes the storage of all of the parameters into the objects. By 
                    typing "self.variable_name" in the object, you can access any of the inputs.

        Inputs:

        
        Objects in Function:


        Outputs:

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
        we = np.einsum("i,is->is",self.w_ss,self.e_ss[:,:-1])

        print self.psi_ss[:,:-1]*self.cvec_ss[:,:-1]**(-self.sigma) - self.beta*(1-self.Mortality_ss[:,:-1])*self.psi_ss[:,1:]*(self.cvec_ss[:,1:]*np.exp(self.g_A))**(-self.sigma)*(1+self.r_ss-self.delta)
        
        print self.cvec_ss[:,:-1] - \
        (we + (1+self.r_ss-self.delta)*self.avec_ss[:,:-1] + self.bqvec_ss[:,:-1] - self.avec_ss[:,1:]*np.exp(self.g_A)) / \
        (1 + we*(self.chi/we)**self.rho)


    #TIMEPATH-ITERATION

    def set_initial_values(self, r_init, bq_init, a_init):
        self.r_init = r_init
        self.bq_init = bq_init
        self.a_init = a_init

    def get_initialguesses(self):

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

        def get_lifetime_decisionsTPI(c_1, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age,lastguy =True):
            
            decisions = self.S - age -1
            cvec_path = np.zeros((self.I,decisions+1))
            avec_path = np.zeros((self.I,decisions+2))
            cvec_path[:,0] = c_1
            avec_path[:,0] = a_current

            for s in range(decisions):

                cvec_path[:,s+1] = (self.beta * (1-mort_life[:,s]) * (1 + r_path[s+1] - self.delta)\
                                   * psi_life[:,s+1]/psi_life[:,s])**(1/self.sigma) * cvec_path[:,s]*np.exp(-self.g_A)

                avec_path[:,s+1] = (w_life[:,s]*e_life[:,s] + (1 + r_life[s] - self.delta)*avec_path[:,s] + \
                        bq_life[:,s] - cvec_path[:,s]*(1+w_life[:,s]*e_life[:,s]*\
                        (self.chi/(w_life[:,s]*e_life[:,s]))**self.rho))*np.exp(-self.g_A)

                avec_path[:,s+2] = (w_life[:,s+1]*e_life[:,s+1] + (1 + r_life[s+1] - self.delta)*avec_path[:,s+1] \
                        - cvec_path[:,s+1]*(1+w_life[:,s+1]*e_life[:,s+1]*(self.chi/(w_path[:,s+1]*e_life[:,s+1]))\
                        **self.rho))*np.exp(-self.g_A)

            return cvec_path, avec_path

        def optc1_Euler_TPI(c1_guess, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age,lastguy = True):

            cpath_indiv, apath_indiv = get_lifetime_decisionsTPI(c1_guess, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age)
            Euler = np.ravel(apath_indiv[:,-1])

            if np.any(cpath_indiv<0):
                print "WARNING! The fsolve for initial optimal consumption guessed a negative number"
                Euler = np.ones(Euler.shape[0])*9999.

            return Euler

        def get_c_a_matrices(w_path, r_path, psi, bqvec_path):

            c_matrix = np.zeros((self.I,self.S,self.T+self.S))
            a_matrix = np.zeros((self.I,self.S+1,self.T+self.S))
            a_matrix[:,:-1,0] = self.a_init
            
            c_matrix[:,self.S-1,0] = (w_path[:,0]*self.e[:,self.S-1,0] + (1 + r_path[0] - self.delta)*self.a_init[:,self.S-1])\
            /(1+w_path[:,0]*self.e[:,self.S-1,0]*(self.chi/(w_path[:,0]*self.e[:,self.S-1,0]))**self.rho)
            #print self.a_init[:,self.S-1], "HERE"
            """
            print c_matrix[0,self.S-1,0]
            test = (self.w_ss*self.e_ss[:,0] + (1 + self.r_ss - self.delta)*self.avec_ss[:,self.S-1])\
            /(1+self.w_ss*self.e_ss[:,0]*(self.chi/(self.w_ss*self.e_ss[:,0]))**self.rho)
            print test[0]
            print self.cvec_ss[0,-1]
            """
            #c1_guess = np.ones(I)*.212
            #yeaaaaas = opt.fsolve(optc1_Euler_TPI, c1_guess, args = (w_path[:,0], r_path[0], self.MortalityRates[:,self.S-1,0], self.e[:,self.S-1,0], 0, self.avec_ss[:,self.S-1], self.S-1))            


            for age in range(self.S-2,0,-1):
                
                c1_guess = (c_matrix[:,age+1,0]*(psi[:,age,0]/psi[:,age+1,1])\
                    /((self.beta*(1+r_path[0]-self.delta))**(1/self.sigma)))/np.exp(self.g_A)
                
                p = self.S-age #Remaining decisions

                w_life = w_path[:,:p]
                r_life = r_path[:p+1]
                mort_life = np.diagonal(self.MortalityRates[:,age:,age:], axis1=1, axis2=2)
                e_life = np.diagonal(self.e[:,age:,:p+1], axis1=1, axis2=2)
                psi_life = np.diagonal(psi[:,age:,:p+2], axis1=1, axis2=2)
                bq_life = np.diagonal(bqvec_path[:,age:,:p+1], axis1=1, axis2=2)
                a_current = self.a_init[:,age]

                opt_c1 = opt.fsolve(optc1_Euler_TPI, c1_guess, args = (w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age))
            
                cpath_indiv, apath_indiv = get_lifetime_decisionsTPI(opt_c1, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age)
                
                for i in xrange(self.I):
                    np.fill_diagonal(c_matrix[i,age:,:], cpath_indiv[i,:])
                    np.fill_diagonal(a_matrix[i,age:,:], apath_indiv[i,:])

            for t in range(self.T):

                c1_guess = c_matrix[:,0,t-1]

                age = 0
                w_life = w_path[:,t:t+self.S+1]
                r_life = r_path[t:t+self.S+1]
                mort_life = np.diagonal(self.MortalityRates[:,:,t:t+self.S+1], axis1=1, axis2=2)
                e_life = np.diagonal(self.e[:,:,t:t+self.S+1], axis1=1, axis2=2)
                psi_life = np.diagonal(psi[:,:,t:t+self.S+1], axis1=1, axis2=2)
                bq_life = np.diagonal(bqvec_path[:,:,t:t+self.S+1], axis1=1, axis2=2)
                a_current = np.zeros(self.I)

                opt_c1 = opt.fsolve(optc1_Euler_TPI, c1_guess, args = (w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age))
            
                cpath_indiv, apath_indiv = get_lifetime_decisionsTPI(opt_c1, w_life, r_life, mort_life, e_life, psi_life, bq_life, a_current, age)
                
                for i in xrange(self.I):
                    np.fill_diagonal(c_matrix[i,:,t:], cpath_indiv[i,:])
                    np.fill_diagonal(a_matrix[i,:,t:], apath_indiv[i,:])

                if self.Print_cabqTimepaths:
                    print "Consumption for year", t
                    print np.round(np.transpose(c_matrix[0,:,:self.T]), decimals=3)
                    #print "c_guess", np.round(c_guess[0], decimals=3)
                    print "Assets for year", t
                    print np.round(np.transpose(a_matrix[0,:,:self.T]), decimals=3)
                    #print "Bequests"
                    #print np.round(np.transpose(bq_timepath[0,:,:p+2]), decimals=3)
                    #print "agent_bq", np.round(agent_bq[0,:], decimals=3)

            print "Euler Household satisfied:", np.isclose(np.max(np.absolute(a_matrix[:,-1,:])), 0)

            return c_matrix[:,:,:self.T], a_matrix[:,:-1,:self.T]

        w_path = np.einsum("it,i->it",np.einsum("i,t->it",self.alpha*self.A,1/r_path)**(self.alpha/(1-self.alpha)),(1-self.alpha)*self.A)

        psi = self.get_Psi(w_path,self.e)

        c_matrix, a_matrix = get_c_a_matrices(w_path, r_path, psi, bqvec_path)

        lhat_path = self.get_lhat(c_matrix, w_path[:,:self.T], self.e[:,:,:self.T])

        n_path = self.get_n(lhat_path)

        kd_path = np.sum(a_matrix*self.Nhat[:,:,:self.T],axis=1)

        y_path = self.get_Y(kd_path,n_path)

        kf_path = np.outer(self.alpha*self.A, 1/r_path[:self.T])**(1/(1-self.alpha)) * n_path-kd_path

        return w_path, c_matrix, a_matrix, kd_path, kf_path, n_path, y_path, lhat_path

    def EulerSystemTPI(self, guess):

        guess = np.expand_dims(guess, axis=1).reshape((self.I+1,self.T))
        r_path = guess[0,:]
        r_path = np.hstack((r_path, np.ones(self.S)*self.r_ss))
        bq_path = guess[1:,:]
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

        Euler_all = np.append(Euler_bq, Euler_kf)

        if self.EulErrors: 
            print "Iteration:", self.Timepath_counter, "Min Euler:", np.min(np.absolute(Euler_all)), "Mean Euler:", np.mean(np.absolute(Euler_all)), "Max Euler_bq:", np.max(np.absolute(Euler_bq)), "Max Euler_kf", np.max(np.absolute(Euler_kf))

        iterations_to_plot = set([1,600])

        if self.Timepath_counter in iterations_to_plot:
            self.plot_iteration(r_path, bq_path, w_path, c_matrix, lhat_path, n_path, kd_path, kf_path)

        self.rpathlist = np.vstack((self.rpathlist, r_path))
        
        self.Timepath_counter += 1
        
        return Euler_all

    def Timepath(self):
        
        rpath_guess, bqpath_guess = self.get_initialguesses()

        print rpath_guess.shape, bqpath_guess.shape
        guess = np.append(rpath_guess, bqpath_guess)

        paths = opt.fsolve(self.EulerSystemTPI, guess, xtol=1e-4)

        print rpath_guess.shape, bqpath_guess.shape
        print paths.shape

        w_path, c_matrix, a_matrix, kd_path, kf_path, n_path, y_path, lhat_path = self.GetTPIComponents(bqvec_path, r_path)

        self.plot_iteration(r_path, bq_path, w_path, c_matrix, lhat_path, n_path, kd_path, kf_path)

    def plot_iteration(self, r_path, bq_path, w_path, c_matrix, lhat_path, n_path, kd_path, kf_path):

        title = str("S = " + str(self.S) + ", T = " + str(self.T))
        plt.suptitle(title)

        plt.subplot(331)
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


        plt.show()

    def GETMAGICINDICESTOUSELATER(self):

        i_indices = []
        j_indices = []

        I=7
        S=8
        T=S*4

        z = np.zeros((I,S,S+T))
        a = np.array(range(1,S+1)).T*1.
        for r in range(1,T):
            a = np.vstack( (a, np.linspace(1+r*S,S+S*r,S).T) )
        a = a.T
        a = np.tile(a, (I,1,1))

        s_indices = np.repeat(range(S), T)
        t_indices = sum([range(s,T+s) for s in range(S)], [])

        z[:,s_indices,t_indices] = a.reshape((I,S*T))
        b = z[:,s_indices,t_indices].reshape(I,S,T)

