from __future__ import division
import sys
import numpy as np
import scipy as sp
import scipy.optimize as opt
import time as time
from matplotlib import pyplot as plt

#DEMOGRAPHICS FUNCTIONS

def getDemographics(params, PrintAges, DiffDemog, countrynames):
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
		PrintAges         = Boolean, prints LeaveHouseAge, FirstFertilityAge, LastFertilityAge, FirstDyingAge, and MaxImmigrantAge if==True
		DiffDemog         = Boolean, Causes data to load from different csv files for each country if==True
		countrynames      = Tuple of length (I), contains the names of each of the I countries, used for pulling data from csv files

	Functions called:
		-None

	Objects in Function:
		f_range            = Scalar >0, Number of years the agents are fertile
		N_matrix           = [I, S+1, T+S+1] Matrix, Population for each country, generation, and year
		Nhat_matrix        = [I, S+1, T+S+1] Matrix, Population share for each country, generation, and year
		KIDs               = [I, S+1, T+S+1] Matrix, Number of kids for each country, generation, and year
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
		KIDs_ss            = [I, S+1] Matrix, Number of kids in the steady-state for each country and generation
		Nhatss_old         = [I, S+1] Matrix, Compared to Nhatss_new for converging to the steady-state
		Nhatss_new         = [I, S+1] Matrix, Compared to Nhatss_old for converging to the steady-state
		pop_old            = [I, S+1] Matrix, Compared to pop_new for converging to the steady-state
		pop_new            = [I, S+1] Matrix, Compared to pop_old for converging to the steady-state
		iter               = Scalar >=0, Counts the number of years beyond T+S+1 until Nhat_matrix converges to the steady-state
		lbar               = [T+S+1,] Vector, Labor endowment for each year

	Returns: MortalityRates, Nhat_matrix, KIDs, Nhatss_new, KIDs_ss, lbar
	"""
	#Unpacks parameters
	I, S, T, T_1, LeaveHouseAge, FirstFertilityAge, LastFertilityAge, FirstDyingAge, MaxImmigrantAge, g_A, tol = params

	#Parameter that represents total number of fertile years. Used for indexing purposes and to make the code easier to read
	f_range = LastFertilityAge+1-FirstFertilityAge

	#Initializes demographics matrices
	N_matrix = np.zeros((I, S+1, T+S+1))
	Nhat_matrix = np.zeros((I, S+1, T+S+1))
	KIDs = np.zeros((I, S+1, T+S+1))
	all_FertilityRates = np.zeros((I, S+1, f_range+T+S+1))
	FertilityRates = np.zeros((I, S+1, T+S+1))
	MortalityRates = np.zeros((I, S+1, T+S+1))
	ImmigrationRates = np.zeros((I, S+1, T+S+1))
	Migrants = np.zeros((I, S+1, T+S+1))
	g_N = np.zeros(T+S+1)
	#QUESTION HOW TO INTERPRET N_temp???
	N_temp = np.ones((I, S+1))/(I*S)

	if PrintAges:
		print "T =", T
		print "T_1", T_1
		print "LeaveHouseAge", LeaveHouseAge
		print "FirstFertilityAge", FirstFertilityAge
		print "LastFertilityAge", LastFertilityAge
		print "FirstDyingAge", FirstDyingAge
		print "MaxImmigrantAge", MaxImmigrantAge

	#Gathers demographic data from the csv files for each I countries in countrynames
	for i in range(I):

		#If we want the countries to have unique demographic data, have the same index
		if DiffDemog:
			index = i
			if I > 7:
				sys.exit("ERROR!!! We can't have more than 7 Countries without unique data. Change either parameter I so it is less than 8 or change DiffDemog to False")
		#If we want the countries to have identical demographic data, have the same index
		else:
			index = 0

		N_matrix[i,:,0] = np.loadtxt(("Data_Files/population.csv"),delimiter=',',skiprows=1, usecols=[index+1])[:S+1]*1000
		all_FertilityRates[i,FirstFertilityAge:LastFertilityAge+1,:f_range+T_1] = np.transpose(np.loadtxt(str("Data_Files/" + countrynames[index] + "_fertility.csv"),delimiter=',',skiprows=1, usecols=range(1,f_range+1))[48-f_range:48+T_1,:])
		MortalityRates[i,FirstDyingAge:-1,:T_1] = np.transpose(np.loadtxt(str("Data_Files/" + countrynames[index] + "_mortality.csv"),delimiter=',',skiprows=1, usecols=range(1,S+1-FirstDyingAge))[:T_1,:])
		Migrants[i,:MaxImmigrantAge,:T_1] = np.einsum("s,t->st",np.loadtxt(("Data_Files/net_migration.csv"),delimiter=',',skiprows=1, usecols=[index+1])[:MaxImmigrantAge]*100, np.ones(T_1))
		print "Got demographics for", countrynames[index]

	#Gets initial population share
	Nhat_matrix[:,:,0] = N_matrix[:,:,0]/np.sum(N_matrix[:,:,0])

	#The last generation dies with probability 1
	MortalityRates[:,-1,:] = np.ones((I, T+S+1))

	#Gets steady-state values for all countries by taking the mean at year T_1-1 across countries
	f_bar = np.mean(all_FertilityRates[:,:,f_range+T_1-1], axis=0)
	rho_bar = np.mean(MortalityRates[:,:,T_1-1], axis=0)

	#Set to the steady state for every year beyond year T_1
	all_FertilityRates[:,:,f_range+T_1:] = np.tile(np.expand_dims(f_bar, axis=2), (I,1,T-T_1+S+1))
	MortalityRates[:,:,T_1:] = np.tile(np.expand_dims(rho_bar, axis=2), (I,1,T-T_1+S+1))

	#FertilityRates is exactly like all_FertilityRates except it begins at time t=0 rather than time t=-f_range
	FertilityRates[:,FirstFertilityAge:LastFertilityAge+1,:] = all_FertilityRates[:,FirstFertilityAge:LastFertilityAge+1,f_range:]

	#Gets initial world population growth rate
	g_N[0] = 0.

	#Calculates population numbers for each country
	for t in range(1,T+S+1):
		#Gets the total number of children and and percentage of children and stores them in generation 0 of their respective matrices
		#See equation 2.1
		N_matrix[:,0,t] = np.sum((N_matrix[:,:,t-1]*FertilityRates[:,:,t-1]), axis=1)
		N_temp[:,0] = np.sum((Nhat_matrix[:,:,t-1]*FertilityRates[:,:,t-1]), axis=1)

		#Finds the immigration rate for each year
		if t <= T_1:
			ImmigrationRates[:,:,t-1] = Migrants[:,:,t-1]/N_matrix[:,:,t-1]
		else:
			ImmigrationRates[:,:,t-1] = np.mean(ImmigrationRates[:,:,T_1-1], axis=0)
			
		#Gets the population distribution for the next year, taking into account immigration and mortality
		#See equation 2.2
		N_matrix[:,1:,t] = N_matrix[:,:-1,t-1]*(1+ImmigrationRates[:,:-1,t-1]-MortalityRates[:,:-1,t-1])
		N_temp[:,1:] = Nhat_matrix[:,:-1,t-1]*(1+ImmigrationRates[:,:-1,t-1]-MortalityRates[:,:-1,t-1])

		#Gets the population share
		Nhat_matrix[:,:,t] = N_matrix[:,:,t]/np.sum(N_matrix[:,:,t])

		#Gets the total number of kids living at home for each generation by looping through each generation that will have kids at home
		for f in range(FirstFertilityAge+1, LastFertilityAge+1+LeaveHouseAge):
			#Gets diagonal of fertility rates and sums them up to get the number of kids. See equation 3.16
			KIDs[:,f,t-1] = np.sum(np.diagonal(all_FertilityRates[:,f-LeaveHouseAge:f,t:], axis1=1, axis2=2), axis=1)
			"""
			print t, f
			print np.round(np.transpose(all_FertilityRates[0,:,0:f_range]), decimals=4)
			print "BREAK", str(str(f_range) + "+" + str(f) + "-" + str(FirstFertilityAge) + "-" + str(LeaveHouseAge))
			print np.round(np.transpose(all_FertilityRates[0,:,f_range:f_range+t+2]), decimals=4)
			print "BREAK"
			print np.round(np.transpose(KIDs[0,:,t-1]), decimals=4)
			print range(f-LeaveHouseAge, f), np.round(np.transpose(agent_kids[0,:]), decimals=4), np.sum(np.round(np.transpose(agent_kids[0,:]), decimals=4))
			"""
		#Gets the growth rate for the next year
		g_N[t] = np.sum(N_temp[:,:])-1

	#Fills in the last entry and gets the steady-state for kids
	KIDs[:,:,-1] = KIDs[:,:,-2]
	KIDs_ss = KIDs[:,:,-1]

	#Gets the steady-state population share by calculating any number of years ahead necessary until it converges
	Nhatss_old = Nhat_matrix[:,:,-2]
	Nhatss_new = Nhat_matrix[:,:,-1]
	pop_old = N_matrix[:,:,-1]
	pop_new = N_matrix[:,:,-1]

	iter = 0

	while np.max(np.abs(Nhatss_old - Nhatss_new)) > tol:
		Nhatss_old = Nhatss_new
		pop_new[:,0] = np.sum((pop_old[:,:]*FertilityRates[:,:,-1]),axis=1)
		pop_new[:,1:] = pop_old[:,:-1]*(1+ImmigrationRates[:,:-1,-1]-MortalityRates[:,:-1,-1])
		Nhatss_new = pop_new/np.sum(pop_new)
		iter+=1

	print "The SS Population Share converged in", iter, "years beyond year T"

	"""
	for i in range(I):
		plt.plot(range(T+S+1+iter), np.sum(Nhat_matrix[i,:,:], axis=0))
	plt.title("World Super Steady State Total Population")
	plt.show()
	plt.clf()

	for i in range(I):
		plt.plot(range(T+S+1), np.sum(ImmigrationRates[i,:,:], axis=0))
	plt.title("ImmigrationRates")
	plt.show()
	plt.clf()

	for i in range(I):
		plt.plot(range(T+S+1), np.sum(MortalityRates[i,:,:], axis=0))
	plt.title("MortalityRates")
	plt.show()
	plt.clf()

	for i in range(I):
		plt.plot(range(T+S+1), np.sum(FertilityRates[i,:,:], axis=0))
	plt.title("FertilityRates")
	plt.show()
	plt.clf()

	plotDemographics((S,T), [0,1], [iter], str("Final"), Nhat_matrix)
	"""
	
	#Gets labor endowment per household. For now it grows at a constant rate g_A
	lbar = np.cumsum(np.ones(T)*g_A)

	return MortalityRates, Nhat_matrix, KIDs, Nhatss_new, KIDs_ss, lbar

def plotDemographics(params, indexes, years, name, N_matrix, countrynames):
	"""
	Description:
		Plots the population distribution of a given country for any number of specified years

	Inputs:
		index: Integer that indicates which country to plot
		years: List that contains each year to plot
		name: String of the country's name. Used in the legend of the plot

	Outputs:
		None
	"""
	S, T = params

	#If we want to plot total population across years...
	if len(years) == 0:
		for i in range(len(indexes)):
			plt.plot(range(S+1), np.sum(N_matrix[i,:,:], axis=0))
		plt.legend(countrynames)
	#If we want to compare years rather than countries...
	else:
		for y in range(len(years)):
			yeartograph = years[y]
			for i in range(len(indexes)):
			#Checks to make sure we haven't requested to plot a year past the max year
				if yeartograph <= N_matrix.shape[2]:
					plt.plot(range(S+1), N_matrix[i,:,yeartograph])
				else:
					print "\nERROR: WE HAVE ONLY SIMULATED UP TO THE YEAR", T
					time.sleep(15)
		if len(years) == 1:
			plt.legend(countrynames)
		else:
			plt.legend(years)

	plt.title(str(name + " Population Distribution"))
	plt.show()
	plt.clf()

#STEADY STATE FUNCTIONS

def get_kd(assets, kf, Nhat):
	"""
	Description: Calculates the amount of domestic capital that remains in a given country

	Inputs:
		assets = [I,S+1] Matrix, Assets for given time period
		kf     = [I,] Vector, Domestic capital held by foreigners
		Nhat   = [I, S+1] Matrix, Population share for a given year

	Functions called:
		-None

	Objects in Function:
		kd = [I,] Vector, Capital that is used in the domestic country

	Returns kd:
	"""

	kd = np.sum(assets[:,1:-1]*Nhat[:,1:-1], axis=1) - kf

	return kd

def get_n(params):
	"""
	Description: Calculates the total labor supply for each country

	Inputs:
		e    = [I,S+1] or [I,S+1,T+S+1] Matrix, Labor productivities
		Nhat = [I,S+1] or [I,S+1,T+S+1] Matrix, Population share for a given year

	Functions called:
		-None

	Objects in Function:
		n = [I,] or [I,T+S+1] Matrix, Total labor supply
	
	Returns: n
	"""

	e, Nhat = params
	n = np.sum(e*Nhat, axis=1)

	return n

def get_Y(params, kd, n):
	"""
	Description: Calculates the output timepath

	Inputs:
		params = Tuple: Contains the parameters alpha and A
		alpha  = Scalar in (0,1), Production share of capital
		A      = [I,] Vector, Technology level for each country
		kd     = [I,] Vector or [I,T+S+1] Matrix, Domestic held capital stock
		n      = [I,] Vector or [I,S+T+1] Matrix, Aggregate labor supply

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
	Description: Calculates the steady state rental rates or rental rates timepath, depending on where this function is called

	Inputs:
		alpha = Scalar in (0,1), Production share of capital
		Y     = [I,] or [T+S+1] Vector, Either steady state output or output timepath
		kd    = [I,] or [T+S+1] Vector, Either steady state domestic-owned capital or domestic-owned capital timepath

	Objects in Function:
		r = [I,] or [T+S+1] Vector, Either steady state rental rate or rental rate timepath

	Returns: r
	"""

	r = alpha * Y / kd

	return r

def get_w(alpha, Y, n):
	"""
	Description: Calculates the wage rate

	Inputs:
		alpha = Scalar in (0,1), Production share of capital
		Y     = [I,] Vector or [I,S+T+1] Matrix, Aggregate output
		n     = [I,] Vector or [I,S+T+1] Matrix, Aggregate labor supplies

	Objects in Function:
		w = [I,] Vector or [I,S+T+1] Matrix, Wage rates

	Returns: w
	"""

	w = (1-alpha) * Y / n
	return w

def getBequests(params, assets):
	"""
	Description:
		-Gets the value of the bequests given to each generation

	Inputs:
		params            = Tuple that contains the parameters I, S, FirstFertilityAge, FirstDyingAge, Nhat_current, Mortality_current
		I                 = Scalar in [1,7], Number of countries
		S                 = Scalar in [10,80], Number of cohorts
		FirstFertilityAge = Scalar in (0,S), First age when agents bear children
		FirstDyingAge     = Scalar in (0,S), First age when agents die
		Nhat_current      = [I, S+1] Matrix, Population share for a given year
		Mortality_current = [I, S+1] Matrix, Mortality rates for a given year
		assets            = [I, S+1] Matrix, Number of assets per agent in each cohort

	Objects in Function:
		bq                    = [I, S+1] Matrix, Number of bequests each agent receives
		BQ                    = [I,] Vector, Total assets of the people who died this year
		num_bequest_receivers = [I,] Vector, Number of people who are eligible to receive bequests this year
		bq_per_agent          = [I,] Vector, Number of bequests each agent receives

	Returns: bq

	"""

	I, S, FirstFertilityAge, FirstDyingAge, Nhat_current, Mortality_current = params

	#Initializes bequests
	bq = np.zeros((I, S+1))

	#Gets the total assets of the people who died this year
	BQ = np.sum(assets[:,FirstDyingAge:]*Mortality_current[:,FirstDyingAge:]*Nhat_current[:,FirstDyingAge:], axis=1)

	#Distributes the total assets equally among the eligible population for each country
	#NOTE: This will likely change as we get a more complex function for distributing the bequests
	num_bequest_receivers = np.sum(Nhat_current[:,FirstFertilityAge:FirstDyingAge], axis=1)
	bq_per_agent = BQ/num_bequest_receivers
	bq[:,FirstFertilityAge:FirstDyingAge] = np.einsum("i,s->is", bq_per_agent, np.ones(FirstDyingAge-FirstFertilityAge))

	return bq

def get_Gamma(params, w, e):
        """
        Description: Calculates the Gamma used for calculating consumption and other variables.

        Inputs:
            -params: tuple of the needed parameters
            -w[I,T+S+1]:Path of wages
            -e[I,S,T+S+1]: Marginal labor productivities

        Objects in Function:
            -rho: parameter
            -sigma: parameter
            -chi: parameter
            -denom: denominator for equation 3.24

        Outputs:
            -Gamma[I,T+S+1]:


        """

        rho, sigma, chi = params

        #Steady state
        if w.ndim==1 and e.ndim==2:
            denom=np.einsum("i,is->is",w, e[:,:-1])

        #Initial guesses
        elif w.ndim==2 and e.ndim==2:
            denom=np.einsum("it,is -> is", w, e[:,:-1])

        #Timepaths
        elif w.ndim==2 and e.ndim==3:
            denom=np.einsum("it,ist -> it", w, e[:,:-1,:])

        Gamma=((1+chi*(chi/(denom)))**rho**((1-rho*sigma)/rho)*(rho/(rho-1)))**(-1/sigma)


        return Gamma

def get_ck(c,Gamma):
        """
        Description: Calculates the children's consumption based on parent's consumption

        Inputs:
            -c[I,S,T+S+1]: Array of consumption (adults)
            -Gamma[I,T+S+1]: Vector of Gamma

        Objects in Function:
            NONE

        Outputs:
            -ckhat[I,S,T+S+1]: Array of children's consumption

        """
        #print "c shape", c.shape
        #print "Gamma shape",Gamma.shape
        ckhat=c*Gamma


        return ckhat

def get_lhat(params,chat,w,e):
        """
        Description: Calculates the lesiure based on wages, consumption and productivity

        Inputs:
            -params: tuple of the necessary parameter
            -chat: Consumption of adults
            -w: wage timepath
            -e: labor productvities

        Objects in Function:
            -chi: parameter
            -rho: parameter
            -denom: Calculated denominator for 3.26

        Outputs:
            -lhat: timepath of leisure decisions

        """
        rho, chi = params

        if w.ndim==1:
        	denom=np.einsum("i,is->is",w,e[:,:-1])
        elif w.ndim==2:
        	denom=np.einsum("it,ist->ist",w,e[:,:-1])   

        lhat=chat*(chi/denom)**rho

        return lhat

def get_chatvecss(params,w,e,r,assets,KIDs,Gamma,lbar=1):
        """
        Description:

        Inputs:
            -w:
            -e:
            -lbar:
            -params:
            -assets:
            -gA:
            -KID:
            -Gamma:

        Objects in Function:
            -denom:
            -delta:
            -chi:

        Outputs:
            -chat:

        """
        delta, chi, rho, gA = params

        denom=np.einsum("i,is->is",w,e[:,:-1])
        part2=KIDs[:,:-1]*Gamma
        part3=np.einsum("i,is->is",(1+r-delta),assets[:,:-1])
        part4=assets[:,1:]

        chat=(denom+part3-part4*np.exp(gA))/(1+part2+(denom)**rho)


        return chat

def check_feasible(K, Y, w, r, c):
	"""
	Description:Checks the feasibility of the inputs.

	Inputs:
		-K[I,T+S+1]: Capital stock timepath
		-y[I,T+S+1]: Output timepath
		-w[I,T+S+1]: Wage timepath
		-r[I,T+S+1]: Rental rate timepath
		-c[I,T+S+1]: consumption timepath

	Objects in Function:
		NONE

	Outputs:
		-Feasible (Boolean): Whether or not the inputs are feasible.

	"""

	Feasible = True

	if np.any(K<0) or np.any(np.isnan(K)):
		Feasible=False
		print "WARNING! INFEASABLE VALUE ENCOUNTERED IN K!"
		print "The following coordinates have values less than 0:"
		print np.argwhere(K<0)
		print "The following coordinates have nan values"
		print np.argwhere(np.isnan(K))

	if np.any(Y<0) or np.any(np.isnan(Y)):
		Feasible=False
		print "WARNING! INFEASABLE VALUE ENCOUNTERED IN Y!"
		print "Y has dimensions", Y.shape
		print "The following coordinates have values less than 0:"
		print np.argwhere(Y<0)
		print "The following coordinates have nan values"
		nan_coords = np.argwhere(np.isnan(Y))
		print nan_coords
		#for i in range(nan_coords.shape[0]):
			#print "Country", i, "for all agents"
			#print Y[nan_coords[i,0]]

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

def SteadyStateSolution(guess,T, I, S, beta, sigma, delta, alpha, e_ss, A, FirstFertilityAge, FirstDyingAge, Nhat_ss, Mortality_ss, g_A, chi, rho, KIDs_ss):
	"""
	Description: 
		-This is the function that will be optimized by fsolve.

	Inputs:
		-guess[I,S+1]: vector that pieced together from assets and kf.

	Objects in Function:
		-kf[I,]:Foreign capital held by foreigners in each country
		-assets[I,S]: Asset path for each country
		-k[I,]:Capital for each country
		-n[I,]:Labor for each country
		-Y[I,]:Output for each country
		-r[I,]:Rental Rate for each country
		-w[I,]:Wage for each country
		-c_vec[I, S]: Consumption by cohort in each country
		-Euler_c[I, S-1]: Corresponds to (1.16)
		-Euler_r[I,]: Corresponds to (1.17)
		-Euler_kf(Scalar): Corresponds to (1.18)

	Output:
	-all_Euler[I*S,]: Similar to guess, it's a vector that's has both assets and kf.
	
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

	gamparams = (rho, sigma, chi)
	Gamma = get_Gamma(gamparams,w,e_ss,)

	chatparams = (delta, chi, rho, g_A)
	chat_vec = get_chatvecss(chatparams, w, e_ss, r, assets, KIDs_ss, Gamma)

	lhatparams = (rho, chi)
	lhat=get_lhat(lhatparams, chat_vec, w, e_ss)

	ck_vec=get_ck(chat_vec,Gamma)

	K = kd+kf

	Feasible = check_feasible(K, Y, w, r, chat_vec)

	if Feasible == False: #Punishes the the poor choice of negative values in the fsolve
		all_Euler=np.ones((I*S))*999.
		print "Punishing fsolve"
	else:
		#Gets Euler equations
		Euler_ck = ck_vec[:,:-1]** (-sigma) - beta * (ck_vec[:,1:]*np.exp(g_A)) ** (-sigma) * (1 + r[0] - delta)
		Euler_r = r[1:] - r[0]
		Euler_kf = np.sum(kf*np.sum(Nhat_ss, axis=1))

		#Makes a new 1D vector of length I*S that contains all the Euler equations
 		all_Euler = np.append(np.append(np.ravel(Euler_ck), np.ravel(Euler_r)), Euler_kf)

		#print "SS Max Euler Error", np.max(np.absolute(all_Euler))

	return all_Euler

def getSteadyState(params, assets_init, kf_init):
	"""
	Description:
        This takes the initial guess for assets and kf. Since the function
	    returns a matrix, this unpacks the individual parts.
	Inputs:
	    -assets_init[I,S-1]:Intial guess for asset path
	    -kf_init[I]:Initial guess on foreigner held capital  

    Objects in Function:
        -guess[I,S]: A combined matrix that has both assets_init and kf_init
        -ss[S*I,]: The result from optimization.

	Outputs:
	    -assets_ss[I,S-1]:Calculated assets steady state
	    -kf_ss[I,]:Calculated domestic capital owned by foreigners steady state
	    -k_ss[I]: Calculated total capital stock steady state
	    -n_ss[I]: Summed labor productivities steady state
	    -y_ss[I]: Calculated output steady state
	    -r_ss[I]: calculated steady state rental rate
	    -w_ss[I]: calculated steady state wage rate
	    -c_vec_ss[I, S]: Calculated steady state counsumption

	"""
	T, I, S, beta, sigma, delta, alpha, e_ss, A, FirstFertilityAge, FirstDyingAge, Nhat_ss, Mortality_ss, g_A , chi, rho, KIDs_ss = params

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
	r_ss = get_r(alpha, Y_ss, kd_ss)
	w_ss = get_w(alpha, Y_ss, n_ss)
	bqparams = (I, S, FirstFertilityAge, FirstDyingAge, Nhat_ss, Mortality_ss)
	bq_ss = getBequests(bqparams, assets_ss)
	Gammaparams = (rho, sigma, chi)
	Gamma_ss= get_Gamma(Gammaparams, w_ss, e_ss)
	chatparams = (delta, chi, rho, g_A)
	chat_vec_ss = get_chatvecss(chatparams, w_ss, e_ss, r_ss, assets_ss, KIDs_ss, Gamma_ss)

	ck_vec_ss = get_ck(chat_vec_ss, Gamma_ss)

	print "\nSteady State Found!\n"

	return assets_ss, kf_ss, kd_ss, n_ss, Y_ss, r_ss[0], w_ss, chat_vec_ss, ck_vec_ss

#TIMEPATH FUNCTIONS

def get_initialguesses(params, assets_ss, kf_ss, w_ss, r_ss, PrintLoc):
	"""
	Description:
		With the parameters and steady state values, this function creates
		initial guesses in a linear path.

	Inputs:
		-Params (Tuple): Tuple of parameters from Main.py
		-Assets_ss[I,S,T+S+1]: Steady state assets value
		-kf_ss[I,]: Steady State value of foreign owned domestic capital
		-w_ss[I,]: Steady state value of wages
		-r_ss[I,]: Steady state value of rental rate

	Objects in Function:
		-othervariable_params (Tuple): A tuple specifically made for GetOtherVariables


	Outputs:
        -assets_init[I,]: Initial Asset path
        -kf_init[I,]: New initial foreign held capital
        -w_initguess[I,T+S+1]: Initial guess wage timepath
        -r_initguess[I,T+S+1]: Initial guess rental rate timepath
        -k_init[I,]: total capital stock initial guess
        -n_init[I,]: total labor initial guess
        -y_init[I,]: output labor initial guess
        -c_init[I,]: consumption initial guess

	"""
	if PrintLoc: print "Getting initial guesses"
	#Unpacks parameters
	I, S, T, delta, alpha, e_init, A, FirstFertilityAge, FirstDyingAge, Nhat_init, \
		Mortality_init, lbar, g_A, chi, rho, sigma, KIDs_init = params

	#Sets initial assets and kf, start with something close to the steady state
	assets_init = assets_ss*.9
	kf_init = kf_ss*0
	w_initguess = np.ones((I, T+S+1))*.25
	r_initguess = np.ones((T+S+1))*.5

	#Gets initial kd, n, y, r, w, and K
	kd_init = get_kd(assets_init, kf_init, Nhat_init)
	nparams = e_init, Nhat_init
	n_init = get_n(nparams)
	Yparams = (alpha, A)
	Y_init = get_Y(Yparams, kd_init, n_init)
	r_init = get_r(alpha, Y_init, kd_init)
	w_init = get_w(alpha, Y_init, n_init)
	bqparams = (I, S, FirstFertilityAge, FirstDyingAge, Nhat_init, Mortality_init)
	bq_init = getBequests(bqparams, assets_init)
	gamparams = (rho, sigma, chi)
	Gamma_init = get_Gamma(gamparams, w_initguess, e_init)

	chatparams = (delta, chi, rho, g_A)
	chat_init = get_chatvecss(chatparams, w_init, e_init, r_init, assets_init, KIDs_init, Gamma_init, lbar)

	lhatparams = (rho, chi)
	lhat_init = get_lhat(lhatparams,chat_init,w_init,e_init)

	ck_init = get_ck(chat_init, Gamma_init)
	#cparams = (e, delta, bq_init, g_A)
	#c_init = get_cvecss(cparams, w_init, r_init, assets_init)

	#Gets initial guess for rental rate path. This is set up to be linear.
	r_initguess[:T+1] = np.linspace(r_init[0], r_ss, T+1)
	r_initguess[T+1:] = r_initguess[T]

	#Gets initial guess for wage path. This is set up to be linear.
	for i in range(I):
		w_initguess[i, :T+1] = np.linspace(w_init[i], w_ss[i], T+1)
		w_initguess[i,T+1:] = w_initguess[i,T]

	return assets_init, kf_init, w_initguess, r_initguess, kd_init, n_init, Y_init, ck_init, lhat_init

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
        kfPath = np.zeros((I,S+T+1))
        kdPath = np.zeros((I,S+T+1))

        #Gets the domestic-owned capital stock for each country except for the first country
        kdPath[1:,:] = (rpath/alpha)**(1/(alpha-1))*np.einsum("i,is->is", A[1:], n[1:,:])

        #This is using equation 1.13 solved for the foreign capital stock to caluclate the foreign capital stock
        #For everyone except the first country
        kfPath[1:,:] = Kpath[1:,:]-kdPath[1:,:]

        #To satisfy 1.18, the first country's assets is the negative of the sum of all the other countries' assets
        kfPath[0,:] = -np.sum(kfPath[1:,:],axis=0)

		#Making every year beyond t equal to the steady-state
        kfPath[:,T:] = np.einsum("i,s->is", kf_ss, np.ones(S+1))
        
        if PrintLoc: print "Leaving get_foreignK_path"
        return kfPath

def get_lifetime_decisions(params, ck_1, wpath_chunk, rpath_chunk, Gamma_chunk, KIDs_chunk, e_chunk, starting_assets, bq, current_age):
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

	I, S, beta, sigma, delta, lbar, g_A, rho, chi  = params

	num_decisions = S-current_age-1# -1 Because we already have (or have guessed) our starting consumption

	#Initializes the cpath and asset path vectors
	ck_path = np.zeros((I, num_decisions+1))
	c_path = np.zeros((I, num_decisions+1))
	asset_path = np.zeros((I, num_decisions+2))

	#For each country, the cpath and asset path vectors' are the initial values provided.
	ck_path[:,0] = ck_1
	asset_path[:,0] = starting_assets

	#Based on the individual chunks, these are the households choices
	for p in range(1,num_decisions+1):
		ck_path[:,p] = ((beta * (1 + rpath_chunk[p] - delta))**(1/sigma) * ck_path[:,p-1])*np.exp(-g_A)
		c_path[:,p-1] = ck_path[:,p-1]/Gamma_chunk[:,p-1]
		asset_path[:,p] = (wpath_chunk[:,p-1]*e_chunk[:,p-1]*lbar[p-1] + (1 + rpath_chunk[p-1] - delta)*asset_path[:,p-1] + bq[:,p-1]\
						  - c_path[:,p-1]*(1+KIDs_chunk[:,p-1]*Gamma_chunk[:,p-1]+(chi/wpath_chunk[:,p-1]/e_chunk[:,p-1])**rho)\
						  )*np.exp(-g_A)
	
	asset_path[:,p+1] = (wpath_chunk[:,p]*e_chunk[:,p]*lbar[p] + (1 + rpath_chunk[p] - delta)*asset_path[:,p] + bq[:,p] - c_path[:,p]*(1+KIDs_chunk[:,p]*Gamma_chunk[:,p]+(chi/wpath_chunk[:,p]/e_chunk[:,p])**rho))*np.exp(-g_A)

	return ck_path, asset_path

def find_optimal_starting_consumptionsK(ck_1, wpath_chunk, rpath_chunk, Gamma_chunk, KIDs_chunk, epath_chunk, starting_assets, bq, current_age, params):
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
	ck_path, assets_path = get_lifetime_decisions(params, ck_1, wpath_chunk, rpath_chunk, Gamma_chunk, KIDs_chunk, epath_chunk, starting_assets, bq, current_age)

        if np.any(ck_path<0):
            print "WARNING! The fsolve for initial optimal consumption guessed a negative number"
            ck_path=np.ones(I)*9999.

	Euler = np.ravel(assets_path[:,-1])
	#print "Max Euler Error in find_optimal_starting_consumptionsK", np.max(np.absolute(Euler))

	return Euler

def get_cons_assets_matrix(params, wpath, rpath, starting_assets, PrintLoc):
	if PrintLoc: print "Entering get_cons_assets_matrix"

	I, S, T, T_1, beta, sigma, delta, e, FirstFertilityAge, FirstDyingAge, Nhat, KIDs, MortalityRates, lbar, g_A, rho, chi = params

	Gammaparams=(rho, sigma, chi)
	Gamma=get_Gamma(Gammaparams,wpath,e)

	#Initializes timepath variables
	c_timepath = np.zeros((I,S,S+T+1))
	ck_timepath = np.zeros((I,S,S+T+1))
	a_timepath = np.zeros((I, S+1, S+T+1)) #I,S+1,S+T+1
	a_timepath[:,:,0]=starting_assets
	bq_timepath = np.zeros((I, S+1, S+T+1))

	c_timepath[:,S-1,0] = wpath[:,0]*e[:,S-1,0] + (1 + rpath[0] - delta)*a_timepath[:,S-1,0] - a_timepath[:,S,1]*np.exp(g_A) \
							/(1 + KIDs[:,S-1,0]*Gamma[:,0])

	ck_timepath[:,S-1,0] = c_timepath[:,S-1,0]*Gamma[:,0]

	#print "Initial matrices"
	#print "Consumption"
	#print np.round(np.transpose(c_timepath[0,:,:2]), decimals=3)
	#print "Assets"
	#print np.round(np.transpose(a_timepath[0,:,:2]), decimals=3)
	#print "Bequests"
	#print np.round(np.transpose(bq_timepath[0,:,:2]), decimals=3)

	#Fills the upper triangle (including the main diagonal) by iterating by number of periods until death
	if PrintLoc: print "Entering upper triangle loop"
	for p in range(1,S):

		#We are only doing this for all generations alive in time t=0
		t = 0

		#Getting the current age of the agent
		current_age = S-p-1

		#Uses the previous generation's consumption at age s to get the value for our guess
		ck_guess = (ck_timepath[:,current_age+1,t]/((beta*(1+rpath[t]-delta))**(1/sigma)))/np.exp(g_A)

		agent_assets = starting_assets[:,current_age]

		agent_KIDs = np.diagonal(KIDs[:,current_age:,t:], axis1=1, axis2=2)

		#Gets the bequests this agent will recieve in his remaining lifetime
		agent_bq = np.diagonal(bq_timepath[:,current_age:,t:], axis1=1, axis2=2)

		#Gets labor productivities this agent will recieve in his remaining lifetime
		agent_e = np.diagonal(e[:,current_age:,t:], axis1=1, axis2=2)

		#Gets optimal initial consumption beginning in the current age of the agent using chunks of w and r that span the lifetime of the given generation
		household_params = (I, S, beta, sigma, delta, lbar, g_A, rho, chi)
		opt_consumpK = opt.fsolve(find_optimal_starting_consumptionsK, ck_guess, args = \
			(wpath[:,t:t+p+1], rpath[t:t+p+2], Gamma[:,t:t+p+1], agent_KIDs, agent_e, agent_assets, agent_bq, current_age, household_params))

		#Gets optimal timepaths beginning initial consumption and starting assets
		ckpath_indiv, apath_indiv = get_lifetime_decisions\
			(household_params, opt_consumpK, wpath[:,t:t+p+1], rpath[t:t+p+2], Gamma[:,t:t+p+1], agent_KIDs, agent_e, agent_assets, agent_bq, current_age)

		print "optimal c_k found! For each country this agent's final assets are", apath_indiv[:,-1]

		for i in xrange(I):
			np.fill_diagonal(ck_timepath[i,current_age:,:], ckpath_indiv[i,:])
			np.fill_diagonal(a_timepath[i,current_age:,:], apath_indiv[i,:])

		bq_params = (I, S, FirstFertilityAge, FirstDyingAge, Nhat[:,:,p-1], MortalityRates[:,:,p-1])
		bq_timepath[:,:,p-1] = getBequests(bq_params, a_timepath[:,:,p-1])

		#print "p =", p, "current_age =", current_age
		#print "Consumption"
		#print np.round(np.transpose(c_timepath[0,:,:p+2]), decimals=3)
		#print "c_guess", c_guess
		#print "Assets"
		#print np.round(np.transpose(a_timepath[0,:,:p+2]), decimals=3)
		#print "Bequests"
		#print np.round(np.transpose(bq_timepath[0,:,:p+2]), decimals=3)
		#print np.round(np.diagonal(bq_timepath[:,current_age:,:], axis1=1, axis2=2), decimals=3).shape
		#print np.round(np.diagonal(bq_timepath[:,current_age:,:], axis1=1, axis2=2), decimals=3)

	if PrintLoc: print "Entering non-upper triangle loop"
	#Fills everything except for the upper triangle (excluding the main diagonal)
	for t in xrange(1,T):
		current_age = 0
		p = S-current_age

		agent_assets = np.zeros((I))

		#Uses the previous generation's consumption at age s to get the value for our guess
		ck_guess = ck_timepath[:,current_age,t-1]

		#Gets the bequests this agent will recieve in his remaining lifetime
		agent_bq = np.diagonal(bq_timepath[:,current_age:,t:], axis1=1, axis2=2)

		agent_KIDs = np.diagonal(KIDs[:,current_age:,t:], axis1=1, axis2=2)

		#Gets labor productivities this agent will recieve in his remaining lifetime
		agent_e = np.diagonal(e[:,current_age:,t:], axis1=1, axis2=2)

		opt_consumpK = opt.fsolve(find_optimal_starting_consumptionsK, ck_guess, args = \
			(wpath[:,t:t+p+1], rpath[t:t+p+2], Gamma[:,t:t+p+1], agent_KIDs, agent_e, agent_assets, agent_bq, current_age, household_params))

		#Gets optimal timepaths beginning initial consumption and starting assets
		ckpath_indiv, apath_indiv = get_lifetime_decisions\
			(household_params, opt_consumpK, wpath[:,t:t+p+1], rpath[t:t+p+2], Gamma[:,t:t+p+1], agent_KIDs, agent_e, agent_assets, agent_bq, current_age)

		print "optimal c_k found! For each country this agent's final assets are", apath_indiv[:,-1]

		for i in range(I):
			np.fill_diagonal(ck_timepath[i,:,t:], ckpath_indiv[i,:])
			np.fill_diagonal(a_timepath[i,:,t:], apath_indiv[i,:])

		if t >= T_1:
			temp_t = T_1
		else:
			temp_t = t

		bq_params = (I, S, FirstFertilityAge, FirstDyingAge, Nhat[:,:,temp_t+S-2], MortalityRates[:,:,temp_t+S-2])
		bq_timepath[:,:,t+S-2] = getBequests(bq_params, a_timepath[:,:,temp_t+S-2])

		#print "t = ", t, "current_age =", current_age
		#print "Consumption"
		#print np.round(np.transpose(c_timepath[0,:,:p+2]), decimals=3)
		#print "c_guess", np.round(c_guess, decimals=3)
		#print "Assets"
		#print np.round(np.transpose(a_timepath[0,:,:p+2]), decimals=3)
		#print "Bequests"
		#print np.round(np.transpose(bq_timepath[0,:,:p+2]), decimals=3)

	c_timepath = np.einsum("ist,it->ist", ck_timepath, 1/Gamma)
	lhatparams = (rho, chi)
	lhat_timepath = get_lhat(lhatparams,c_timepath,wpath,e)

	if PrintLoc: print "Leaving get_cons_assets_matrix"
	return c_timepath, ck_timepath, a_timepath, lhat_timepath

def get_wpathnew_rpathnew(params, wpath, rpath, starting_assets, kd_ss, kf_ss, w_ss, r_ss, PrintLoc):
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

	I, S, T, T_1, beta, sigma, delta, alpha, e, A, FirstFertilityAge, FirstDyingAge, Nhat, KIDs, MortalityRates, lbar, g_A, rho, chi = params

	ca_params = (I, S, T, T_1, beta, sigma, delta, e, FirstFertilityAge, FirstDyingAge, Nhat, KIDs, MortalityRates, lbar, g_A, rho, chi)

	c_timepath, ck_timepath, a_timepath, lhat_timepath = get_cons_assets_matrix(ca_params, wpath, rpath, starting_assets, PrintLoc)

	#Calculates the total amount of capital in each country
	Kpath=np.sum(a_timepath*Nhat, axis=1)

	#Calculates Aggregate Consumption
	Cpath=np.sum(c_timepath, axis=1)

	#After time period T, the total capital stock and total consumption is forced to be the steady state
	Kpath[:,T:] = np.einsum("i,t->it", kd_ss+kf_ss, np.ones(S+1))
	Cpath[:,T:] = np.einsum("i,t->it", Cpath[:,T-1], np.ones(S+1))

	#Gets the foriegned owned capital
	kf_params = (I, S, T, alpha, e, A, Nhat)
	kfpath = get_foreignK_path(kf_params, Kpath, rpath, kf_ss, PrintLoc)

	#Based on the overall capital path and the foreign owned capital path, we get new w and r paths.
	kdpath = Kpath - kfpath
	nparams = (e, Nhat)
	npath = get_n(nparams)
	Yparams = (alpha, A)
	Ypath = get_Y(Yparams, kdpath, npath)
	rpath_new = get_r(alpha, Ypath[0], kdpath[0])
	wpath_new = get_w(alpha, Ypath, npath)


	#Checks to see if any of the timepaths have negative values or nans
	Feasible = check_feasible(Kpath, Ypath, wpath, rpath, c_timepath)

	if PrintLoc: print "Leaving get_wpathnew_rpathnew"
	return wpath_new, rpath_new, Cpath, Kpath, Ypath

def get_Timepath(params, wstart, rstart, assets_init, kd_ss, kf_ss, w_ss, r_ss, PrintLoc):

    I, S, T, T_1, beta, sigma, delta, alpha, e, A, FirstFertilityAge, FirstDyingAge, Nhat, KIDs, MortalityRates, lbar, g_A, distance, diff, xi, MaxIters, rho, chi = params

    Iter=1 #Serves as the iteration counter
    wr_params = (I, S, T, T_1, beta, sigma, delta, alpha, e, A, FirstFertilityAge, FirstDyingAge, Nhat, KIDs, MortalityRates, lbar, g_A, rho, chi)

    while distance>diff and Iter<MaxIters: #The timepath iteration runs until the distance gets below a threshold or the iterations hit the maximum
            wpath_new, rpath_new, Cpath_new, Kpath_new, Ypath_new = \
            get_wpathnew_rpathnew(wr_params, wstart, rstart, assets_init, kd_ss, kf_ss, w_ss, r_ss, PrintLoc)

            dist1=sp.linalg.norm(wstart-wpath_new,2) #Norm of the wage path
            dist2=sp.linalg.norm(rstart-rpath_new,2) #Norm of the intrest rate path
            distance=max([dist1,dist2]) #We take the maximum of the two norms to get the distance

            print "Iteration:",Iter,", Norm Distance: ", distance#, "Euler Error, ", EError
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

def CountryLabel(Country): #Activated by line 28
    '''
    Description: 
        Converts the generic country label given for the graphs and converts it to a proper name

    Inputs:
        -Country (String): This is simply the generic country label

    Objects in Function:
        -NONE
    Outputs:
        -Name (String): The proper name of the country which you decide. Make sure the number of country names lines
            up with the number of countries, otherwise, the function will not proceed.
    '''

    #Each country is given a number
    if Country=="Country 0":
        Name="United States"
    if Country=="Country 1":
        Name="Europe"
    if Country=="Country 2":
        Name="Japan"
    if Country=="Country 3":
        Name="China"
    if Country=="Country 4":
        Name="India"
    if Country=="Country 5":
        Name="Russia"
    if Country=="Country 6":
        Name="Korea"

    #Add More Country labels here

    return Name

def plotTimepaths(I, S, T, wpath, rpath, cpath, kpath, Ypath, CountryNamesON):

    for i in xrange(I): #Wages
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,T),wpath[i,:T], label=label1)
    plt.title("Time path for Wages")
    plt.ylabel("Wages")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()

    #Rental Rates
    label1='Global Interest Rate'    
    plt.plot(np.arange(0,T),rpath[:T], label=label1)
    plt.title("Time path for Rental Rates")
    plt.ylabel("Rental Rates")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()

    for i in xrange(I): #Aggregate Consumption
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,S+T+1),cpath[i,:],label=label1)
    plt.title("Time Path for Aggregate Consumption")
    plt.ylabel("Consumption Level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()

    for i in xrange(I): #Aggregate Capital Stock
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,T),kpath[i,:T],label=label1)
    plt.title("Time path for Capital Path")
    plt.ylabel("Capital Stock level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()

    for i in xrange(I):
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,T),Ypath[i,:T],label=label1)
    plt.title("Time path for Output")
    plt.ylabel("Output Stock level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()
