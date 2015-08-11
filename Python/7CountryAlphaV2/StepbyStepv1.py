from __future__ import division
import numpy as np
import scipy as sp
import scipy.optimize as opt
from matplotlib import pyplot as plt

#DEMOGRAPHICS FUNCTIONS

def getDemographics(params):
	"""
	Description:
		-Imports data from csv files for initial populations, fertility rates, mortality rates, and net migrants. 
		-Stores these data sets in their respective matrices, and calculates population distribuitons through year T.
		-NOTE: FOR NOW THIS FUNCTION ONLY USES DATA FOR THE USA. NEEDS TO EVENTUALLY ADD MORE COUNTRIES

	Inputs:
		-None, but uses the global variables T, T_1, StartFertilityAge, EndFertilityAge, StartDyingAge, and MaxImmigrantAge

	Objects in Function:
		-USAPopdata: (S+1) vector that has the initial population of the U.S straight from the csv
		-USAFertdata: (T_1,EndFertilityAge+2-StartFertilityAge) vector that has U.S. fertility straight from the csv
		-USAMortdata: (T_1,S+1-StartDyingAge) vector that has U.S. mortality straight from the csv
		-USAMigdata: (MaxImmigrantAge) vector that contains the number of net U.S. migrants straight from the csv
		-g_N: (T) vector that contains the exogenous population growth rates
		-g_A: Constant that represents the technical growth rate
		-l_endowment: (T) vector labor endowment per household
		-f_bar: (I) vector that represents the fertility rate after period T_1
		-p_bar: (I) vector that represents the mortality rate after period T_1
		-m_bar: (I) vector that represents the immigration rate after period T_1

	Output:
		-FertilityRates: Numpy array that contains fertilty rates for all countries, ages, and years
		-MortalityRates: Numpy array that contains mortality rates for all countries, ages, and years
		-Migrants: Numpy array that contains net migration for all countries and ages
		-N_matrix: Numpy array that contains population numbers for all countries, ages, and years
		-Nhat matrix: Numpy array that contains the population percentage for all countries, ages, and years
	"""

	I, S, T, T_1, StartFertilityAge, EndFertilityAge, StartDyingAge, MaxImmigrantAge, g_A, PrintAges = params
	if PrintAges:
		print "T =", T
		print "T_1", T_1
		print "StartFertilityAge", StartFertilityAge
		print "EndFertilityAge", EndFertilityAge
		print "StartDyingAge", StartDyingAge
		print "MaxImmigrantAge", MaxImmigrantAge

	#Imports and scales data for the USA. Imports a certain number of generations according to the value of S
	USAPopdata = np.loadtxt(("Data_Files/population.csv"),delimiter=',',skiprows=1, usecols=[1])[:S+1]*1000
	USAFertdata = np.loadtxt(("Data_Files/usa_fertility.csv"),delimiter=',',skiprows=1, usecols=range(1,EndFertilityAge+2-StartFertilityAge))[48:48+T_1,:]
	USAMortdata = np.loadtxt(("Data_Files/usa_mortality.csv"),delimiter=',',skiprows=1, usecols=range(1,S+1-StartDyingAge))[:T_1,:]
	USAMigdata = np.loadtxt(("Data_Files/net_migration.csv"),delimiter=',',skiprows=1, usecols=[1])[:MaxImmigrantAge]*100

	#Initializes demographics matrices
	N_matrix = np.zeros((I, S+1, T))
	Nhat_matrix = np.zeros((I, S+1, T))
	FertilityRates = np.zeros((I, S+1, T))
	MortalityRates = np.zeros((I, S+1, T))
	Migrants = np.zeros((I, S+1, T))
	g_N = np.zeros(T)

	#NOTE: For now we set fertility, mortality, number of migrants, and initial population the same for all countries. 
	
	#Sets initial total population (N_matrix), percentage of total world population (Nhat_matrix)
	N_matrix[:,:,0] = np.tile(USAPopdata, (I, 1))
	Nhat_matrix[:,:,0] = N_matrix[:,:,0]/np.sum(N_matrix[:,:,0])

	#Fertility Will be equal to 0 for all ages that don't bear children
	FertilityRates[:,StartFertilityAge:EndFertilityAge+1,:T_1] = np.einsum("ts,i->ist", USAFertdata, np.ones(I))

	#Mortality be equal to 0 for all young people who aren't old enough to die
	MortalityRates[:,StartDyingAge:-1,:T_1] = np.einsum("ts,it->ist", USAMortdata, np.ones((I,T_1)))

	#The last generation dies with probability 1
	MortalityRates[:,-1,:] = np.ones((I, T))

	#The number of migrants is the same for each year
	Migrants[:,:MaxImmigrantAge,:T_1] = np.einsum("s,it->ist", USAMigdata, np.ones((I,T_1)))

	#Gets steady-state values
	f_bar = FertilityRates[:,:,T_1-1]
	p_bar = MortalityRates[:,:,T_1-1]
	m_bar = Migrants[:,:,T_1-1]

	#Set to the steady state for every year beyond year T_1
	FertilityRates[:,:,T_1:] = np.tile(np.expand_dims(f_bar, axis=2), (1,1,T-T_1))
	MortalityRates[:,:,T_1:] = np.tile(np.expand_dims(p_bar, axis=2), (1,1,T-T_1))
	Migrants[:,:,T_1:] = np.tile(np.expand_dims(m_bar, axis=2), (1,1,T-T_1))

	#Gets the initial immigration rate
	ImmigrationRate = Migrants[:,:,0]/N_matrix[:,:,0]

	#Gets initial world population growth rate
	g_N[0] = np.sum(Nhat_matrix[:,:,0]*(FertilityRates[:,:,0] + ImmigrationRate - MortalityRates[:,:,0]))

	#Calculates population numbers for each country
	for t in range(1,T):
		#Gets the total number of children and and percentage of children and stores them in generation 0 of their respective matrices
		#See equations 2.1 and 2.10
		N_matrix[:,0,t] = np.sum((N_matrix[:,:,t-1]*FertilityRates[:,:,t-1]),axis=1)
		Nhat_matrix[:,0,t] = np.exp(-g_N[t-1])*np.sum((Nhat_matrix[:,:,t-1]*FertilityRates[:,:,t-1]),axis=1)

		#Finds the immigration rate for each year
		ImmigrationRate = Migrants[:,:,t-1]/N_matrix[:,:,t-1]

		#Gets the population and percentage of population for the next year, taking into account immigration and mortality
		#See equations 2.2 and 2.11
		N_matrix[:,1:,t] = N_matrix[:,:-1,t-1]*(1+ImmigrationRate[:,:-1]-MortalityRates[:,:-1,t-1])
		Nhat_matrix[:,1:,t] = np.exp(-g_N[t-1])*Nhat_matrix[:,:-1,t-1]*(1+ImmigrationRate[:,:-1]-MortalityRates[:,:-1,t-1])
		
		#Gets the growth rate for the next year
		g_N[t] = np.sum(Nhat_matrix[:,:,t]*(FertilityRates[:,:,t] + ImmigrationRate - MortalityRates[:,:,t]), axis=(0, 1))
		#print np.sum(Nhat_matrix[:,:,t])#This should be equal to 1.0

	#Gets labor endowment per household. For now it grows at a constant rate g_A
	l_endowment = np.cumsum(np.ones(T)*g_A)

	return FertilityRates, MortalityRates, Migrants, N_matrix, Nhat_matrix

def plotDemographics(params, index, years, name, N_matrix):
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

	for y in range(len(years)):
		yeartograph = years[y]
		#Checks to make sure we haven't requested to plot a year past the max year
		if yeartograph <= T:
			plt.plot(range(S+1), N_matrix[index,:,yeartograph])
		else:
			print "\nERROR: WE HAVE ONLY SIMULATED UP TO THE YEAR", T
			time.sleep(15)

	plt.title(str(name + " Population Distribution"))

	plt.legend(years)
	plt.show()
	plt.clf()

def getBequests(params, assets_old):
	"""
	Description:
		-Gets the value of the bequests given to each generation

	Inputs:
		-assets: Assets for each generation in a given year
		-current_t: Integer that indicates the current year. Used to pull information from demographics global matrices like FertilityRates

	Objects in Function: 
		-BQ: T
		-num_bequest_receivers:
		-bq_Distribution:

	Output:
		-bq: Numpy array that contains the number of bequests for each generation in each country.

	"""

	I, S, T, StartFertilityAge, StartDyingAge, pop_old, pop_working, current_mort = params

	#Initializes bequests
	bq = np.zeros((I, S+1))

	#Gets the total assets of the people who died this year
	BQ = np.sum(assets_old*current_mort*pop_old, axis=1)

	#Distributes the total assets equally among the eligible population for each country
	#NOTE: This will likely change as we get a more complex function for distributing the bequests
	num_bequest_receivers = np.sum(pop_working, axis=1)
	bq_Distribution = BQ/num_bequest_receivers
	bq[:,StartFertilityAge:StartDyingAge+1] = np.einsum("i,s->is", bq_Distribution, np.ones(StartDyingAge+1-StartFertilityAge))

	return bq

def hatvariables(Kpathreal, kfpathreal, Nhat_matrix):

	#THIS FUNCTION HAS EQUATIONS 2.13-2.16 AND 2.19-2.20, BUT STILL NEEDS TO BE INCORPORATED INTO THE REST OF THE MODEL TO COMPLETELY TEST

	#We are only using up until T periods rather than T+S+1 since Nhat only goes out to T
	Kpath = Kpathreal[:,:T]
	kfpath = kfpathreal[:,:T]
	temp_e = np.ones((I, S+1, T))#THIS SHOULD ONLY BE UNTIL WE GET S GENERATIONS RATHER THAN S-1

	n = np.sum(temp_e[:,:,:T]*Nhat_matrix, axis=1)
	ypath = (Kpath**alpha) * (np.einsum("i,it->it", A, n)**(1-alpha))
	rpath = alpha * ypath / Kpath
	wpath = (1-alpha) * ypath / n
	"""
	#NOTE:This goes in the get_householdchoices_path function

	c_path = np.zeros((I, S))
	asset_path = np.zeros((I, S+1))

	c_path[:,0] = c_1
	asset_path[:,0] = starting_assets

	for s in range(1,S):
		c_path[:,s] = ((beta * (1 + rpath_chunk[:,s] - delta))**(1/sigma) * c_path[:,s-1])/np.exp(g_A)
		asset_path[:,s] = (wpath_chunk[:,s]*e[:,0,s-1] + (1 + rpath_chunk[:,s-1] - delta)*asset_path[:,s-1] + bq_chunk - c_path[:,s-1])/np.exp(g_A)

	asset_path[:,s+1] = wpath_chunk[:,s]*e_chunk[:,s] + (1 + rpath_chunk[:,s] - delta)*asset_path[:,s] - c_path[:,s]

	"""

#STEADY STATE FUNCTIONS

def getOtherVariables(params, assets, kf):
	"""
	Description:
		-Based on the assets and capital held by foreigners, we calculate the other variables.

	Inputs:
		-assets [I,S+1]: Matrix of assets
		-kf[I, ]: Domestic capital held by foreigners

	Objects in function:
		-NONE that aren't already listed

	Output:
		-k[I,]: Capital (1.10)
		-n[I,]: Sum of labor productivities (1.11)
		-y[I,]: Output (1.12)
		-r[I,]: Rental Rate (1.13)
		-w[I,]: Wage (1.14)
		-c_vec[I,S]: Vector of consumptions (1.15)
	"""
        Feasible=True

	I, delta, alpha, e, A = params

	#You have 0 assets when you're born, and 0 when you die
	assets = np.column_stack((np.zeros(I), assets, np.zeros(I)))

	k = np.sum(assets[:,1:-1], axis=1) - kf
	n = np.sum(e[:,:,0], axis=1)
	y = (k**alpha) * ((A*n)**(1-alpha))
	r = alpha * y / k
	w = (1-alpha) * y / n

        if np.any(k<0):
            Feasible=False
            print "WARNING! INFEASABLE VALUE ENCOUNTERED IN k!"
            print "The following coordinates have infeasible values:"
            print np.argwhere(k<0)

        if np.any(n<0):
            Feasible=False
            print "WARNING! INFEASABLE VALUE ENCOUNTERED IN n!"
            print "The following coordinates have infeasible values:"
            print np.argwhere(n<0)

        if np.any(y<0):
            Feasible=False
            print "WARNING! INFEASABLE VALUE ENCOUNTERED IN y!"
            print "The following coordinates have infeasible values:"
            print np.argwhere(y<0)

        if np.any(r<0):
            Feasible=False
            print "WARNING! INFEASABLE VALUE ENCOUNTERED IN r!"
            print "The following coordinates have infeasible values:"
            print np.argwhere(r<0)

        if np.any(w<0):
            Feasible=False
            print "WARNING! INFEASABLE VALUE ENCOUNTERED IN w!"
            print "The following coordinates have infeasible values:"
            print np.argwhere(w<0)



	c_vec = np.einsum("i, is -> is", w, e[:,:,0]) + np.einsum("i, is -> is",(1 + r - delta) \
	, assets[:,:-1]) - assets[:,1:]

        if np.any(c_vec<0):
            Feasible=False
            print "WARNING! INFEASABLE VALUE ENCOUNTERED IN c_vec!"
            print "The following coordinates have infeasible values:"
            print np.argwhere(c_vec<0)


	return k, n, y, r, w, c_vec

def SteadyStateSolution(guess, I, S, beta, sigma, delta, alpha, e, A):
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
		-y[I,]:Output for each country
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

	#Sets kf as the last element of the guess vector for each country and assets as everything else
	assets = guess[:,:-1]
	kf = guess[:,-1]

	#Based on the assets and kf, we get the other vectors
	params = I, delta, alpha, e, A
	k, n, y, r, w, c_vec = getOtherVariables(params, assets, kf)

	#Gets Euler equations
	Euler_c = c_vec[:,:-1] ** (-sigma) - beta * c_vec[:,1:] ** (-sigma) * (1 + r[0] - delta)
	Euler_r = r[1:] - r[0]
	Euler_kf = np.sum(kf)

        if np.any(Euler_c<0): #Punishes the the poor choice of negative values in the fsolve
            Euler_c=np.ones(I)*9999.

	#Makes a new 1D vector of length I*S that contains all the Euler equations
	all_Euler = np.append(np.append(np.ravel(Euler_c), np.ravel(Euler_r)), Euler_kf)

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
	    -kf_ss[I,]:Calculated foreign capital
	    -k_ss[I]: ASK JEFF
	    -n_ss[I]: steady-state labor something
	    -y_ss[I]: steady-state labor something
	    -y_ss[I]: steady-state labor something
	    -y_ss[I]: steady-state labor something
	    -y_ss[I, S]: steady-state consumption vector
	"""
	I, S, beta, sigma, delta, alpha, e, A = params

	#Merges the assets and kf together into one matrix that can be inputted into the fsolve function
	guess = np.column_stack((assets_init, kf_init))

	#Solves for the steady state
	solver_params = (I, S, beta, sigma, delta, alpha, e, A)
	ss = opt.fsolve(SteadyStateSolution, guess, args=solver_params)

	print "\nSteady State Found!\n"

	#Reshapes the ss code
	ss = np.array(np.split(ss, I))

    #Breaks down the steady state matrix into the two separate assets and kf matrices.
	assets_ss = ss[:,:-1]
	kf_ss = ss[:,-1]

	#Gets the other steady-state values using assets and kf
	othervariable_params = I, delta, alpha, e, A
	k_ss, n_ss, y_ss, r_ss, w_ss, c_vec_ss = getOtherVariables(othervariable_params, assets_ss, kf_ss)

	return assets_ss, kf_ss, k_ss, n_ss, y_ss, r_ss, w_ss, c_vec_ss

#TIMEPATH FUNCTIONS

def get_initialguesses(params, assets_ss, kf_ss, w_ss, r_ss):

	I, S, T, delta, alpha, e, A = params

	#Sets initial assets and kf, start with something close to the steady state
	assets_init = assets_ss*.95
	kf_init = kf_ss*.95
	w_initguess = np.zeros((I, T+S+1))
	r_initguess = np.ones((I, T+S+1))*.5

	#Gets initial k, n, y, r, w, and c
	othervariable_params = 	I, delta, alpha, e, A
	k_init, n_init, y_init, r_init, w_init, c_init = getOtherVariables(othervariable_params, assets_init, kf_init)

	#Gets initial guess for w and r paths. This is set up to be linear.
	for i in range(I):
		w_initguess[i, :T+1] = np.linspace(w_init[i], w_ss[i], T+1)
		r_initguess[i, :T+1] = np.linspace(r_init[i], r_ss[i], T+1)
		w_initguess[i,T+1:] = w_initguess[i,T]
		r_initguess[i,T+1:] = r_initguess[i,T]

	return assets_init, kf_init, w_initguess, r_initguess, k_init, n_init, y_init, c_init

def get_prices(params, Kpath, kf_tpath, w_ss, r_ss):
	"""
	Description:
		Based on the given paths, the paths for wages and rental rates are figured
		out based on equations 1.4-1.5

	Inputs:
		-assets_tpath: Asset timepath
		-kf_tpath: Foreign held capital timepath.

	Objects in Functions:
		-Kdpath[I, S+T+1]:Path of domestic owned capital stock

	Outputs:
	-wpath[I, S+T+1]: Wage path
	-rpath[I, S+T+1]: Rental rate path
	-ypath[I, S+T+1]: Output path

	"""
        Feasible=True
	S, T, alpha, e, A = params

	Kdpath=Kpath-kf_tpath

	#Gets non-price variables needed to caluclate prices
	n = np.sum(e, axis=1) #Sum of the labor productivities

	#Gets the path for output, y
	ypath = (Kpath**alpha) * (np.einsum("i,is->is", A, n)**(1-alpha))

	#Gets prices
	rpath = alpha * ypath / Kpath
	wpath = (1-alpha) * ypath / n

	#Tiles the steady-state for each year beyond the steady state
	rpath[:,T:] = np.einsum("i,s->is", r_ss, np.ones(S+1))
	wpath[:,T:] = np.einsum("i,s->is", w_ss, np.ones(S+1))
        if np.any(rpath<0):
            Feasible=False
            print "WARNING! INFEASABLE VALUE ENCOUNTERED IN rpath!"
            print "The following coordinates have infeasible values:"
            print np.argwhere(rpath<0)

        if np.any(wpath<0):
            Feasible=False
            print "WARNING! INFEASABLE VALUE ENCOUNTERED IN wpath!"
            print "The following coordinates have infeasible values:"
            print np.argwhere(wpath<0)

        if np.any(ypath<0):
            Feasible=False
            print "WARNING! INFEASABLE VALUE ENCOUNTERED IN ypath!"
            print "The following coordinates have infeasible values:"
            print np.argwhere(ypath<0)


	#Returns only the first country's interest rate, since they should be the same in theory
	return wpath, rpath, ypath

def get_foreignK_path(params, Kpath, rpath, k_ss, kf_ss):
        """
        Description:
           This calculates the timepath of the foreign capital stock. This is based on equation (1.12 and 1.13).
        Inputs:
            apath: Asset path, from our calculations
            rpath: Rental Rate path, also from our calculation
        
        Objects in Function:
            kDpath[I,S+T+1]: Path of domestic owned capital
            n[I,S+T+1]: Path of total labor
            kf_ss[I,]: Calculated from the steady state. 
            A[I,]: Parameters from above

        Outputs:
            kfPath[I,S+T+1]: Path of domestic capital held by foreigners.
        """
        I, S, T, alpha, e, A = params

        #Sums the labor productivities across cohorts
        n = np.sum(e, axis=1)

        #Declares the array that will later be used.
        kfPath=np.zeros((I,S+T+1))
        kDPath=np.zeros((I,S+T+1))

        #Gets the domestic-owned capital stock for each country except for the first country
        kDPath[1:,:]=(rpath[1:,:]/alpha)**(1/(alpha-1))*np.einsum("i,is->is", A[1:], n[1:,:])

        #This is using equation 1.13 solved for the foreign capital stock to caluclate the foreign capital stock
        kfPath=Kpath-kDPath

        #To satisfy 1.18, the first country's assets is the negative of the sum of all the other countries' assets
        kfPath[0,:]=-np.sum(kfPath,axis=0)

		#Making every year beyond t equal to the steady-state
        kfPath[:,T:] = np.einsum("i,s->is", kf_ss, np.ones(S+1))
        
        return kfPath

def get_lifetime_decisions(params, c_1, wpath_chunk, rpath_chunk, e_chunk, starting_assets, current_s):
	"""
	Description:
		This solves for equations 1.15 and 1.16 in the StepbyStep pdf for a certain generation
	Inputs:
		-c_1: Initial consumption (not necessarily for the year they were born)
		-wpath_chunk: Wages of an agents lifetime, a section of the timepath
		-rpath_chunk: Rental rate of an agents lifetime, a section of the timepath
		-e_chunk: Worker productivities of an agents lifetime, a section of the global matrix
		-starting_assets: Initial assets of the agent. Will be 0s if we are beginning in the year the agent was born
		-current_s: Current age of the agent

		Objects in Function:
			-NONE

	Outputs:
		-c_path[I, S]: Path of consumption until the agent dies
		-asset_path[I, S+1]: Path of assets until the agent dies
	"""

	I, S, beta, sigma, delta = params

	#Initializes the cpath and asset path vectors
	c_path = np.zeros((I, S))
	asset_path = np.zeros((I, S+1))

	#For each country, the cpath and asset path vectors' are the initial values provided.
	c_path[:,0] = c_1
	asset_path[:,0] = starting_assets

	#Based on the individual chunks, these are the households choices
	for s in range(1,S):
		c_path[:,s] = (beta * (1 + rpath_chunk[:,s] - delta))**(1/sigma) * c_path[:,s-1]
		asset_path[:,s] = wpath_chunk[:,s]*e_chunk[:,s-1] + (1 + rpath_chunk[:,s-1] - delta)*asset_path[:,s-1] - c_path[:,s-1]
	
	asset_path[:,s+1] = wpath_chunk[:,s]*e_chunk[:,s] + (1 + rpath_chunk[:,s] - delta)*asset_path[:,s] - c_path[:,s]

	#Returns the relevant part of c_path and asset_path for all countries 
	return c_path[:,0:S-current_s], asset_path[:,0:S+1-current_s]

def find_optimal_starting_consumptions(c_1, wpath_chunk, rpath_chunk, epath_chunk, starting_assets, current_s, params):
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
		-current_s: Current age of the agent


	Objects in Function:
		-cpath: Path of consumption based on chunk given.
		-assets_path: Path of assets based on the chunks given

	Outputs:
		-Euler:A flattened version of the assets_path matrix

	"""
	#Executes the get_household_choices_path function. Sees above.
	c_path, assets_path = get_lifetime_decisions(params, c_1, wpath_chunk, rpath_chunk, epath_chunk, starting_assets, current_s)

        if np.any(c_path<0):
            c_path=np.ones(I)*999.

	Euler = np.ravel(assets_path[:,-1])

	return Euler

def get_cons_assets_matrix(params, wpath, rpath, starting_assets):

	I, S, T, T_1, beta, sigma, delta, e, StartFertilityAge, StartDyingAge, N_matrix, MortalityRates = params

	#Initializes timepath variables
	c_timepath = np.zeros((I,S,S+T+1))
	a_timepath = np.zeros((I, S+1, S+T+1)) #I,S+1,S+T+1
	a_timepath[:,:,0]=starting_assets
	bq_timepath = np.zeros((I, S+1, S+T+1)) #Is this too big?

	c_timepath[:,S-1,0] = wpath[:,0]*e[:,S-1,0] + (1 + rpath[:,0] - delta)*a_timepath[:,S-1,0]

	#Fills the upper triangle
	for s in range(S-2,-1, -1):
		agent_assets = starting_assets[:,s]

		#We are only doing this for all generations alive in time t=0
		t = 0
		#We are iterating through each generation in time t=0
		current_s = s

		#Uses the previous generation's consumption at age s to get the value for our guess
		c_guess = c_timepath[:,s+1,t]/((beta*(1+rpath[:,t]-delta))**(1/sigma))

		#Gets optimal initial consumption beginning in the current age of the agent using chunks of w and r that span the lifetime of the given generation
		household_params = (I, S, beta, sigma, delta)

		opt_consump = opt.fsolve(find_optimal_starting_consumptions, c_guess, args = \
			(wpath[:,t:t+S], rpath[:,t:t+S], e[:,0,t:t+S],agent_assets, current_s, household_params))

		#Gets optimal timepaths beginning initial consumption and starting assets
		cpath_indiv, apath_indiv = get_lifetime_decisions\
			(household_params, opt_consump, wpath[:,t:t+S], rpath[:,t:t+S], e[:,0,t:t+S], agent_assets, current_s)

		for i in xrange(I):
			np.fill_diagonal(c_timepath[i,s:,:], cpath_indiv[i,:])
			np.fill_diagonal(a_timepath[i,s:,:], apath_indiv[i,:])

		bq_params = (I, S, T, StartFertilityAge, StartDyingAge, N_matrix[:,StartDyingAge:,s], N_matrix[:,StartFertilityAge:StartDyingAge+1,s], MortalityRates[:,StartDyingAge:,s])
		bq_timepath[:,:,S-s-2] = getBequests(bq_params, a_timepath[:,StartDyingAge:,S-s-2])

		#print np.round(cpath_indiv[0,:], decimals=3), opt_consump[0]
		#print np.round(np.transpose(c_timepath[0,:,:T_1-s+3]), decimals=3)
		#print np.round(starting_assets[0,:], decimals=3)
		#print np.round(assetpath_indiv[0,:], decimals=3), agent_assets[0]
		#print np.round(np.transpose(a_timepath[0,:,:T_1]), decimals=3)

	#Fills everything except for the upper triangle
	for t in xrange(1,T):
		current_s = 0 #This is always zero because this section deals with people who haven't been born yet in time T=0
		agent_assets = np.zeros((I))

		#Uses the previous generation's consumption at age s to get the value for our guess
		c_guess = c_timepath[:,s+1,t]/((beta*(1+rpath[:,t+1]-delta))**(1/sigma))

		optimalconsumption = opt.fsolve(find_optimal_starting_consumptions, c_guess, args = \
			(wpath[:,t:t+S], rpath[:,t:t+S], e[:,0,t:t+S], agent_assets, current_s, household_params))

		cpath_indiv, assetpath_indiv = get_lifetime_decisions\
			(household_params, optimalconsumption, wpath[:,t:t+S], rpath[:,t:t+S], e[:,0,t:t+S], agent_assets, current_s)

		for i in range(I):
			np.fill_diagonal(c_timepath[i,:,t:], cpath_indiv[i,:])
			np.fill_diagonal(a_timepath[i,:,t:], assetpath_indiv[i,:])

		if t >= T_1:
			temp_t = T_1
		else:
			temp_t = t
		bq_params = (I, S, T, StartFertilityAge, StartDyingAge, N_matrix[:,StartDyingAge:,temp_t+S-2], N_matrix[:,StartFertilityAge:StartDyingAge+1,temp_t+S-2], MortalityRates[:,StartDyingAge:,temp_t+S-2])
		bq_timepath[:,:,t+S-2] = getBequests(bq_params, a_timepath[:,StartDyingAge:,temp_t+S-2])

		#bq_timepath[:,:,t+S-2] = getBequests(a_timepath[:,:,t+S-2], t+S-2)

	return c_timepath, a_timepath

def get_wpathnew_rpathnew(params, wpath, rpath, starting_assets, k_ss, kf_ss, w_ss, r_ss):
	"""
	Description:
		Takes initial paths of wages and rental rates, gives the consumption path and the the wage and rental paths that are implied by that consumption path.

	Inputs:
		-w_path0[I, S+T+1]: initial w path
		-r_path0[I, S+T+1]: initial r path

	Objects in Function:
	Note that these vary in dimension depending on the loop.
		-current_s: The age of the cohort at time 0
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
		-ypath1[I, S+T+1]: timepath of assets implied from initial guess

	"""
	I, S, T, T_1, beta, sigma, delta, alpha, e, A, StartFertilityAge, StartDyingAge, N_matrix, MortalityRates = params
	ca_params = (I, S, T, T_1, beta, sigma, delta, e, StartFertilityAge, StartDyingAge, N_matrix, MortalityRates)
	c_timepath, a_timepath = get_cons_assets_matrix(ca_params, wpath, rpath, starting_assets)

	#Calculates the total amount of capital in each country
	Kpath=np.sum(a_timepath,axis=1)

	#Calculates Aggregate Consumption
	Cpath=np.sum(c_timepath,axis=1)

	#After time period T, the total capital stock and total consumption is forced to be the steady state
	Kpath[:,T:] = np.einsum("i,t->it", k_ss, np.ones(S+1))
	Cpath[:,T:] = np.einsum("i,t->it", Cpath[:,T-1], np.ones(S+1))

	#Gets the foriegned owned capital
	kf_params = (I, S, T, alpha, e, A)
	kfpath=get_foreignK_path(kf_params, Kpath, rpath, k_ss, kf_ss)

	#Based on the overall capital path and the foreign owned capital path, we get new w and r paths.
	prices_params = (S, T, alpha, e, A)
	wpath_new, rpath_new, Ypath = get_prices(prices_params, Kpath, kfpath, w_ss, r_ss)

	return wpath_new, rpath_new, Cpath, Kpath, Ypath

def get_Timepath(params, wstart, rstart, assets_init, k_ss, kf_ss, w_ss, r_ss):

    I, S, T, T_1, beta, sigma, delta, alpha, e, A, StartFertilityAge, StartDyingAge, N_matrix, MortalityRates, distance, diff, xi, MaxIters = params

    Iter=1 #Serves as the iteration counter
    wr_params = (I, S, T, T_1, beta, sigma, delta, alpha, e, A, StartFertilityAge, StartDyingAge, N_matrix, MortalityRates)

    while distance>diff and Iter<MaxIters: #The timepath iteration runs until the distance gets below a threshold or the iterations hit the maximum

            wpath_new, rpath_new, cpath_new, kpath_new, ypath_new = \
            get_wpathnew_rpathnew(wr_params, wstart,rstart, np.column_stack((np.zeros(I), assets_init, np.zeros(I))), k_ss, kf_ss, w_ss, r_ss)

            dist1=sp.linalg.norm(wstart-wpath_new,2) #Norm of the wage path
            dist2=sp.linalg.norm(rstart-rpath_new,2) #Norm of the intrest rate path
            distance=max([dist1,dist2]) #We take the maximum of the two norms to get the distance

            print "Iteration:",Iter,", Norm Distance: ", distance#, "Euler Error, ", EError
            Iter+=1 #Updates the iteration counter
            if distance<diff or Iter==MaxIters: #When the distance gets below the tolerance or the maximum of iterations is hit, then the TPI finishes.
                wend=wpath_new
                rend=rpath_new
                cend=cpath_new
                kend=kpath_new
                yend=ypath_new
            if Iter==MaxIters: #In case it never gets below the tolerance, it will throw this warning and give the last timepath.
                print "Doesn't converge within the maximum number of iterations"
                print "Providing the last iteration"

            wstart=wstart*xi+(1-xi)*wpath_new #Convex conjugate of the wage path
            rstart=rstart*xi+(1-xi)*rpath_new #Convex conjugate of the intrest rate path

    return wend, rend, cend, kend, yend

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
        Name="China"
    if Country=="Country 3":
        Name="Japan"
    if Country=="Country 4":
        Name="Korea"
    if Country=="Country 5":
        Name="Russia"
    if Country=="Country 6":
        Name="India"

    #Add More Country labels here

    return Name

def plotTimepaths(I, S, T, wpath, rpath, cpath, kpath, ypath, CountryNamesON):

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

    for i in xrange(I): #Rental Rates
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,T),rpath[i,:T], label=label1)
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
        plt.plot(np.arange(0,T),ypath[i,:T],label=label1)
    plt.title("Time path for Output")
    plt.ylabel("Output Stock level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()
	