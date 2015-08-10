from __future__ import division
import numpy as np
import scipy as sp
import scipy.optimize as opt
from matplotlib import pyplot as plt

#Parameters Zone
S = 20 #Upper bound of age for agents
Countries = 2 #Number of Countries
T = int(round(2.5*S)) #Number of time periods to convergence, based on Rick Evans' function.
beta = .95 #Future consumption discount rate
sigma = 1 #Leave it at 1, fsolve struggles with the first sigma we used (3).
delta = .1 #Depreciation Rate
alpha = .3 #Capital Share of production
e = np.ones((Countries, S, T+S+1)) #Labor productivities
A = np.ones(Countries) #Techonological Change, used for idential countries
#A=np.array([1.25,1.35,1,1.65,1.1]) #Techonological Change, used for when countries are different
diff=1e-6 #Convergence Tolerance
distance=10 #Used in taking the norm, arbitrarily set to 10
xi=.8 #Parameter used to take the convex conjugate of paths
MaxIters=300 #Maximum number of iterations on TPI.

#Program Levers

PrintSS=True #Prints the result of the Steady State functions
CalcTPI=True #Activates the calculation of Time Path Iteration
#NOTE: Graphing only works if CalcTPI is activated.
Graphs=True #Activates graphing the graphs.
CountryNamesON=True #Turns on labels for the graphs. Replaces "Country x" with proper names.

#BEGIN NEW CODE
T_1 = S #This is like TransYear in the FORTRAN I think

if S > 50:
	T_1 = 50

StartFertilityAge = int(S/80.*23)#The age when agents have their first children
EndFertilityAge = int(S/80.*45)#The age when agents have their last children
StartDyingAge = int(S/80.*68)#The first age agents can begin to die
MaxImmigrantAge = int(S/80.*65)#All immigrants are between ages 0 and MaxImmigrantAge
g_A = 0.001#Technical growth rate

print "T =", T
print "T_1", T_1
print "StartFertilityAge", StartFertilityAge
print "EndFertilityAge", EndFertilityAge
print "StartDyingAge", StartDyingAge
print "MaxImmigrantAge", MaxImmigrantAge

#DEMOGRAPHICS FUNCTIONS

def getDemographics():
	"""
	Description:
		-Imports data from csv files for initial populations, fertility rates, mortality rates, and net migrants. 
		-Stores these data sets in their respective matrices, and calculates population distribuitons through year T.
		-NOTE: FOR NOW THIS FUNCTION ONLY USES DATA FOR THE USA. NEEDS TO EVENTUALLY ADD MORE COUNTRIES
	Inputs:
		-None, but uses the global variables T, T_1, StartFertilityAge, EndFertilityAge, StartDyingAge, and MaxImmigrantAge
	Output:
		-FertilityRates: Numpy array that contains fertilty rates for all countries, ages, and years
		-MortalityRates: Numpy array that contains mortality rates for all countries, ages, and years
		-Migrants: Numpy array that contains net migration for all countries and ages
		-N_matrix: Numpy array that contains population numbers for all countries, ages, and years
	"""

	#Imports and scales data for the USA. Imports a certain number of generations according to the value of S
	USAPopdata = np.loadtxt(("Data_Files/population.csv"),delimiter=',',skiprows=1, usecols=[1])[:S+1]*1000
	USAFertdata = np.loadtxt(("Data_Files/usa_fertility.csv"),delimiter=',',skiprows=1, usecols=range(1,EndFertilityAge+2-StartFertilityAge))[48:48+T_1,:]
	USAMortdata = np.loadtxt(("Data_Files/usa_mortality.csv"),delimiter=',',skiprows=1, usecols=range(1,S+1-StartDyingAge))[:T_1,:]
	USAMigdata = np.loadtxt(("Data_Files/net_migration.csv"),delimiter=',',skiprows=1, usecols=[1])[:MaxImmigrantAge]*100

	#Initializes demographics matrices
	N_matrix = np.zeros((Countries, S+1, T))
	Nhat_matrix = np.zeros((Countries, S+1, T))
	FertilityRates = np.zeros((Countries, S+1, T))
	MortalityRates = np.zeros((Countries, S+1, T))
	Migrants = np.zeros((Countries, S+1, T))
	g_N = np.zeros(T)

	#NOTE: For now we set fertility, mortality, number of migrants, and initial population the same for all countries. 
	
	#Sets initial total population (N_matrix), percentage of total world population (Nhat_matrix)
	N_matrix[:,:,0] = np.tile(USAPopdata, (Countries, 1))
	Nhat_matrix[:,:,0] = N_matrix[:,:,0]/np.sum(N_matrix[:,:,0])

	#Fertility Will be equal to 0 for all ages that don't bear children
	FertilityRates[:,StartFertilityAge:EndFertilityAge+1,:T_1] = np.einsum("ts,it->ist", USAFertdata, np.ones((Countries,T_1)))

	#Mortality be equal to 0 for all young people who aren't old enough to die
	MortalityRates[:,StartDyingAge:-1,:T_1] = np.einsum("ts,it->ist", USAMortdata, np.ones((Countries,T_1)))

	#The last generation dies with probability 1
	MortalityRates[:,-1,:] = np.ones((Countries, T))

	#The number of migrants is the same for each year
	Migrants[:,:MaxImmigrantAge,:T_1] = np.einsum("s,it->ist", USAMigdata, np.ones((Countries,T_1)))

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
		N_matrix[:,0,t] = np.sum((N_matrix[:,:,t-1]*FertilityRates[:,:,t-1]),axis=1)
		Nhat_matrix[:,0,t] = np.exp(-g_N[t-1])*np.sum((Nhat_matrix[:,:,t-1]*FertilityRates[:,:,t-1]),axis=1)

		#Finds the immigration rate for each year
		ImmigrationRate = Migrants[:,:,t-1]/N_matrix[:,:,t-1]

		#Gets the population and percentage of population for the next year, taking into account immigration and mortality
		N_matrix[:,1:,t] = N_matrix[:,:-1,t-1]*(1+ImmigrationRate[:,:-1]-MortalityRates[:,:-1,t-1])
		Nhat_matrix[:,1:,t] = np.exp(-g_N[t-1])*Nhat_matrix[:,:-1,t-1]*(1+ImmigrationRate[:,:-1]-MortalityRates[:,:-1,t-1])
		#print "For t =", t, "\ng_N[t] = np.sum(\n", Nhat_matrix[:,:,t], "\n*(\n", FertilityRates[:,:,t], "\n+\n", ImmigrationRate, "\n-\n", MortalityRates[:,:,t],"\n))"
		
		#Gets the growth rate for the next year
		g_N[t] = np.sum(Nhat_matrix[:,:,t]*(FertilityRates[:,:,t] + ImmigrationRate - MortalityRates[:,:,t]), axis=(0, 1))
		#print np.sum(Nhat_matrix[:,:,t])#This should be equal to 1.0

	#Gets labor endowment per household. For now it grows at a constant rate g_A
	l_endowment = np.cumsum(np.ones(T)*g_A)

	return FertilityRates, MortalityRates, Migrants, N_matrix, Nhat_matrix

def plotDemographics(index, years, name):
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

def getBequests(assets, current_t):
	"""
	Description:
		-Gets the value of the bequests given to each generation
	Inputs:
		-assets: Assets for each generation in a given year
		-current_t: Integer that indicates the current year. Used to pull information from demographics global matrices like FertilityRates
	Output:
		-bq: Numpy array that contains the number of bequests for each generation in each country.
	"""

	#Initializes bequests
	bq = np.zeros((Countries, S+1))

	#Gets the total assets of the people who died this year
	BQ = np.sum(assets[:,StartDyingAge:]*MortalityRates[:,StartDyingAge:,current_t]*N_matrix[:,StartDyingAge:,current_t], axis=1)

	#Distributes the total assets equally among the eligible population for each country
	#NOTE: This will likely change as we get a more complex function for distributing the bequests
	num_bequest_receivers = np.sum(N_matrix[:,StartFertilityAge:StartDyingAge+1,current_t], axis=1)
	bq_Distribution = BQ/num_bequest_receivers
	bq[:,StartFertilityAge:StartDyingAge+1] = np.einsum("i,s->is", bq_Distribution, np.ones(StartDyingAge+1-StartFertilityAge))

	return bq

FertilityRates, MortalityRates, Migrants, N_matrix, Nhat_matrix = getDemographics()

#plotDemographics(0,[0,19],"USA")

def hatvariables(Kpathreal, kfpathreal, Nhat_matrix):

	#THIS FUNCTION HAS EQUATIONS 2.13-2.16 AND 2.19-2.20, BUT STILL NEEDS TO BE INCORPORATED INTO THE REST OF THE MODEL TO COMPLETELY TEST

	#We are only using up until T periods rather than T+S+1 since Nhat only goes out to T
	Kpath = Kpathreal[:,:T]
	kfpath = kfpathreal[:,:T]
	temp_e = np.ones((Countries, S+1, T))#THIS SHOULD ONLY BE UNTIL WE GET S GENERATIONS RATHER THAN S-1

	n = np.sum(temp_e[:,:,:T]*Nhat_matrix, axis=1)
	ypath = (Kpath**alpha) * (np.einsum("i,it->it", A, n)**(1-alpha))
	rpath = alpha * ypath / Kpath
	wpath = (1-alpha) * ypath / n
	"""
	#NOTE:This goes in the get_householdchoices_path function

	c_path = np.zeros((Countries, S))
	asset_path = np.zeros((Countries, S+1))

	c_path[:,0] = c_1
	asset_path[:,0] = starting_assets

	for s in range(1,S):
		c_path[:,s] = ((beta * (1 + rpath_chunk[:,s] - delta))**(1/sigma) * c_path[:,s-1])/np.exp(g_A)
		asset_path[:,s] = (wpath_chunk[:,s]*e[:,0,s-1] + (1 + rpath_chunk[:,s-1] - delta)*asset_path[:,s-1] + bq_chunk - c_path[:,s-1])/np.exp(g_A)

	asset_path[:,s+1] = wpath_chunk[:,s]*e_chunk[:,s] + (1 + rpath_chunk[:,s] - delta)*asset_path[:,s] - c_path[:,s]

	"""

def get_hh_Rick_way(c_guess, w_chunk, r_chunk, e_chunk, agent_assets, current_s):

	print "Hello"


#STEADY STATE FUNCTIONS

def getOtherVariables(assets, kf):
<<<<<<< HEAD
    	"""
    	Description:
        	-Based on the assets and capital held by foreigners, we calculate the other variables.
    	Inputs:
        	-assets [Countries,S+1]: Matrix of assets
        	-kf[Countries, ]: Domestic capital held by foreigners
        Objects in function:
                -NONE that aren't Already listed

    	Output:
        	-k[Countries,]: Capital (1.10)
        	-n[Countries,]: Sum of labor productivities (1.11)
        	-y[Countries,]: Output (1.12)
        	-r[Countries,]: Rental Rate (1.13)
        	-w[Countries,]: Wage (1.14)
        	-c_vec[Countries,S]: Vector of consumptions (1.15)
    	"""
	
	k = np.sum(assets[:,1:-1], axis=1) - kf
	n = np.sum(e[:,:,0], axis=1)
	y = (k**alpha) * ((A*n)**(1-alpha))
	r = alpha * y / k
	w = (1-alpha) * y / n

	c_vec = np.einsum("i, is -> is", w, e[:,:,0]) + np.einsum("i, is -> is",(1 + r - delta) \
	, assets[:,:-1]) - assets[:,1:]

	return k, n, y, r, w, c_vec

def SteadyStateSolution(guess):
	"""
    Description: 
	    -This is the function that will be optimized by fsolve.
	Inputs:
   	    -guess[Countries,S+1]: vector that pieced together from assets and kf.
	Objects in Function:
	-kf[Countries,]:
            -assets[Countries,S]:
            -k[Countries,]:
            -n[Countries,]:
            -y[Countries,]:
            -r[Countries,]:
            -w[Countries,]:
            -c_vec[Countries, S]:
            -Euler_c:
            -Euler_r:
            -Euler_kf:
        Output:
            -all_Euler: Similar to guess, it's a vector that's has both assets and kf.    	
	"""
	#Takes a 1D guess of length Countries*S and reshapes it to match what the original input into the fsolve looked like since fsolve flattens numpy arrays
	guess = np.reshape(guess[:,np.newaxis], (Countries, S))

        #Sets kf as the last element of the guess vector for each country and assets as everything else
	assets = guess[:,:-1]
	kf = guess[:,-1]

	#You have 0 assets when you're born, and 0 when you die
	assets = np.column_stack((np.zeros(Countries), assets, np.zeros(Countries)))

	#Based on the assets and kf, we get the other vectors
	k, n, y, r, w, c_vec = getOtherVariables(assets, kf)

	#Gets Euler equations
	Euler_c = c_vec[:,:-1] ** (-sigma) - beta * c_vec[:,1:] ** (-sigma) * (1 + r[0] - delta)
	Euler_r = r[1:] - r[0]
	Euler_kf = np.sum(kf)

        print "Euler C", Euler_c.shape
        print "Euler r",Euler_r.shape
        print "Euler kf",Euler_kf.shape

	#Makes a new 1D vector of length Countries*S that contains all the Euler equations
	all_Euler = np.append(np.append(np.ravel(Euler_c), np.ravel(Euler_r)), Euler_kf)

        print "all_Euler",all_Euler.shape

	return all_Euler

def getSteadyState(assets_init, kf_init):
	"""
	Description:
        This takes the initial guess for assets and kf. Since the function
	    returns a matrix, this unpacks the individual parts.
	Inputs:
	    -assets_init:Intial guess for asset path
	    -kf_init:Initial guess on foreigner held capital  

        Objects in Function:
            -guess:
            -ss:

	Outputs:
	    -assets_ss:Calculated assets steady state
	    -kf_ss:Calculated foreign capital
	"""
        #Merges the assets and kf together into one matrix that can be inputted into the fsolve function
	guess = np.column_stack((assets_init, kf_init))

        #Solves for the steady state
	ss = opt.fsolve(SteadyStateSolution, guess)

	print "\nSteady State Found!\n"
	ss = np.array(np.split(ss, Countries))
        #Breaks down the steady state into the two separate assets and kf matrices.
	assets_ss = ss[:,:-1]
	kf_ss = ss[:,-1]

        return assets_ss, kf_ss

#TIMEPATH FUNCTIONS

def get_householdchoices_path(c_1, wpath_chunk, rpath_chunk, e_chunk, starting_assets, current_s):
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
	    -c_path[Countries, S]: Path of consumption until the agent dies
	    -asset_path[Countries, S+1]: Path of assets until the agent dies
	"""
        #Initializes the cpath and asset path vectors
	c_path = np.zeros((Countries, S))
	asset_path = np.zeros((Countries, S+1))
        #print c_1.shape
        #For each country, the cpath and asset path vectors' are the initial values provided.
	c_path[:,0] = c_1
	asset_path[:,0] = starting_assets

        #Based on the individual chunks, these are the households choices
	for s in range(1,S):
		c_path[:,s] = (beta * (1 + rpath_chunk[:,s] - delta))**(1/sigma) * c_path[:,s-1]
                asset_path[:,s] = wpath_chunk[:,s]*e[:,s,s-1] + (1 + rpath_chunk[:,s-1] - delta)*asset_path[:,s-1] - c_path[:,s-1]
	
        asset_path[:,s+1] = wpath_chunk[:,s]*e_chunk[:,s] + (1 + rpath_chunk[:,s] - delta)*asset_path[:,s] - c_path[:,s]

	#Returns the relevant part of c_path and asset_path for all countries 

	return c_path[:,0:S-current_s], asset_path[:,0:S+1-current_s]

def find_optimal_starting_consumptions(c_1, wpath_chunk, rpath_chunk, epath_chunk, starting_assets, current_s):
	"""
	Description:
		Takes the assets path from the get_householdchoices_path function and creates  
	Inputs:
		c_1: Initial consumption (not necessarily for the year they were born)
		wpath_chunk: Wages of an agents lifetime, a part of the timepath
		rpath_chunk: Rental rate of an agents lifetime, another part of the timeparth.
		epath_chunk: Worker productivities of an agents lifetime, another part.
		starting_assets: Initial assets of the agent. It's 0 at the beginning of life.
		current_s: Current age of the agent

        Objects in Function:
            -cpath:
            -assets_path:

        Outputs:
            -Euler:A flattened version of the assets_path matrix

	"""

	#Executes the get_household_choices_path function. Sees above.
	c_path, assets_path = get_householdchoices_path(c_1, wpath_chunk, rpath_chunk, epath_chunk, starting_assets, current_s)

	Euler = np.ravel(assets_path[:,-1])

	return Euler

def get_prices(Kpath, kf_tpath):
    	"""
        Description:
            Based on the given paths, the paths for wages and rental rates are figured
            out based on equations 1.4-1.5

        Inputs:
            -assets_tpath: Asset timepath
            -kf_tpath: Foreign held capital timepath.

        Objects in Function:
            -ypath:
            -Kdpath:

        Outputs:
            -wpath: Wage path
            -rpath: Rental rate path

    	"""
<<<<<<< HEAD
	#Gets non-price variables needed to caluclate prices
	n = np.sum(e, axis=1)
	ypath=np.zeros((Countries,S+T+1))
	for i in xrange(Countries):
		ypath[i,:] = (Kpath[i,:]**alpha) * ((A[i]*n[i,:])**(1-alpha))
=======

        #Initializes the path for output, y.
        ypath=np.zeros((Countries,S+T+1))

        Kdpath=Kpath-kf_tpath

	#Gets non-price variables needed to caluclate prices
        n = np.sum(e, axis=1) #Sum of the labor productivities
        for i in xrange(Countries):
            ypath[i,:] = (Kpath[i,:]**alpha) * ((A[i]*n[i,:])**(1-alpha))

>>>>>>> upstream/master
	#Gets prices
	rpath = alpha * ypath / Kpath
	wpath = (1-alpha) * ypath / n

<<<<<<< HEAD
	#Stacks extra year of steady state for later code that tries to use outside index
	#rpath = np.column_stack((rpath, rpath[:,-1]))
	#wpath = np.column_stack((wpath, wpath[:,-1]))

	#Returns only the first country's interest rate, since they should be the same in theory
	return wpath, rpath
=======
        for i in xrange(Countries):
            rpath[i,T:]=r_ss[i]
            wpath[i,T:]=w_ss[i]

	#Returns only the first country's interest rate, since they should be the same in theory
	return wpath, rpath, ypath
>>>>>>> upstream/master

def get_foreignK_path(Kpath,rpath):
        """
        Description:
           This calculates the timepath of the foreign capital stock. This is based on equation (1.12 and 1.13).
<<<<<<< HEAD
        Inputs:
            apath-Asset path, from our calculations
            rpath-Rental Rate path, also from our calculation
        Outputs:
            kfPath-Path of domestic capital held by foreigners.
        """
        #Sums the labor productivities across cohorts
        n = np.sum(e[:,:,:], axis=1)
        #Sums the assets across cohorts to give the domestic capital stock

        #Declares the array that will later be used.
        kfPath=np.zeros((Countries,S+T+1))
        kDPath=np.zeros((Countries,S+T+1))

        for i in xrange(1,Countries):
            kDPath[i,:]=(rpath[i,:]/alpha)**(1/(alpha-1))*n[i,:]*A[i]

        #This is using equation 1.13 solved for the foreign capital stock to caluclate the foreign capital stock
        kfPath=Kpath-kDPath

        #To satisfy 1.18, the first country's assets is the negative of the sum of all the other countries' assets
        kfPath[0,:]=-np.sum(kfPath,axis=0)
        
        return kfPath

def get_wpath1_rpath1(w_path0, r_path0, starting_assets):
	"""
	Description:
		Takes initial paths of wages and rental rates, gives the consumption path and the the wage and rental paths that are implied by that consumption path.

	Inputs:
		w_path0: initial w path
		r_path0: initial r path:
	Outputs:
		w_path1: calculated w path
		r_path1: calculated r path
		c_timepath: timepath of consumption implied from initial guess
		assets_timepath: timepath of assets implied from initial guess

	"""
	#Initializes timepath variables
	c_timepath = np.zeros((Countries,S,T+1))
	assets_timepath = np.zeros((Countries,S+1, S+T+1))

	#NEW CODE
	bequests_timepath = np.zeros((Countries, S+1, T+1))
	bequests_timepath[:,:,0] = getBequests(np.tile(np.arange(S+1)/S, (Countries,1)), 0)
	assets_timepath[:,:,0] = starting_assets
=======
        Inputs:
            apath:Asset path, from our calculations
            rpath:Rental Rate path, also from our calculation
        
        Objects in Function:
            kDpath:
            n:
            kf_ss:
            A:

        Outputs:
            kfPath-Path of domestic capital held by foreigners.
        """
        #Sums the labor productivities across cohorts
        n = np.sum(e, axis=1)
        #Sums the assets across cohorts to give the domestic capital stock

        #Declares the array that will later be used.
        kfPath=np.zeros((Countries,S+T+1))
        kDPath=np.zeros((Countries,S+T+1))

        #Goes through each country and individually calculates the domestic owned capital stock.
        for i in xrange(1,Countries):
            kDPath[i,:]=(rpath[i,:]/alpha)**(1/(alpha-1))*n[i,:]*A[i]

        #This is using equation 1.13 solved for the foreign capital stock to caluclate the foreign capital stock
        kfPath=Kpath-kDPath

        #To satisfy 1.18, the first country's assets is the negative of the sum of all the other countries' assets
        kfPath[0,:]=-np.sum(kfPath,axis=0)

        for i in xrange(Countries):
            kfPath[i,T:]=kf_ss[i]

        #print kfPath[0,:]
        
        return kfPath

def get_wpath1_rpath1(w_path0, r_path0, starting_assets):
        """
        Description:
            Takes initial paths of wages and rental rates, gives the consumption path and the the wage and rental paths that are implied by that consumption path.

        Inputs:
            -w_path0: initial w path
            -r_path0: initial r path
        
        Objects in Function:
            -current_s:
            -opt_consump:
            -starting_assets:
            -cpath_indiv:
            -assetpath_indiv:
            -optimalconsumption:
            -C_timepath:
            -assets_timepath:
            -kfpath:
            -agent assets:


        Outputs:
            -w_path1: calculated w path
            -r_path1: calculated r path
            -CPath:
            -Kpath:
            -ypath1: timepath of assets implied from initial guess

        """
        #Initializes timepath variables
	c_timepath = np.zeros((Countries,S,T+1))
	#test = np.zeros((Countries,S,T+1))
	assets_timepath = np.zeros((Countries,S+1, S+T+1)) #Countries,S+1,S+T+1
        assets_timepath[:,:,0]=starting_assets

        c_timepath[:,S-1,0]=w_path0[:,0]*e[:,S-1,0]+(1+r_path0[:,0]-delta)*assets_timepath[:,S-1,0]

	#Makes an initial guess for fsolve. Will be used in each generation
	#c_guess = np.ones(Countries)*.9
>>>>>>> upstream/master

	c_timepath[:,S-1,0] = w_path0[:,0]*e[:,S-1,0] + (1 + r_path0[:,0] - delta)*assets_timepath[:,S-1,0]
	#c_guess = np.ones(Countries)*.9
	#Fills the upper triangle
<<<<<<< HEAD
	for s in range(S-2,-1, -1):
		agent_assets = starting_assets[:,s]
=======
	for s in xrange(S-2, -1, -1):
                #c_guess=np.diagonal(c_timepath[:,S-s:,:s-1])
                #print c_guess.shape
                #starting_assets=np.ones((Countries))*.1*s
                agent_assets=starting_assets[:,s]

>>>>>>> upstream/master
		#We are only doing this for all generations alive in time t=0
		t = 0
		#We are iterating through each generation in time t=0
		current_s = s
		#test = np.arange(Countries*S*(T+1)).reshape(Countries,S,T+1)

		#NOTE: TELL JEFF THAT WE ADDED THE SUBTRACT AT THE END

		c_guess = c_timepath[:,s+1,t]/((beta*(1+r_path0[:,t]-delta))**(1/sigma))
		c_guessnew= np.diagonal(c_timepath[:, s+1:, :], axis1=1, axis2=2)
		c_guessnewfirst = c_timepath[:,s+1,t]/((beta*(1+r_path0[:,t]-delta))**(1/sigma))
		#print "New guess", np.round(c_guessnew, decimals=3), c_guessnewfirst[0]
		#print np.round(np.transpose(c_timepath[0,:,:T_1-s+3]), decimals=3)

		c_guessnew = np.column_stack((c_guessnewfirst,c_guessnew))


                c_guess=c_timepath[:,s+1,t]/(beta*(1+r_path0[:,t+1]-delta))**(1/sigma)

		#Gets optimal initial consumption beginning in the current age of the agent using chunks of w and r that span the lifetime of the given generation
<<<<<<< HEAD
		opt_consump = opt.fsolve(find_optimal_starting_consumptions, c_guess, args = (w_path0[:,t:t+S], r_path0[:,t:t+S], e[:,0,t:t+S], agent_assets, current_s))

		#Gets optimal timepaths beginning initial consumption and starting assets
		cpath_indiv, assetpath_indiv = get_householdchoices_path(opt_consump, w_path0[:,t:t+S], r_path0[:,t:t+S], e[:,0,t:t+S], agent_assets, current_s)

		for i in range(Countries):
			np.fill_diagonal(c_timepath[i,s:,:], cpath_indiv[i,:])
			np.fill_diagonal(assets_timepath[i,s:,:], assetpath_indiv[i,:])

		#NEW CODE
		bequests_timepath[:,:,S-s-1] = getBequests(assets_timepath[:,:,S-s-1], s)

		
		#print np.round(cpath_indiv[0,:], decimals=3), opt_consump[0]
		#print np.round(np.transpose(c_timepath[0,:,:T_1-s+3]), decimals=3)
		#print np.round(starting_assets[0,:], decimals=3)
		#print np.round(assetpath_indiv[0,:], decimals=3), agent_assets[0]
		#print np.round(np.transpose(assets_timepath[0,:,:T_1]), decimals=3)
		

	#Fills everything except for the upper triangle
	for t in range(1,T):
		current_s = 0
		starting_assets = np.zeros((Countries))

		c_guess = c_timepath[:,s+1,t]/((beta*(1+r_path0[:,t]-delta))**(1/sigma))

		optimalconsumption = opt.fsolve(find_optimal_starting_consumptions, c_guess, args = (w_path0[:,t:t+S], r_path0[:,t:t+S], e[:,0,t:t+S], starting_assets, current_s))
=======
		
		opt_consump_1 = opt.fsolve(find_optimal_starting_consumptions, c_guess, args = \
                        (w_path0[:,t:t+S], r_path0[:,t:t+S], e[:,s,t:t+S],agent_assets, current_s))
		
		#Gets optimal timepaths beginning initial consumption and starting assets
		cpath_indiv, assetpath_indiv = get_householdchoices_path\
                        (opt_consump_1, w_path0[:,t:t+S], r_path0[:,t:t+S], e[:,s,t:t+S], agent_assets, current_s)

		for i in xrange(Countries):
			np.fill_diagonal(c_timepath[i,s:,:], cpath_indiv[i,:])
			np.fill_diagonal(assets_timepath[i,s:,:], assetpath_indiv[i,:])

	#Fills everything except for the upper triangle
	for t in xrange(1,T):
		current_s = 0 #This is always zero because this section deals with people who haven't been born yet in time T=0
		agent_assets = np.zeros((Countries))
                c_guess=c_timepath[:,s+1,t]/(beta*(1+r_path0[:,t+1]-delta))**(1/sigma)
>>>>>>> upstream/master

		optimalconsumption = opt.fsolve(find_optimal_starting_consumptions, c_guess, args = \
                        (w_path0[:,t:t+S], r_path0[:,t:t+S], e[:,s,t:t+S], agent_assets, current_s))

<<<<<<< HEAD
		for i in range(Countries):
			np.fill_diagonal(c_timepath[i,:,t:], cpath_indiv[i,:])
			np.fill_diagonal(assets_timepath[i,:,t:], assetpath_indiv[i,:])

		#NEW CODE
		bequests_timepath[:,:,t] = getBequests(assets_timepath[:,:,t], t)

	Kpath=np.sum(assets_timepath,axis=1)
	for i in xrange(Countries):
		Kpath[i,T:]=k_ss[i]

	kfpath=get_foreignK_path(Kpath, r_path0)
	w_path1, r_path1 = get_prices(Kpath,kfpath)

	return w_path1, r_path1, c_timepath, assets_timepath

def CountryLabel(Country):
=======
		cpath_indiv, assetpath_indiv = get_householdchoices_path\
                        (optimalconsumption, w_path0[:,t:t+S], r_path0[:,t:t+S], e[:,s,t:t+S], agent_assets, current_s)

                for i in xrange(Countries):
			np.fill_diagonal(c_timepath[i,:,t:], cpath_indiv[i,:])
			np.fill_diagonal(assets_timepath[i,:,t:], assetpath_indiv[i,:])

        #Calculates the total amount of capital in each country
        Kpath=np.sum(assets_timepath,axis=1)
        #After time period T, the total capital stock is forced to be the steady state
        for i in xrange(Countries):
            Kpath[i,T:]=k_ss[i]

        #Calculates Aggregate Consumption
        Cpath=np.sum(c_timepath,axis=1)
        for i in xrange(Countries):
            Cpath[i,T:]=Cpath[i,T-1]

        #Gets the foriegned owned capital
        kfpath=get_foreignK_path(Kpath, r_path0)

        #Based on the overall capital path and the foreign owned capital path, we get new w and r paths.
 	w_path1, r_path1, ypath1 = get_prices(Kpath,kfpath)

        #kdpath=Kpath-kfpath

        return w_path1, r_path1, Cpath, Kpath, ypath1

def CountryLabel(Country): #Activated by line 28
>>>>>>> upstream/master
    '''
    Description: 
        Converts the generic country label given for the graphs and converts it to a proper name

    Inputs:
        -Country: This is simply the generic country label

<<<<<<< HEAD
=======
    Objects in Function:
        -NONE


>>>>>>> upstream/master
    Outputs:
        -Name: The proper name of the country which you decide. Make sure the number of country names lines
            up with the number of countries, otherwise, the function will not proceed.
    '''
<<<<<<< HEAD
=======
    #Each country is given a number
>>>>>>> upstream/master
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

<<<<<<< HEAD
    #Add Country labels here
=======
    #Add More Country labels here
>>>>>>> upstream/master

    return Name
	
#MAIN CODE

#Initalizes initial guesses
assets_guess = np.ones((Countries, S-1))*.15
kf_guess = np.zeros((Countries))
<<<<<<< HEAD
w_initguess = np.zeros((Countries, T+S+1)) #T+S+1
r_initguess = np.ones((Countries, T+S+1))*.5 #T+S+1
=======
w_initguess = np.zeros((Countries, T+S+1)) 
r_initguess = np.ones((Countries, T+S+1))*.5 
>>>>>>> upstream/master

#Gets the steady state variables
assets_ss, kf_ss = getSteadyState(assets_guess, kf_guess)
k_ss, n_ss, y_ss, r_ss, w_ss, c_vec_ss = getOtherVariables(np.column_stack((np.zeros(Countries), assets_ss, np.zeros(Countries))), kf_ss)


<<<<<<< HEAD
if PrintSS==True:
=======
if PrintSS==True: #Prints the results of the steady state, line 23 activates this
>>>>>>> upstream/master
        print "assets steady state", assets_ss
        print "kf steady state", kf_ss
        print "k steady state", k_ss
        print "n steady state",n_ss
        print "y steady state", y_ss
        print "r steady state",r_ss
        print "w steady state", w_ss
        print "c_vec_ss steady state",c_vec_ss


<<<<<<< HEAD
#Sets initial assets and kf
assets_init = assets_ss/2
kf_init = kf_ss/2
=======
#Sets initial assets and kf, start with something close to the steady state
assets_init = assets_ss*.95
kf_init = kf_ss*.95
>>>>>>> upstream/master

#Gets initial k, n, y, r, w, and c
k_init, n_init, y_init, r_init, w_init, c_vec_init = getOtherVariables(np.column_stack((np.zeros(Countries), assets_init, np.zeros(Countries))), kf_init)

<<<<<<< HEAD
#NEW CODE
GAMMA1 = c_vec_init

#Gets initial guess for w and r paths
=======
#Gets initial guess for w and r paths. This is set up to be linear.
>>>>>>> upstream/master
for c in range(Countries):
	w_initguess[c, :T+1] = np.linspace(w_init[c], w_ss[c], T+1)
	r_initguess[c, :T+1] = np.linspace(r_init[c], r_ss[c], T+1)
	w_initguess[c,T+1:] = w_initguess[c,T]
	r_initguess[c,T+1:] = r_initguess[c,T]

<<<<<<< HEAD
#Gets new paths for wages and rental rates to see how close it is to the original code
wstart=w_initguess
#print "wstart dimensions",wstart.shape
rstart=r_initguess
#print "rstart dimensions",rstart.shape
Iter=0
dist1=10
dist2=10

if CalcTPI==True:
    while distance>diff and Iter<MaxIters:

    	    print "Iter",Iter

            wpath0,rpath0, cpath0, apath0=get_wpath1_rpath1(wstart,rstart, np.column_stack((np.zeros(Countries), assets_init, np.zeros(Countries))))

            dist1=sp.linalg.norm(wstart[:,:T]-wpath0[:,:T],2)
            dist2=sp.linalg.norm(rstart[:,:T]-rpath0[:,:T],2)
            distance=(dist1+dist2)/2

            wstart=wstart*xi+(1-xi)*wpath0
            rstart=rstart*xi+(1-xi)*rpath0

            Iter+=1
            if distance<diff or Iter==MaxIters:
                wend=wstart
                rend=rstart


if GraphWPath==True and CalcTPI==True:
    for i in xrange(Countries):
=======

#Gets new paths for wages and rental rates to see how close it is to the original code
wstart=w_initguess
rstart=r_initguess

Iter=1 #Serves as the iteration counter
if CalcTPI==True: #Time Path Iteration, activated by line 24
    while distance>diff and Iter<MaxIters: #The timepath iteration runs until the distance gets below a threshold or the iterations hit the maximum
            wpath0,rpath0, cpath0, kpath0, ypath0=get_wpath1_rpath1(wstart,rstart,np.column_stack((np.zeros(Countries),assets_init,np.zeros(Countries))))
            dist1=sp.linalg.norm(wstart-wpath0,2) #Norm of the wage path
            dist2=sp.linalg.norm(rstart-rpath0,2) #Norm of the intrest rate path
            distlist=[dist1,dist2] #Lists the two norms taken
            distance=max(distlist) #We take the maximum of the two norms to get the distance

            #print distance

            wstart=wstart*xi+(1-xi)*wpath0 #Convex conjugate of the wage path
            rstart=rstart*xi+(1-xi)*rpath0 #Convex conjugate of the intrest rate path

            print "Iteration:",Iter,", Norm Distance: ", distance#, "Euler Error, ", EError
            Iter+=1 #Updates the iteration counter
            if distance<diff or Iter==MaxIters: #When the distance gets below the tolerance or the maximum of iterations is hit, then the TPI finishes.
                wend=wstart
                rend=rstart
            if Iter==MaxIters: #In case it never gets below the tolerance, it will throw this warning and give the last timepath.
                print "Doesn't converge within the maximum number of iterations"
                print "Providing the last iteration"


#GRAPHING ZONE

if Graphs==True and CalcTPI==True:
    for i in xrange(Countries): #Wages
>>>>>>> upstream/master
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,T),wend[i,:T], label=label1)
    plt.title("Time path for Wages")
    plt.ylabel("Wages")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()

<<<<<<< HEAD
if GraphRPath==True and CalcTPI==True:
    for i in xrange(Countries):
=======
    for i in xrange(Countries): #Rental Rates
>>>>>>> upstream/master
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,T),rend[i,:T], label=label1)
    plt.title("Time path for Rental Rates")
    plt.ylabel("Rental Rates")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()

<<<<<<< HEAD

if GraphConsPath==True and CalcTPI==True:
    for i in xrange(Countries):
        plt.plot(np.sum(cpath0[i,:,:], axis=0))
    plt.title("Time Path for Aggregate Consumption")
    plt.show()
    
=======
    for i in xrange(Countries): #Aggregate Consumption
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,T+1),cpath0[i,:],label=label1)
    plt.title("Time Path for Aggregate Consumption")
    plt.ylabel("Consumption Level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()

    for i in xrange(Countries): #Aggregate Capital Stock
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,T),kpath0[i,:T],label=label1)
    plt.title("Time path for Capital Path")
    plt.ylabel("Capital Stock level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()

    for i in xrange(Countries):
        label1='Country '+str(i)
        if CountryNamesON==True:
            label1=CountryLabel(label1)
        plt.plot(np.arange(0,T),ypath0[i,:T],label=label1)
    plt.title("Time path for Output")
    plt.ylabel("Output Stock level")
    plt.xlabel("Time Period")
    plt.legend(loc="upper right")
    plt.show()


>>>>>>> upstream/master
