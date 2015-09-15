import numpy as np
import StepbyStepv1 as Stepfuncs
#import RobustChecker as rc

def TheWholeSmack(S,I,sigma):
	#Parameters Zone
	#I = 10 #Number of countries
	#S = 80 #Upper bound of age for agents
	#I=rc.I
	#S=rc.S
	T = int(round(2.5*S)) #Number of time periods to convergence, based on Rick Evans' function.

	T_1 = S #This is like TransYear in the FORTRAN I think
	if S > 50:
		T_1 = 50
	StartFertilityAge = int(S/80.*23)#The age when agents have their first children
	EndFertilityAge = int(S/80.*45)#The age when agents have their last children
	StartDyingAge = int(S/80.*68)#The first age agents can begin to die
	MaxImmigrantAge = int(S/80.*65)#All immigrants are between ages 0 and MaxImmigrantAge
	g_A = 0.001#Technical growth rate

	beta_ann=.95 #Starting future consumption discount rate
	delta_ann=.08 #Starting depreciation rate
	beta = beta_ann**(70/S) #Future consumption discount rate
	#sigma = 2 #Utility curvature parameter
	#sigma=rc.Sigma
	delta = 1-(1-delta_ann)**(70/S) #Depreciation Rate
	alpha = .3 #Capital Share of production
	e = np.ones((I, S, T+S+1)) #Labor productivities
	#A = np.ones(I) #Techonological Change, used for idential countries
	A=np.linspace(1,2,num=I) #Techonological Change, used for when countries are different

	diff=1e-6 #Convergence Tolerance
	distance=10 #Used in taking the norm, arbitrarily set to 10
	xi=.99 #Parameter used to take the convex conjugate of paths
	MaxIters=3000 #Maximum number of iterations on TPI.

	#Program Levers
	PrintAges = False #Prints the different key ages in the demographics
	PrintSS = False #Prints the result of the Steady State functions
	CalcTPI = True #Activates the calculation of Time Path Iteration
	#NOTE: Graphing only works if CalcTPI is activated.
	Graphs = False #Activates graphing the graphs.
	CountryNamesON = False #Turns on labels for the graphs. Replaces "Country x" with proper names.

	#MAIN CODE

	#Gets demographic data
	demog_params = (I, S, T, T_1, StartFertilityAge, EndFertilityAge, StartDyingAge, MaxImmigrantAge, g_A, PrintAges)
	FertilityRates, MortalityRates, Migrants, N_matrix, Nhat_matrix = Stepfuncs.getDemographics(demog_params)
	#Stepfuncs.plotDemographics((S,T),0,[0,19],"USA", N_matrix)

	#Initalizes initial guesses
	assets_guess = np.ones((I, S-1))*.15
	kf_guess = np.zeros((I))

	#Gets the steady state variables
	params_ss = (I, S, beta, sigma, delta, alpha, e, A)
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

	if CalcTPI==True: #Time Path Iteration, activated by line 24
		print "Beginning TPI"
		initialguess_params = (I, S, T, delta, alpha, e, A)
		assets_init, kf_init, w_initguess, r_initguess, kd_init, n_init, y_init, c_init = \
			Stepfuncs.get_initialguesses(initialguess_params, assets_ss, kf_ss, w_ss, r_ss)

		tp_params = (I, S, T, T_1, beta, sigma, delta, alpha, e, A, StartFertilityAge, StartDyingAge, N_matrix, MortalityRates, distance, diff, xi, MaxIters)
		wpath, rpath, cpath, Kpath, ypath = Stepfuncs.get_Timepath(tp_params, w_initguess, r_initguess, assets_init, kd_ss, kf_ss, w_ss, r_ss)
                #print "Cohorts:",S,"Countries:", I,"Curvature:", sigma
	
		if Graphs==True:
			Stepfuncs.plotTimepaths(I, S, T, wpath, rpath, cpath, Kpath, ypath, CountryNamesON)

if __name__=='__main__':
	sys.exit(main(sys.argv[1]))
