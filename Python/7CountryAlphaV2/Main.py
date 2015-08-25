import numpy as np
import StepbyStepv1 as Stepfuncs

np.set_printoptions(linewidth=150)

#Parameters Zone
I = 3 #Number of countries
S = 80 #Upper bound of age for agents
T = int(round(2.5*S)) #Number of time periods to convergence, based on Rick Evans' function.

T_1 = S #This is like TransYear in the FORTRAN I think
if S > 50:
	T_1 = 50
FirstFertilityAge = int(np.round(S/80.*23))#The age when agents have their first children
LastFertilityAge = int(np.round(S/80.*45))#The age when agents have their last children
StartDyingAge = int(np.round(S/80.*68))#The first age agents can begin to die
MaxImmigrantAge = int(np.round(S/80.*65))#All immigrants are between ages 0 and MaxImmigrantAge
LeaveHouseAge = int(np.round(S/80.*21))

g_A = 0.015 #Technical growth rate

beta_ann=.95 #Starting future consumption discount rate
delta_ann=.08 #Starting depreciation rate
beta = beta_ann**(70/S) #Future consumption discount rate
sigma = 3 #Utility curvature parameter
delta = 1-(1-delta_ann)**(70/S) #Depreciation Rate
alpha = .3 #Capital Share of production
e = np.ones((I, S+1, T+S+1)) #Labor productivities
A = np.ones(I) #Techonological Change, used for idential countries
chi = 1.
rho = 1.5 #This can't be one or else you'll get a divide by zero.

diff = 1e-8 #Convergence Tolerance
demog_ss_tol = 1e-8 #Used in getting ss for population share
distance = 10 #Used in taking the norm, arbitrarily set to 10
xi = .7 #Parameter used to take the convex conjugate of paths
MaxIters = 500 #Maximum number of iterations on TPI.

#Program Levers
PrintAges = False #Prints the different key ages in the demographics
PrintLoc = False #Displays the current locations of the program inside key TPI functions
PrintSS = True #Prints the result of the Steady State functions
CalcTPI = False #Activates the calculation of Time Path Iteration
#NOTE: Graphing only works if CalcTPI is activated.
Graphs = False #Activates graphing the graphs.
CountryNamesON = False #Turns on labels for the graphs. Replaces "Country x" with proper names.
DiffDemog = True #Turns on different demographics over countries.

if DiffDemog == True:
	A = np.ones(I)+np.cumsum(np.ones(I)*.08)-.08 #Techonological Change, used for when countries are different

#MAIN CODE

#Gets demographic data
demog_params = (I, S, T, T_1, FirstFertilityAge, LastFertilityAge, StartDyingAge, MaxImmigrantAge, LeaveHouseAge, g_A, demog_ss_tol)
FertilityRates, MortalityRates, Migrants, N_matrix, Nhat_matrix, KIDs, Nhat_ss, KIDs_ss = Stepfuncs.getDemographics(demog_params, PrintAges, DiffDemog)

#Stepfuncs.plotDemographics((S,T),[0,1,2,3,4,5,6],[0],str("Initial"), Nhat_matrix)

#Initalizes initial guesses
assets_guess = np.ones((I, S-1))*.15
kf_guess = np.zeros((I))

#Gets the steady state variables
params_ss = (T, I, S, beta, sigma, delta, alpha, e, A, FirstFertilityAge, StartDyingAge, Nhat_ss, MortalityRates[:,:,-1], g_A, chi, rho,KIDs_ss)
assets_ss, kf_ss, kd_ss, n_ss, y_ss, r_ss, w_ss, c_vec_ss, ck_vec_ss = Stepfuncs.getSteadyState(params_ss, assets_guess, kf_guess)

if PrintSS==True: #Prints the results of the steady state, line 23 activates this
	print "assets steady state", assets_ss
	print "kf steady state", kf_ss
	print "kd steady state", kd_ss
	print "n steady state",n_ss
	print "y steady state", y_ss
	print "r steady state",r_ss
	print "w steady state", w_ss
	print "c_vec_ss steady state",c_vec_ss
        print "ck_vec_ss steady state", ck_vec_ss

if CalcTPI==True: #Time Path Iteration, activated by line 24
	print "Beginning TPI"
	initialguess_params = (I, S, T, delta, alpha, e, A, FirstFertilityAge, StartDyingAge, Nhat_matrix[:,:,0], MortalityRates[:,:,0], g_A)
	assets_init, kf_init, w_initguess, r_initguess, kd_init, n_init, y_init, c_init = \
		Stepfuncs.get_initialguesses(initialguess_params, assets_ss, kf_ss, w_ss, r_ss)

	tp_params = (I, S, T, T_1, beta, sigma, delta, alpha, e, A, FirstFertilityAge, StartDyingAge, Nhat_matrix, MortalityRates, g_A, distance, diff, xi, MaxIters)
	wpath, rpath, cpath, Kpath, ypath = Stepfuncs.get_Timepath(tp_params, w_initguess, r_initguess, assets_init, kd_ss, kf_ss, w_ss, r_ss, PrintLoc)
	
	if Graphs==True:
		Stepfuncs.plotTimepaths(I, S, T, wpath, rpath, cpath, Kpath, ypath, CountryNamesON)

