import numpy as np
import StepbyStepv1 as Stepfuncs
import time as time

np.set_printoptions(threshold = 3000, linewidth=2000, suppress=True)

#Parameters Zone
I_all = ["usa","eu","japan","china","india","russia","korea"]
I = 7 #Number of countries
S = 20 #Upper bound of age for agents
T = int(round(2.5*S)) #Number of time periods to convergence, based on Rick Evans' function.
I_touse = ["usa","eu","japan","china","russia","korea","india"]

T_1 = S #This is like TransYear in the FORTRAN I think
if S > 50:
	T_1 = 50

g_A = 0.015 #Technical growth rate

beta_ann=.95 #Starting future consumption discount rate
delta_ann=.08 #Starting depreciation rate
beta = beta_ann**(70/S) #Future consumption discount rate
sigma = 4 #Utility curvature parameter corresponds to 1/gamma
delta = 1-(1-delta_ann)**(70/S) #Depreciation Rate
alpha = .3 #Capital Share of production
chi = 1.5
rho = .4 #This can't be one or else you'll get a divide by zero. Elasticity of substitution of consumption of leisure and labor (less than 1, complementary)

diff = 1e-8 #Convergence Tolerance
demog_ss_tol = 1e-8 #Used in getting ss for population share
distance = 10 #Used in taking the norm, arbitrarily set to 10
xi = .9 #Parameter used to take the convex conjugate of paths
MaxIters = 500 #Maximum number of iterations on TPI.

#Program Levers
PrintAges = True #Prints the different key ages in the demographics
PrintLoc = False #Displays the current locations of the program inside key TPI functions
PrintSS = False #Prints the result of the Steady State functions
CalcTPI = True #Activates the calculation of Time Path Iteration
#NOTE: Graphing only works if CalcTPI is activated.
Graphs = True #Activates graphing the graphs.
CountryNamesON = False #Turns on labels for the graphs. Replaces "Country x" with proper names.
DiffDemog = True #Turns on different demographics over countries.
UseSSDemog = True
DiffProductivities = False
UseTape = False

LeaveHouseAge, FirstFertilityAge, LastFertilityAge, MaxImmigrantAge, FirstDyingAge, agestopull = Stepfuncs.getkeyages(S, PrintAges)

if len(I_touse) < I:
	print "WARNING: We are changing I from", I, "to", len(I_touse), "to fit the length of I_touse"
	I = len(I_touse)	
	time.sleep(2)

if DiffDemog:
	A = np.ones(I)+np.cumsum(np.ones(I)*.08)-.08 #Techonological Change, used for when countries are different
else:
	A = np.ones(I) #Techonological Change, used for idential countries

if DiffProductivities:
	e = np.ones((I, S, T+S))
	e[:,FirstDyingAge:,:] = 0.1
	e[:,:LeaveHouseAge,:] = 0.01
else:
	e = np.ones((I, S, T+S)) #Labor productivities

#MAIN CODE

#Gets demographic data
demog_params = (I, S, T, T_1, LeaveHouseAge, FirstFertilityAge, LastFertilityAge, FirstDyingAge, MaxImmigrantAge, agestopull, g_A, demog_ss_tol)
MortalityRates, Nhat_matrix, KIDs, Nhat_ss, KIDs_ss, lbar = Stepfuncs.getDemographics(demog_params, DiffDemog, I_all, I_touse)

Stepfuncs.plotDemographics((S,T),range(I),[-1], Nhat_matrix, I_touse)

#Initalizes initial guesses
assets_guess = np.ones((I, S-1))*.1
kf_guess = np.zeros((I))

#Gets the steady state variables
params_ss = (T, I, S, beta, sigma, delta, alpha, e[:,:,-1], A, FirstFertilityAge, FirstDyingAge, Nhat_ss, MortalityRates[:,:,-1], g_A, chi, rho, KIDs_ss)
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

if UseSSDemog == True:
	print "NOTE: USING SS DEMOGRAPHICS FOR TIMEPATH\n"
	Nhat_matrix = np.einsum("is,t->ist", Nhat_matrix[:,:,-1],np.ones(T+S))
	MortalityRates = np.einsum("is,t->ist", MortalityRates[:,:,-1],np.ones(T+S))
	KIDs = np.einsum("is,t->ist", KIDs[:,:,-1],np.ones(T+S))
	lbar[:] = 1
	time.sleep(2)


if CalcTPI==True: #Time Path Iteration, activated by line 24
	print "Beginning TPI"

	initialguess_params = (I, S, T, delta, alpha, e[:,:,0], A, FirstFertilityAge, FirstDyingAge, Nhat_matrix[:,:,0], MortalityRates[:,:,0], lbar[0], g_A, chi, rho, sigma, KIDs[:,:,0])
	assets_init, kf_init, wpath_initguess, rpath_initguess = \
		Stepfuncs.get_initialguesses(initialguess_params, assets_ss, kf_ss, w_ss, r_ss, PrintLoc)

	tp_params = (I, S, T, T_1, beta, sigma, delta, alpha, e, A, FirstFertilityAge, FirstDyingAge, Nhat_matrix, KIDs, MortalityRates, lbar, g_A, distance, diff, xi, MaxIters, rho, chi)

	wpath, rpath, cpath, Kpath, ypath = Stepfuncs.get_Timepath(tp_params, wpath_initguess, rpath_initguess, assets_init, kd_ss, kf_ss, w_ss, r_ss, PrintLoc, UseTape)
	
	if Graphs==True:
		Stepfuncs.plotTimepaths(I, S, T, wpath, rpath, cpath, Kpath, ypath, I_touse)

