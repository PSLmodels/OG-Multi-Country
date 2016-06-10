import numpy as np

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
        - self.g_A               = Scalar: Growth rate of technology
        - self.rho               = Scalar: The intratemporal elasticity of substitution between consumption and leisure
        - self.chi               = Scalar: Leisure preference parameter

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

def oldLoopMethod():
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
        if Print_caTimepaths:
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