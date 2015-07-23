import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import WorldModule as WM


#TO Call any of the global variables/arrays, call on the world module, WM. I created this module to keep all of the golbal variables separate from the rest of the code. It was making the central code really messy to look at.


def ProbofSurvival2input(a, t):
    '''
    Gets the probabilities of an agent age a surviving to every age in year t

    Inputs:
        a -- Agent's current age (integer)
        t -- The current year (integer)

    Returns:
        probabilities -- vector of length 23 that contains the probability of reaching each age 68-90
        
    '''

    #The year the agent turns 68.
    birthday68th = 68 - a + t

    #Index of where to start the diagonal. >0 indicat:es it will start along the top of the matrix, and <0 indicates it will start along the left side of the matrix
    index = birthday68th - 2008

    #Rates of not dying for each year. May be less than size 23 if we don't have data for those years. Will be an empty list if index > 50
    survival_rates = np.diag((1-Russia.get_mortality_rate("Blah", True)), index)

    #Steady-state survival rates
    ss_survival_rates = (1-Russia.get_mortality_rate("Blah", True))[:,-1]

    #Stacks Probability of survival (1) for each year 0-68, and probability of survival for each age afterwards until we (maybe) ran out of data, and the steady state value for each year for which we don't have data 
    survival_rates = np.hstack((np.ones(68), survival_rates, ss_survival_rates[survival_rates.size:]))

    #Gets the probability of surviving to that year by taking the cumulative product
    probabilities = np.cumprod(survival_rates)

    return probabilities



def ProbofSurvival3input(a, i, t):
    '''
    Gets the probability of an agent age a at time t surviving to age i

    Inputs:
        a -- Agent's current age (integer between 0-90)
        i -- Age for which we are calculated the probability of surviving until (integer between 0-90)
        t -- The current year (integer greater than 2008)

    Returns:
        prob -- Probability of surviving to age i
        
    '''

    #%100 chance of surviving until 67
    if i < 68:
        return 1

    #The year when agent is age i
    yeartoget = t - a + i
    
    #The corresponding index of age and yeartoget in the mortailityrates matrix
    yearindex = yeartoget - 2008
    ageindex = i - 68
    
    #Returns steady-state value if the index is greater than the number of years for which we have data
    if yearindex > 50:
        yearindex = 50

    #Takes the product of all age from 68 to i for t
    prob = np.prod((1-Russia.get_mortality_rate("Blah", True))[0:ageindex+1, yearindex])
    print prob
    return prob

### Utility Function Section

def GetYearBecomingJ(Year,Gen,JJ): #Fortran 3704-3716

    GetYearJ=Year+JJ-Gen
    if (GetYearJ > WM.Years):
        GetYearJ=WM.Years
    if (GetYearJ < 0):
        GetYearJ=0
    if (Year==0):
        GetYearJ=0
    if (Year==WM.Years):
        GetYearJ=WM.Years

    return GetYearJ

def GetLifetimeUtility(Gen,Year,YClass,Country):
    '''
    Inputs:
        Gen-
        Year-
        YClass-
        Country-
    Returns:
        Lifetime utility for a particular agent
    '''
    lifetimeutil=0.0
    for j in xrange(Gen,WM.Gens):
        IK=GetYearBecomingJ(Year,Gen,j)
        lifetimeutil+=(1./(1.+WM.Delta[Gen,Year,Country]))**(j-Gen)*\
                WM.Survival_Probability[j,IK,Country]*((WM.Consump[j,IK,Yclass,Country]**(1.-(1./WM.Rho))+\
                WM.Alp*WM.Leisure[j,IK,YClass,Country]**(1.-(1./WM.Rho)))**((1.-(1./WM.Gamma))/(1.-(1./WM.Rho))))
        if j >=23 and j <=65:
            lifetimeutil+=(1./(1.+WM.Delta[Gen,Year,Country]))**(j-Gen)*WM.Survival_Probability[j,IK,Country]*\
                    (WM.Theta*Kids(j,IK,YClass,Country)*WM.Consump_Kids[j,IK,Yclass,Country]**(1.-(1./WM.Gamma)))

    lifetimeutil=lifetimeutil/(1.-(1./WM.Gamma))
    return lifetimeutil


####

def KIDS4input(a, i, t, k): #Corresponds to Fortran 3068-3102
    '''
    Gets the number of children of an agent age a at time t in skill class k when the agent is age i

    Inputs:
        a -- Agent's current age (integer between 0-90)
        i -- Some age > a
        t -- The current year (integer greater than 2008)
        k -- Agent's skill class

    Returns:
        kids -- Number of kids
        
    '''

    #No kids before age 23
    if i < 23:
        return 0

    #The year when agent is age i
    yeartoget = t - a + i
    
    #The corresponding indexes of age and yeartoget in the KIDS matrix
    yearindex = yeartoget - 2008
    

    kids = Russia.get_kids_mat(yeartoget)[i, k]
    print kids.shape
    print kids

    return kids

def KIDS3input(a, t):
    '''
    Gets the number of children of an agent age a at time t in skill class k when the agent is age i

    Inputs:
        a -- Agent's current age (integer between 0-90)
        t -- The current year (integer greater than 2008)

    Returns:
        Kids -- Array of size (GENS, SKILL_GROUPS) that contains numbers of kids for each age 0 to 90
        
    '''

    #The year the agent was born
    yearborn = t - a

    #The corresponding index to get the correct number of kids
    yearindex = yearborn - 2008

    #Gets the number of kids this agent will have over the course of his life. May cause an ERROR if it tries to go beyond the Max_Year
    Kids = np.transpose(np.diagonal(Russia.get_kids_mat("Blah", True), yearindex))

    return Kids

def Taxes(Year, Country): #Corresponds to Fortran 2200-2439
    '''
    Calculate Tax rates for Year and Country
    Inputs:
        Year - Inputs the year you wish to compute taxes for
        Country- Choose the country
    Returns:
        This changes the global variables        
    '''
    #First, we need to clear the variables
    WM.Total_Expenditures[Year,Country]=0.0
    TDAMP=.1
    RevPro=0.
    Tax_Base=np.zeros((3,WM.countries))
    Endogenous_Tax_Base=np.zeros((3,WM.countries))
    Tax_Base_Squared=np.zeros(WM.countries)
    Base=np.ones((WM.Yclasses+1,WM.countries))

    if WM.Iter >=20:
        TDAMP=.1
    IPL1=GetYearBecomingJ(Year,1,2)
    #WM.Education_Expenditures[Year,Country]=0.0

    #Calcualte Education Expenditures

    for i in xrange(WM.lasteducation+1):
        if Year>0:
            WM.Education_Expenditures[Year,Country]+=WM.Education_Expenditures_Ind[i,Country]*\
                    (WM.YY_0[Year,Country]/WM.YY[0,Country]*WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]/\
                    WM.POP[WM.startworkingage,0,WM.Yclasses,WM.firstcountry]*WM.POP_Efficient[0,Country]/\
                    WM.POP_Efficient[Year,Country])*(1.+WM.Tech)**(WM.startworkingage-i)*\
                    WM.POP[i,Year,WM.Yclasses,Country]/WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]
        elif Year<=0:
            WM.Education_Expenditures+=WM.Education_Expenditures_Ind[i,Country]*(1.+WM.Tech)**(WM.startworkingage-i)*\
                    WM.POP[i,Year,WM.Yclasses,Country]/WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]

    #2235 is where to pick up.
    #Calculating the total expenditures in year 301 and 0.
    if Year==WM.Years+1 or Year==0:
        WM.Government_Expenditures[Year,Country]=WM.Govs[Country]*WM.POP[91,Year,WM.Yclasses,Country]/\
                WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]+WM.Education_Expenditures[Year,Country]+\
                WM.Mu2_gov[Year,Country]*WM.Health_Benefits[Year,Country]
    #Calculate the discretionary spending
        WM.Government_Discretionary_Spending=WM.Govs[Country]*WM.POP[91,Year,WM.Yclasses,Country]/\
                WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]

    elif Year>0 and Year<WM.Years+1:

        #WM.Government_Expenditures[Year,Country]=0.0
        WM.Government_Discretionary_Spending[Year,Country]=0.0
        for j in xrange(WM.gens):
            WM.Government_Expenditures[Year,Country]+=WM.Govs[Country]*(WM.YY_0[Year,Country]/WM.YY[0,Country]*\
                    WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]/\
                    WM.POP[WM.startworkingage,0,WM.Yclasses,WM.firstcountry]*\
                    WM.POP_Efficient[0,Country]/WM.POP_Efficient[Year,Country])*\
                    WM.POP[j,Year,WM.Yclasses,Country]/WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]

        WM.Government_Discretionary_Spending[Year,Country]=WM.Government_Expenditures[Year,Country]

        WM.Government_Expenditures[Year,Country]+=WM.Education_Expenditures[Year,Country]+WM.Mu2_gov[Year,Country]*\
                WM.Health_Benefits[Year,Country]

    #Endogenous Tax rates are calculated to keep debt to GDP Constant
    WM.Debt[Year,Country]=WM.Debt_Level[Year,Country]*WM.YY[Year,Country]
    Revenues_Needed=WM.Government_Expenditures[Year,Country]+WM.RG[Year]*WM.Debt[Year,Country]+WM.Mu1*\
            WM.Pension_Benefits[Year,Country]+(WM.Mu2_tax[Year,Country]-WM.Mu2_gov[Year,Country])*\
            WM.Health_Benefits[Year,Country]+WM.Mu3[Year,Country]*WM.Disability_Benefits[Year,Country]-\
            ((1.+WM.NPOP[IPL1,WM.firstcountry])*(1.+WM.Tech)*WM.Debt[IPL1,Country]-WM.Debt[Year,Country])

    #Inheritance Tax Revenues
    WM.Tax_Revenues[3,Year,Country]=0.0
    for i in xrange(WM.Yclasses):
        for j in xrange(WM.WorkToDeath):
            WM.Tax_Revenues[3,Year,Country]+=WM.Inheritance_Tax_Rate[Year,Country]*Inheritances(j,Year,i,Country)*\
                    GetEfficientPopulation(j,Year,i,Country)

    WM.Tax_Revenues[5,Year,Country]=WM.Tax_Revenues[3,Year,Country]

    #Corporate Tax Revenues
    WM.Tax_Revenues[4,Year,Country]=WM.YY[Year,Country]-WM.Agg_Assets_Migrants[Year,WM.Yclasses,Country]
    for i in xrange(WM.Yclasses):
        WM.Tax_Revenues[4,Year,Country]-=WM.Wage_Index_Class[Year,i,Country]*WM.Labor[Year,i,Country]

    WM.Tax_Revenues[4,Year,Country]=(WM.Tax_Revenues[5,Year,Country]-WM.Del*WM.Capital[Year,Country])*\
            WM.Corp_Tax[Year,Country]

    #Calculate Corporate Tax Transfers
    if Country==WM.firstcountry:
        for i in xrange(WM.countries):
            for j in xrange (WM.Yclasses):
                for k in xrange(WM.WorkToDeath):
                    WM.H_Transfer[k,Year,j,i,Country]=WM.Assets[k,Year,j,i]/WM.Agg_Assets_World[Year]/\
                            (1.-WM.Mortality_rates[k,Year,i])*WM.Tax_Revenues[4,Year,Country]*WM.Mu4[Year,Country]
                    WM.Transfer[k,Year,j,i]+=WM.H_Transfer[k,Year,j,i,Country]

    WM.Tax_Revenues[4,Year,Country]=WM.Tax_Revenues[4,Year,Country]*(1.0-WM.Mu4[Year,Country])
    WM.Tax_Revenues[5,Year,Country]=WM.Tax_Revenues[5,Year,Country]+WM.Tax_Revenues[4,Year,Country]
    
    #Add Resource Endowment Revenues

    WM.Tax_Revenues[5,Year,Country]+=WM.Gov_Endow_Share[Year,Country]*WM.Endowment[Year]

    #Calculate Tax Bases
    Tax_Base[0,Country]=WM.CC[Year,Country]
    Tax_Base[1,Country]=0.

    for i in xrange(WM.Yclasses):
        Tax_Base[1,Country]+=WM.Wage_Index_Class[Year,i,Country]*WM.Labor[Year,i,Country]

    Tax_Base_Squared[Country]=Sum_Individual_Variables(Year,4,Country)

    Tax_Base[2,Country]=WM.Agg_Assets[Year,Country]*WM.RG[Year]

    #for i in xrange(WM.Yclasses+1):
        #WM.Agg_Marg_Wage_Tax_Rate[Year, i, Country]=0.
        #WM.Agg_Avg_Wage_Tax_Rate[Year,i,Country]=0.
        #Base[i,Country]=0.

    for i in xrange(3):
        Endogenous_Tax_Base[i,Country]=0.

    for i in xrange(2): #IPASS
        for j in xrange(3): #KK
            #Note that Revpro is the proportion of needed revenue covered by the endogenous tax
            WM.Aleph[j,Year,Country]=WM.Tax_Rate[j,0,Year,Country]
            WM.Beth[j,Year,Country]=WM.Tax_Rate[j,1,Year,Country]
            RevPro=WM.Endogenous_Tax_Ratio[j,Year,Country]
            if (i==0 and WM.Aleph[j,Year,Country]==-1) or (WM.Beth[j,Year,Country]==-1):
                Endogenous_Tax_Base[j,Country]=Tax_Base[j,Country]

            if ((i==0 and WM.Aleph[j,Year,Country] != 1) and (WM.Beth[j,Year,Country] != -1.)) or\
                    ((i==1 and WM.Aleph[j,Year,Country]==-1) or WM.Beth[j,Year,Country]==-1.):
                        if i==1 and WM.Aleph[j,Year,Country]==-1.:
                            ABBACAB=1
                            #print "Fix this"
                            #WM.Aleph[j,Year,Country]=(Revx*RevPro-WM.Beth[j,Year,Country]*\
                                    #Tax_Base_Squared[Country]/2.)/Endogenous_Tax_Base[j,Country]
                        if i==1 and WM.Beth[j,Year,Country]==-1.:
                            WM.Beth[j,Year,Country]=2.*(Revx*RevPro-WM.Aleph[j,Year,Country]*Endogenous_Tax_Base[j,Country])/\
                                    Tax_Base_Squared[Country]

                        if j==0:
                            WM.Consump_Price[Year,Country]=WM.Consump_Price[Year,Country]*(1.-TDAMP)+\
                                    (1.+WM.Aleph[j,Year,Country])*TDAMP
                            WM.Tax_Revenues[j,Year,Country]=(WM.Consump_Price[Year,Country]-1.)*WM.CC[Year,Country]
                            WM.Tax_Revenues[5,Year,Country]+=WM.Tax_Revenues[j,Year,Country]

                        elif j==1:
                            WM.Tax_Revenues[j,Year,Country]=0.
                            for k in xrange(WM.Yclasses):
                                for l in xrange(WM.WorkToDeath):
                                    Z=(WM.Hours-WM.Leisure[l,Year,k,Country])*GetWage(l,Year,k,Country)

                                    WM.Marg_Wage_Tax_Rate[l,Year,k,Country]=TDAMP*(1.-WM.Aleph[j,Year,Country]-\
                                            WM.Beth[j,Year,Country]*Z)+(1.0*TDAMP)*WM.Marg_Wage_Tax_Rate[l,Year,k,Country]

                                    WM.Avg_Wage_Tax_Rate[l,Year,k,Country]=TDAMP*(1.-WM.Aleph[j,Year,Country]\
                                            -WM.Beth[j,Year,Country]*Z/2.)+(1.-TDAMP)*WM.Avg_Wage_Tax_Rate[l,Year,k,Country]

                                    WM.Agg_Marg_Wage_Tax_Rate[Year,k,Country]+=(WM.Aleph[j,Year,Country]+WM.Beth[j,Year,Country]*\
                                            Z)*Z*GetEfficientPopulation(l,Year,k,Country)
                                    
                                    WM.Agg_Marg_Wage_Tax_Rate[Year,WM.Yclasses,Country]+=\
                                            WM.Agg_Marg_Wage_Tax_Rate[Year,k,Country]

                                    WM.Agg_Avg_Wage_Tax_Rate[Year,k,Country]+=(WM.Aleph[j,Year,Country]+WM.Beth[j,Year,Country]\
                                            *Z/2.)*Z*GetEfficientPopulation(l,Year,k,Country)
                                    
                                    WM.Agg_Avg_Wage_Tax_Rate[Year,WM.Yclasses,Country]+=WM.Agg_Avg_Wage_Tax_Rate[Year,k,Country]

                                    Base[k,Country]+=Z*GetEfficientPopulation(l,Year,k,Country)

                                    Base[WM.Yclasses,Country]+=Base[k,Country]

                                    WM.Tax_Revenues[j,Year,Country]+=(1.-WM.Avg_Wage_Tax_Rate[l,Year,k,Country])*Z*\
                                            GetEfficientPopulation(l,Year,k,Country)

                            WM.Tax_Revenues[5,Year,Country]+=WM.Tax_Revenues[j,Year,Country]
                        elif j==2:
                            WM.Capital_Tax_Rate[Year, Country]=WM.Capital_Tax_Rate[Year,Country]*(1.-TDAMP)+(1-\
                                    WM.Aleph[j,Year,Country])*TDAMP
                            WM.Tax_Revenues[j,Year,Country]=(1.-WM.Capital_Tax_Rate[Year,Country])*WM.Agg_Assets[Year,Country]*\
                                    WM.RG[Year]
                            WM.Tax_Revenues[5,Year,Country]+=WM.Tax_Revenues[j,Year,Country]


        if i==0:
            Revx=Revenues_Needed-WM.Tax_Revenues[5,Year,Country]

    if Base[0,Country]<=0:
        WM.Agg_Marg_Wage_Tax_Rate[Year,:,Country]=0
        WM.Agg_Avg_Wage_Tax_Rate[Year,:,Country]=0
    if Base[1,Country]<=0:
        WM.Agg_Marg_Wage_Tax_Rate[Year,:,Country]=0
        WM.Agg_Avg_Wage_Tax_Rate[Year,:,Country]=0
    if Base[2,Country]<=0:
        WM.Agg_Marg_Wage_Tax_Rate[Year,:,Country]=0
        WM.Agg_Avg_Wage_Tax_Rate[Year,:,Country]=0

    else:
        for i in xrange(WM.Yclasses+1):
            WM.Agg_Marg_Wage_Tax_Rate[Year,i,Country]=WM.Agg_Marg_Wage_Tax_Rate[Year,i,Country]/Base[i,Country]
            WM.Agg_Avg_Wage_Tax_Rate[Year,i,Country]=WM.Agg_Avg_Wage_Tax_Rate[Year,i,Country]/Base[i,Country]


    return None

def annualMargProd(Year, Country): #Corresponds to Fortran 2167-2193
    '''
    Compute wages for all countries, capital productivity only for the US from production function
    Inputs:
        Year- Current year
        Country- Country number
    Returns:
        NONE        
    '''
    for i in xrange(WM.Yclasses):
        HLabor=1.
        for j in xrange(1,WM.Yclasses):
            if j != i:
                HLabor=HLabor*WM.Labor[Year,j,Country]**WM.Beta[j]
        WM.Wage_Index_Class[Year,i,Country]=WM.Beta[i]*WM.AP[Year]*WM.Captial[Year,Country]**WM.Alpha*WM.Labor[Year,i,Country]**(WM.Beta[i]-1.)*HLabor

    if Country==WM.firstcountry:
        HLabor=1.
        for i in xrange(1, WM.Yclasses):
            HLabor=HLabor*WM.Labor[Year,i,Country]**WM.Beta[i]

        WM.R[Year,Country]=WM.Alpha*WM.AP[Year]*WM.Capital[Year,Country]**(WM.Alpha-1.)*HLabor

    return None

def Get_Fraction_PVE_Consumed(Gen,Year,YClass,Country): #Russia code: 2266-2315
    '''
    It computes the fraction of PVE consumed at age startworking age
    Input:
        Gen-
        Year-
        YClass-
        Country-
    Output:
        Get Fraction PVE Consumed.
    '''
    FractionPVEConsumed=0.0
    S=1.
    for j in xrange (Gen, WM.gens+1):
        IK=GetYearBecomingJ(Year,Gen,j)
        #print "IK", IK
        IKP1=GetYearBecomingJ(Year,Gen,j+1)
        if j==Gen:
            PC1=WM.Consump_Price[IK,Country]
        H=(1.+WM.Alp**WM.Rho*(Sum_Wage(j,IK,YClass,Country)*(WM.Marg_Wage_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Pension_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Health_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Disability_Tax_Rate[j-WM.startworkingage,IK,YClass,Country])\
                /WM.Consump_Price[IK,Country])**(1.-WM.Rho))**\
                ((WM.Rho-WM.Gamma)/(1.-WM.Rho))
        #print "H", H

        if j==Gen:
            H1=H

        #The Amount of consumption in a given year is importantly a function of DELTA, the time preference parameter.
        WM.ZXI[j-WM.startworkingage,YClass,Country]=(S/(1+WM.Delta[Gen,Year,Country])**(j-Gen)*PC1/WM.Consump_Price[IK,Country]*\
                WM.Survival_Probability[j,IK,Country]/WM.Survival_Probability[Gen,Year,Country])**WM.Gamma*H/H1

        
        FractionPVEConsumed+=WM.ZXI[j-WM.startworkingage,YClass,Country]*\
                WM.Consump_Price[IK,Country]*(1.+WM.Theta*KIDS(j,IK,YClass,Country)/\
                H+Sum_Wage(j,IK,YClass,Country)*(WM.Avg_Wage_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Avg_Pension_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Avg_Health_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Avg_Disability_Tax_Rate[j-WM.startworkingage,IK,YClass,Country])/WM.Consump_Price[IK,Country]*\
                (Sum_Wage(j,IK,YClass,Country)*(WM.Marg_Wage_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Pension_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Health_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Disability_Tax_Rate[j-WM.startworkingage,IK,YClass,Country])/\
                (WM.Alp*WM.Consump_Price[IK,Country]))**(-WM.Rho))/S
        
        #print "Fraction",FractionPVEConsumed

    
        S*=(1.+WM.RG[IKP1]*WM.Capital_Tax_Rate[IKP1,Country])

        #print "S", S



    #FractionPVEConsumed=1./FractionPVEConsumed

    return FractionPVEConsumed
    #return None

def IndividualDemands(Gen, Year, YClass,Country): #Corresponds to Fortran 1837-1898
    '''
    Computes the lifetime consumption,leiure and asset path for person age GEn in Year with parents age GG
    at the person's birth for the rest of his life
    Inputs:
        Gen-
        Year-
        Y Class-
        Country-
    Returns:
        Alters the global arrays        
    '''
    H=0.
    #Sets up initial assets at beginning of year IYR.
    if Gen >WM.startworkingage:
        WM.Assets[Gen,Year,YClass,Country]=WM.Base_Assets[Gen,YClass,Country]
    if Gen==WM.startworkingage and Year>=0:
        WM.Assets[Gen,Year,YClass,Country]-0.0

    #C1 is consumption in Year for this cohort. The nature of the utility function is such that
    #All future consumption can be derived fro C1 by a factor, ZXIN(Gen) (Gen is age), a function of tax rates,
    #interest rates and wages.
    C1=WM.Consump[Gen,Year,YClass,Country]

    for j in xrange(Gen, WM.gens):
        IK=GetYearBecomingJ(Year,Gen,j)

        WM.Consump[j-WM.startworkingage,IK,YClass,Country]=WM.ZXI[j-WM.startworkingage,YClass,Country]*C1

        WM.Leisure[j-WM.startworkingage,IK,YClass,Country]=WM.Consump[j-WM.startworkingage,IK,YClass,Country]*WM.Alp**WM.Rho*\
                (Sum_Wage(j,IK,YClass,Country)*(WM.Marg_Wage_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Pension_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]- \
                WM.Marg_Health_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Disability_Tax_Rate[j-WM.startworkingage,IK,YClass,Country])/WM.Consump_Price[IK,Country])**(-WM.Rho)

        H=(1.+WM.Alp**WM.Rho*(Sum_Wage(j,IK,YClass,Country)*(WM.Marg_Wage_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Pension_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Health_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Marg_Disability_Tax_Rate[j-WM.startworkingage,IK,YClass,Country])/\
                WM.Consump_Price[IK,Country])**(1.-WM.Rho))**\
                ((WM.Rho-WM.Gamma)/(1.-WM.Rho))
        if j>= WM.startfertilityage and j<=65:
            WM.Consump_Kids[j-WM.startworkingage,IK,YClass,Country]=WM.Theta*WM.Consump[j-WM.startworkingage,IK,YClass,Country]/H
        if j>Gen:
            IKM1=GetYearBecomingJ(Year,Gen,j-1)-WM.startworkingage
            IM1=j-WM.startworkingage-1
            #New Assets= Old Assets + Interest + Net Earnings - Consumption
            WM.Assets[j-WM.startworkingage,IK,YClass,Country]=(1.+WM.RG[IKM1]*WM.Capital_Tax_Rate[IKM1,Country])*\
                    (WM.Assets[IM1,IKM1,YClass,Country]+Inheritances(IM1,IKM1,YClass,Country))-\
                    WM.Inheritance_Tax_Rate[IKM1,Country]*Inheritances(IM1,IKM1,YClass,Country)+\
                    (WM.Hours-WM.Leisure[IM1,IKM1,YClass,Country])*Sum_Wage(IM1,IKM1,YClass,Country)*\
                    (WM.Avg_Wage_Tax_Rate[IM1,IKM1,YClass,Country]-WM.Avg_Pension_Tax_Rate[IM1,IKM1,YClass,Country]-\
                    WM.Avg_Health_Tax_Rate[IM1,IKM1,YClass,Country]-WM.Avg_Disability_Tax_Rate[IM1,IKM1,YClass,Country])+\
                    WM.Pension_Benefits_Ind[IM1,IKM1,YClass,Country]-WM.Consump[IM1,IKM1,YClass,Country]*\
                    KIDS(IM1,IKM1,YClass,Country)*WM.Consump_Kids[IM1,IKM1,YClass,Country]

            WM.Assets[j-WM.startworkingage,IK,YClass,Country]+=(WM.Health_Benefits_Ind[IM1,IKM1,Country]+\
                    Kids_Health_Benefits(IM1,IKM1,YClass,Country))*(1.-WM.Mu2_gov[IKM1,Country])+\
                    WM.Transfer[IM1,IKM1,YClass,Country]

            if j>= 22 and j<=65:
                WM.Assets[j-WM.startworkingage,IK,YClass,Country]+=WM.Disability_Benefits_Ind[IKM1,Country]

    return None

def TransitionPath(): #Fortran 1432-1728
    '''
    Use the power method to compute the dominant eigenvalue.

    Inputs:
        A -- numpy array, assume it has a dominant eigenvalue
        tol -- float tolerance
    Returns:
        eigenvalue -- dominant eigenvalue of A
        eigenvector -- corresponding eigenvector of A
        
    '''
    return None

def Get_PV_Earnings(Gen,Year,YClass,Country): #222-12264 RUSSIA CODE
    '''
    Gets the present value of lifetime earnings
    Inputs:
        Gen-
        Year-
        YClass-
        Country-
    Output:
        PresentValue earnings-
    '''
    PVEarnings=0.
    #Set initial assets
    if Gen>WM.startworkingage:
        PVEarnings=WM.Base_Assets[Gen-WM.startworkingage,YClass,Country]*(1.+WM.RG[Year]*WM.Capital_Tax_Rate[Year,Country])

    #Discount Rate
    S=1.
    #Compute the present value for year corresponding to age j
    #This is calculated by adding in everything the person will earn and transferred over their lifetime
    for j in xrange(Gen, WM.gens+1):
        IK=GetYearBecomingJ(Year,Gen,j)
        IKP1=GetYearBecomingJ(Year,Gen,j+1)

        PVEarnings+=(WM.Hours*Sum_Wage(j,IK,YClass,Country)*(WM.Avg_Wage_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Avg_Pension_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Avg_Health_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]-\
                WM.Avg_Disability_Tax_Rate[j-WM.startworkingage,IK,YClass,Country]))/S

        PVEarnings+=((1.-WM.Inheritance_Tax_Rate[IK,Country]+WM.RG[IK]*WM.Capital_Tax_Rate[IK,Country])*\
                Inheritances(j,IK,YClass,Country))/S

        PVEarnings+=(WM.Health_Benefits_Ind[j-WM.startworkingage,IK,Country]+\
                Kids_Health_Benefits(j-WM.startworkingage,IK,YClass,Country))*\
                (1.-WM.Mu2_gov[IK,Country])/S

        PVEarnings+=WM.Transfer[j-WM.startworkingage,IK,YClass,Country]

        #Pay attention here if you're changing the retirement age for developed countries
        if j>=WM.startworkingage and j<=91:
            PVEarnings+=WM.Disability_Benefits_Ind[IK,Country]/S
        if (j>= WM.RetirementAge[IK,Country]):
            PVEarnings+=WM.Pension_Benefits_Ind[j-WM.startworkingage,IK,YClass,Country]/S

        S*=(1+WM.RG[IKP1]*WM.Capital[IKP1,Country])

    return PVEarnings

def SteadyState(): #Fortran 1310-1425, there's a lot of other functions that need to happen first, I'll come back to this later.
    '''
    Use the power method to compute the dominant eigenvalue.

    Inputs:
        None. It takes the global variables and calculates the steady state. 
    Returns:
        Steady State        
    '''
    Hlabor=0.
    j=0
    #Sets the initial guesses for interest rate, wages and producer prices
    WM.RG[-1]=.2
    WM.Wage_Index_Class[-1,0,:]=1.
    WM.Wage_Index_Class[-1,1,:]=2.

    #Now that the initial guesses are made, we then iterate. If the model doesn't converge by the maximum iterations
    #Then the model moves on to the next step using the last guess.
    for j in xrange(WM.MaxIter+1):
        #Calculate AP such that the wage rate of the first income class in the US is one
        #AP is essentially a total factor productivity term
        WM.AP[-1]=1.0/(WM.Beta[0]*WM.Capital[-1,WM.firstcountry]**WM.Alpha*\
                WM.Labor[-1,0,WM.firstcountry]**(WM.Beta[0]-1.)*WM.Labor[-1,1,WM.firstcountry]**WM.Beta[1])

        #Global Interest Rate
        WM.RG[-1]=(1.-WM.Corp_Tax[-1,WM.firstcountry])*(WM.R[-1,WM.firstcountry]-WM.Del)

        #Calcualte interest rates in all non US countries by taking the global interest rate
        #and adding in taxes and depreciation
        if WM.lastcountry != WM.firstcountry:
            for k in xrange(WM.firstcountry+1,WM.lastcountry+1):
                WM.R[-1,k]=WM.RG[-1]/(1.-WM.Corp_Tax[-1,k]+WM.Del)
                Hlabor=1.
                for p in xrange(WM.Yclasses):
                    Hlabor*=WM.Labor[-1,p,k]**WM.Beta[p]


        #Then we calculate capital demand in these countries using the normal cobb-douglas FOC.
                WM.Capital[-1,k]=(WM.R[-1,k]/(WM.Alpha*WM.AP[-1]*Hlabor))**(1./(WM.Alpha-1.))

        #HOUSEHOLDS
        #Calculate individual demands for all people in all countries, the exact calculations are done
        #in separate subroutines
        for k in xrange(WM.countries):
            for p in xrange(WM.Yclasses):
                WM.Consump[WM.startworkingage,-1,p,k]=Get_PV_Earnings(WM.startworkingage,-1,p,k)*\
                        Get_Fraction_PVE_Consumed(WM.startworkingage,-1,p,k)
                IndividualDemands(WM.startworkingage,-1,p,k)
                ShadowWage(WM.startworkingage,-1,p,k)
        #Having Calculated individual demands, call the function which aggregates them
        Get_Aggregate_Variables(-1)

        #Production Side
        for k in xrange(WM.countries):
            Get_Marginal_Products(-1,k)

        #Government Sector
        for k in xrange(WM.countries):
            Taxes(-1,k)

        #Get the Health System, remember that developing countries have no disability insurance
        for k in xrange(WM.countries):
            Health_System(-1,k)
            Disability_Insurance(-1,k)
            PensionSystem(-1,k)
        #Prints, so we can see the incremental change
        #print j
        WM.Iter=j
        if (GoodsMarket(-1)==1):
            break
        if (j==WM.MaxIter): #WM.MaxIter):
            print "Caution! Convergence not reached!!"

    for q in xrange(WM.countries):
        if q==0:
            Country_Name='MURICA'
        if q==1:
            Country_Name='Europe'
        if q==2:
            Country_Name='Rising Sun'
        if q==3:
            Country_Name='Middle Kingdom'
        if q==4:
            Country_Name='India'
        if q==5:
            Country_Name='Mother Russia'
        if q==6:
            Country_Name='South Korea'
        print Country_Name, "\n"
        #print output(-1,1,WM.YClasses,q)




    return None

def Disability(year, country): #Fortran 2645-2772
    '''
    Use the power method to compute the dominant eigenvalue.

    Inputs:
        A -- numpy array, assume it has a dominant eigenvalue
        tol -- float tolerance
    Returns:
        eigenvalue -- dominant eigenvalue of A
        eigenvector -- corresponding eigenvector of A
        
    '''
    TDAMP = 1.
    TDAM = 1.
    t = year
    r = country

    WM.Disability_Benefits[:,:] = 0
    Earnings_Base_Payroll_Ind = np.zeros((WM.gens+1-WM.startworkingage, WM.Yclasses))
    Earnings_Base_Payroll = np.zeros(WM.Yclasses+1)

    #Google translate says: Only turn on when HEALTH_SYSTEM is missing
    #CALL GET_AVG_LABOR_EARNINGS(YEAR, COUNTRY)

    Max_Taxable_Earnings = WM.Contribution_Ceiling[t, r]*WM.Avg_Labor_Earnings[t, r]

    if t > 0:
        WM.Disability_Benefits_Ind[t, r] = WM.Disability_Benefits_Ind[t, r]   \
                                                    *WM.YY_0[t, r]   \
                                                    /WM.YY[0, r]   \
                                                    *WM.POP[WM.startworkingage, t, -1, WM.firstcountry]   \
                                                    *WM.POP_Efficient[0, r]   \
                                                    /WM.POP_Efficient[t, r]
    for g in range(WM.startworkingage, 64+1):
        WM.Disability_Benefits[t, r] = WM.Disability_Benefits[t, r] + WM.Disability_Benefits_Ind[t, r]   \
                                                                                    *(1+WM.Tech)**(WM.startworkingage-g)   \
                                                                                    *WM.POP[WM.startworkingage, t, -1, WM.firstcountry]
    for k in range(WM.Yclasses):
        Earnings_Base_Payroll[k] = 0
        for g in range(WM.startworkingage, WM.RetirementAge[t, r]):
            Earnings_Base_Payroll_Ind[g, k] = dummy_Sum_Wage(g, t, k, r) * (WM.Hours - WM.Leisure[g, t, k, r])
            if Earnings_Base_Payroll_Ind[g, k] < 0: Earnings_Base_Payroll_Ind[g, k] < 0
            if Earnings_Base_Payroll_Ind[g, k] <= Max_Taxable_Earnings:
                Earnings_Base_Payroll[k] = Earnings_Base_Payroll[k] + Earnings_Base_Payroll_Ind[g, k]*GetEfficientPopulation(g, t, k, r)
            else:
                Earnings_Base_Payroll[k] = Earnings_Base_Payroll[k]+Max_Taxable_Earnings*GetEfficientPopulation(g, t, k, r)
        Earnings_Base_Payroll[-1] = Earnings_Base_Payroll[-1] + Earnings_Base_Payroll[k]

    if Earnings_Base_Payroll[-1] !=0:
        WM.Agg_Disability_Tax_Rate[t, r] = (1-WM.Mu3[t, r]) * Disability_Benefits[t, r] / Earnings_Base_Payroll[-1]
    else:
        print "\nDisability_Benefits Function wants me to print something here when Earnings_Base_Payroll[-1] == 0. I don't know what, but look at line 2700 of the Fortran\n"
        time.sleep(4)

    WM.Agg_Marg_Disability_Tax_Rate[t, -1, r] = 0
    WM.Agg_Avg_Disability_Tax_Rate[t, -1, r] = 0
    for k in range(WM.Yclasses):
        WM.Agg_Marg_Disability_Tax_Rate[t, k, r] = 0
        WM.Agg_Avg_Disability_Tax_Rate[t, k, r] = 0
        for g in range(WM.startworkingage, WM.RetirementAge[t,r]):
            if Earnings_Base_Payroll[k] <= Max_Taxable_Earnings:
                if k < WM.Yclasses:
                    WM.Marg_Disability_Tax_Rate[g, t, k, r] * (1-TDAMP) + TDAMP * WM.Agg_Disability_Tax_Rate[t, r]
                else:
                    WM.Marg_Disability_Tax_Rate[g, t, k, r] = 0
                WM.Avg_Disability_Tax_Rate[g, t, k, r] = WM.Avg_Disability_Tax_Rate[g, t, k, r] * \
                        (1-TDAMP) + TDAMP * WM.Agg_Disability_Tax_Rate[t, r]
                WM.Agg_Marg_Disability_Tax_Rate[t, k, r] += WM.Marg_Disability_Tax_Rate[g, t, k, r]*\
                        Earnings_Base_Payroll_Ind[g, k] * GetEfficientPopulation(g, t, k, r)
                WM.Agg_Avg_Disability_Tax_Rate[t, k, r] += WM.Avg_Disability_Tax_Rate[g, t, k, r] *\
                        Earnings_Base_Payroll_Ind[g, k] * GetEfficientPopulation(g, t, k, r) 
            else:
                WM.Marg_Disability_Tax_Rate[g, t, k, r] = 0
                WM.Avg_Disability_Tax_Rate[g, t, k, r] = Avg_Disability_Tax_Rate[g, t, k, r] * \
                        (1-TDAM) + TDAM * Agg_Disability_Tax_Rate[t, r] * Max_Taxable_Earnings /\
                        Earnings_Base_Payroll_Ind[g, k]
                WM.Agg_Marg_Disability_Tax_Rate[t, k, r] += WM.Marg_Disability_Tax_Rate[g, t, k, r] *\
                        Max_Taxable_Earnings * GetEfficientPopulation(g, t, k, r)
                WM.Agg_Avg_Disability_Tax_Rate[t, k, r] += WM.Avg_Disability_Tax_Rate[g, t, k, r] *\
                        Max_Taxable_Earnings * GetEfficientPopulation(g, y, k, r)
            if WM.Agg_Disability_Tax_Rate[t, r] == 0:
                WM.Marg_Disability_Tax_Rate[g, t, k, r] = 0
                WM.Avg_Disability_Tax_Rate[g, t, k, r] = 0
                WM.Agg_Marg_Disability_Tax_Rate[t, k, r] = 0
                WM.Agg_Avg_Disability_Tax_Rate[t, k, r] = 0

        WM.Agg_Marg_Disability_Tax_Rate[t, -1, r] += WM.Agg_Marg_Disability_Tax_Rate[t, k, r]
        WM.Agg_Avg_Disability_Tax_Rate[t, -1, r] += WM.Agg_Avg_Disability_Tax_Rate[t, k, r]
        WM.Agg_Marg_Disability_Tax_Rate[t, k, r] = WM.Agg_Marg_Disability_Tax_Rate[t, k, r]/Earnings_Base_Payroll[k]
        WM.Avg_Disability_Tax_Rate[t, k, r] = WM.Agg_Avg_Disability_Tax_Rate[t, k, r]/Earnings_Base_Payroll[k]

    WM.Agg_Marg_Disability_Tax_Rate[t, -1, r] = WM.Agg_Marg_Disability_Tax_Rate[t, -1, r]/Earnings_Base_Payroll[-1]
    WM.Agg_Avg_Disability_Tax_Rate[t, -1, r] = WM.Agg_Avg_Disability_Tax_Rate[t, -1, r]/Earnings_Base_Payroll[-1]

    return None

def Hicksian(): #4402-4457
    '''
    Use the power method to compute the dominant eigenvalue.

    Inputs:
        A -- numpy array, assume it has a dominant eigenvalue
        tol -- float tolerance
    Returns:
        eigenvalue -- dominant eigenvalue of A
        eigenvector -- corresponding eigenvector of A
        
    '''
    return None

def PrintToFile(): #4051-4397 (Want to merge with 3927-4044)
    '''
    Use the power method to compute the dominant eigenvalue.

    Inputs:
        A -- numpy array, assume it has a dominant eigenvalue
        tol -- float tolerance
    Returns:
        eigenvalue -- dominant eigenvalue of A
        eigenvector -- corresponding eigenvector of A
        
    '''
    return None



def PopulationShare(IYR,IAU,IAO,Country): #3680-3695 #Not sure what it does. The Fortran code doesn't really explain much
    '''
    Use the power method to compute the dominant eigenvalue.

    Inputs:
        IYR-
        IAU-
        IAO-
        Country
    Returns:
        PopulationShare 
    '''
    IAYR=IYR-2000 #Model is normalized to 2000, 2000=0
    PopulationShare=0.0
    for i in range (IAU, IAO):
        PopulationShare+= WM.POP[I,IAYR,Yclasses+1,Country]

    PopulationShare = PopulationShare/ WM.POP[91,IAYR,Yclasses+1,Country]*100
    
    
    return PopulationShare


def ReadData(): #3208-3525 (This will likely be where James' code will be put in)
    '''
    Use the power method to compute the dominant eigenvalue.

    Inputs:
        IYR-
        IAU-
        IAO-
        Country
    Returns:
        PopulationShare        
    '''
    return None

def KidsHealth(Gen, Year, YClass, Country): #3164-3201
    '''
    Computes health benefits for kids of parents age GEN in YEAR
    Inputs:
        Gen-
        Year-
        YClass-
        Country-
    Returns:
        Kids Health Benefits-        
    '''
    HealthBenefits=0.
    UU=1.
    if Gen <= 46:
        UU=1
    else:
        UU=Gen-45

    if Gen- WM.startfertilityage < 20:
        MM=Gen-WM.startfertilityage
    else:
        MM=20

    if Gen <= 22 or Gen >= 66:
        HealthBenefits=0.
    if Gen >= WM.startfertilityage and Gen <=65:
        if Gen <=45 and Year >0:
            HealthBenefits=WM.Fertility[Gen,Country]*WM.Health_Benefits_Ind[0,Year,Country]
        if Gen <=45 and Year ==0:
            HealthBenefits=WM.Fertility_ss[Gen,Country]*WM.Health_Benefits_Ind[0,Year,Country]
        for j in xrange(UU,MM):
            if Year > 0:
                HealthBenefits=HealthBenefits+WM.POP[j,Year,YClass,Country]\
                *WM.Fertility[Gen-j,Year-j,Country]/WM.Total_Fertility[Year-j,Country]\
                /WM.POP[Gen,Year,YClass,Country]*WM.Health_Benefits_Ind[j,Year,Country]
            if Year ==0:
                HealthBenefits=HealthBenefits+WM.POP[j,Year,YClass,Country]*\
                WM.Fertility_ss[Gen-j,Country]/WM.Total_Fertility_ss[Country]\
                /WM.POP[Gen,Year,YClass,Country]*WM.Health_Benefits_Ind[j,Year,Country]


    return HealthBenefits


def Inheritances(Gen, Year, YClass, Country): #3140-3157
    '''
    Use the power method to compute the dominant eigenvalue.

    Inputs:
        Gen-
        Year-
        YClass-
        Country
    Returns:
        eigenvalue -- dominant eigenvalue of A        
    '''

    if Gen >= 50:
        Inheritance=0.
    else : #Age can only be between 21-49
        Inheritance= WM.shr[Gen]/100 * WM.Agg_Assets_For_Bequests[Year,YClass,Country]/\
                WM.POP[Gen,Year,YClass,Country]*(1.+WM.Tech)**(Gen-WM.startworkingage)

    return Inheritance

def Sum_Wage(J, IK, YClass, Country): #3109-3118
    '''
    
    Inputs:
        J-
        IK-
        YClass-
        Country-
    Returns:
        Adds the Wage and Shadow Wage together        
    '''
    K=J-WM.startworkingage
    #This is for the purposes of converting the age into a python compatable index

    SWage=GetWage(J,IK,YClass,Country)+WM.Shadow_Wage[K,IK,YClass,Country]
   
    return SWage

def PensionSystem(Year, Country): #2780-2914
    '''
    Calculates the compoents of the countriy's pension system
    Inputs:
        Country -- Indicates which country you're calculating
        Year-- Indicates the year
    Returns:
       Doesn't return a value, rather, it calculates the changes.
    '''

    TDAMP=.2
    TDAM=.2

    #First you need to compute average indexed labor earnings
    AvgIndexedEarnings(Year,Country)

    EarningsBasePayroll=np.zeros(WM.Yclasses+1)
    EarningsBasePayrollIND=np.zeros((WM.WorkToDeath,WM.Yclasses))

    #The Do loops that the fortran code uses are considered redunant since they're already all zeros

    #WM.Pension_Benefits[Year,Country]=0.
    #EarningsBasePayroll[WM.Yclasses]=0.
    MaxTaxableEarnings=WM.Contribution_Ceiling[Year,Country]*WM.Avg_Labor_Earnings[Year,Country]
    #print "TEST", WM.Avg_Labor_Earnings[Year,Country]

    for i in xrange(WM.Yclasses):
        #EarningsBasePayroll=0.
        for j in xrange(WM.startworkingage,WM.RetirementAge[Year,Country]-1):
            EarningsBasePayrollIND[j,i]=Sum_Wage(j,Year,i,Country)*(WM.Hours-WM.Leisure[j,Year,i,Country])

            if EarningsBasePayrollIND[j,i]<0.:
                EarningsBasePayroll[j,i]=0.

            if EarningsBasePayrollIND[j,i] <= MaxTaxableEarnings:
                EarningsBasePayroll[i]+=EarningsBasePayrollIND[j,i]*GetEfficientPopulation(j,Year,i,Country)

            else:
                EarningsBasePayroll[i]+=MaxTaxableEarnings*GetEfficientPopulation(j,Year,i,Country)

        EarningsBasePayroll[WM.Yclasses]+=EarningsBasePayroll[i]

    for i in xrange(WM.Yclasses):
        for j in xrange(WM.gens-WM.RetirementAge[Year,Country]):
            #print "j plug-in",j
            IK=Year_of_Retirement(j,Year,Country)
            IK=IK-WM.startworkingage
            if Year>= 0 and IK <0:
                IK=0

            WM.Pension_Replacement_Rate[j,IK,i,Country]=WM.Tax_Rate[3,0,IK,Country]-WM.Tax_Rate[3,1,IK,Country]*\
                    WM.Avg_Indexed_Earnings[IK,i,Country]
            if Country==0 and j <63:
                WM.Pension_Replacement_Rate[j,IK,i,Country]=0.0
            #print WM.Pension_Replacement_Rate[j,IK,i,Country]
            if WM.Pension_Replacement_Rate[j,IK,i,Country]<0.: #CHECK THIS AGAIN THERE'S SOME WEIRDNESS HERE, TOO
                WM.Pension_Replacement_Rate[j,IK,i,Country]=0.

            WM.Pension_Benefits_Ind[j,Year,i,Country]=WM.Pension_Replacement_Rate[j,IK,i,Country]*\
                    WM.Avg_Indexed_Earnings[IK,i,Country]*WM.YY_0[IK,Country]/WM.YY[IK,Country]
               
            WM.Pension_Benefits[Year,Country]+=WM.Pension_Benefits_Ind[j-WM.startworkingage,Year,i,Country]*\
                    GetEfficientPopulation(j,Year,i,Country)

    if EarningsBasePayroll[WM.Yclasses] != 0.:
        WM.Agg_Pension_Tax_Rate[Year,Country]=(1.-WM.Mu1[Year,Country])*WM.Pension_Benefits[Year,Country]/\
                 EarningsBasePayroll[WM.Yclasses]
    #else:
        #print " PY=0, Country=", Country, ", Year=", Year

    WM.Agg_Marg_Pension_Tax_Rate[Year,WM.Yclasses,Country]=0.
    WM.Agg_Avg_Pension_Tax_Rate[Year,WM.Yclasses,Country]=0.
    for i in xrange(WM.Yclasses):
        WM.Agg_Marg_Pension_Tax_Rate[Year,i,Country]=0.
        WM.Agg_Avg_Pension_Tax_Rate[Year,i,Country]=0.
        for j in xrange(WM.startworkingage,WM.RetirementAge[Year,Country]-1):
            if EarningsBasePayrollIND[j,i] <= MaxTaxableEarnings:
                if i < WM.Yclasses:
                    WM.Marg_Pension_Tax_Rate[j,Year,i,Country]=WM.Marg_Pension_Tax_Rate[j,Year,i,Country]*(1-TDAMP)+\
                            TDAMP*WM.Agg_Pension_Tax_Rate[Year,Country]
                else:
                    Marg_Pension_Tax_Rate[j,Year,i,Country]=0.
                WM.Avg_Pension_Tax_Rate[j,Year,i,Country]=WM.Avg_Pension_Tax_Rate[j,Year,i,Country]*(1.-TDAMP)+\
                        TDAMP*WM.Agg_Pension_Tax_Rate[Year,Country]
                WM.Agg_Marg_Pension_Tax_Rate[Year,i,Country]+=WM.Marg_Pension_Tax_Rate[j,Year,i,Country]*\
                        EarningsBasePayrollIND[j,i]*GetEfficientPopulation(j,Year,i,Country)
                WM.Agg_Avg_Pension_Tax_Rate[Year,i,Country]+=WM.Avg_Pension_Tax_Rate[j,Year,i,Country]*\
                        EarningsBasePayrollIND[j,i]*GetEfficientPopulation(j,Year,i,Country)
            else:
                WM.Marg_Pension_Tax_Rate[j,Year,i,Country]=0.
                if EarningsBasePayrollIND[j,i]>0:
                    WM.Avg_Pension_Tax_Rate[j,Year,i,Country]=WM.Avg_Pension_Tax_Rate[j,Year,i,Country]*(1.-TDAM)+\
                            TDAM*WM.Agg_Pension_Tax_Rate[Year,Country]*MaxTaxableEarnings/EarningsBasePayrollIND[j,i]
                else:
                    WM.Avg_Pension_Tax_Rate[j,Year,i,Country]=0.
                WM.Agg_Marg_Pension_Tax_Rate[Year,i,Country]+=WM.Marg_Pension_Tax_Rate[j,Year,i,Country]*\
                        MaxTaxableEarnings*GetEfficientPopulation(j,Year,i,Country)
                WM.Agg_Avg_Pension_Tax_Rate[Year,i,Country]+=WM.Avg_Pension_Tax_Rate[j,Year,i,Country]*MaxTaxableEarnings*\
                        GetEfficientPopulation(j,Year,i,Country)

            if WM.Agg_Pension_Tax_Rate[Year,Country]==0.:
                WM.Marg_Pension_Tax_Rate[j,Year,i,Country]=0.
                WM.Avg_Pension_Tax_Rate[j,Year,i,Country]=0.
                WM.Agg_Marg_Pension_Tax_Rate[Year,i,Country]=0.
                WM.Agg_Avg_Pension_Tax_Rate[Year,i,Country]=0.

        WM.Agg_Marg_Pension_Tax_Rate[Year,WM.Yclasses,Country]+=WM.Agg_Marg_Pension_Tax_Rate[Year,i,Country]
        WM.Agg_Avg_Pension_Tax_Rate[Year,WM.Yclasses,Country]+=WM.Agg_Avg_Pension_Tax_Rate[Year,i,Country]
        WM.Agg_Marg_Pension_Tax_Rate[Year,i,Country]=WM.Agg_Marg_Pension_Tax_Rate[Year,i,Country]/\
                EarningsBasePayroll[i]
        WM.Agg_Avg_Pension_Tax_Rate[Year,i,Country]=WM.Agg_Avg_Pension_Tax_Rate[Year,i,Country]/\
                EarningsBasePayroll[i]

    WM.Agg_Marg_Pension_Tax_Rate[Year,WM.Yclasses,Country]=WM.Agg_Marg_Pension_Tax_Rate[Year,WM.Yclasses,Country]/\
            EarningsBasePayroll[WM.Yclasses]
    WM.Agg_Avg_Pension_Tax_Rate[Year,WM.Yclasses,Country]=WM.Agg_Avg_Pension_Tax_Rate[Year,WM.Yclasses,Country]/\
            EarningsBasePayroll[WM.Yclasses]
    

    return None

def GetWage(Gen,Year,Yclass,Country): #3125-3133
    '''
    Computes the wage profile

    Inputs:
        Gen-
        Year-
        Yclass-
        Country-
    Outputs:
        Get_Wage-

    '''
    Gen-=WM.startworkingage #This is for the purposes of converting the age into a python compatable index

    Get_Wage=WM.Wage_Index_Class[Year,Yclass,Country]*WM.Age_Efficiency[Gen,Year,Yclass,Country]

    return Get_Wage

def ShadowWage(Gen, Year, Yclass, Country): #Fortran 1905-1942
    '''
    Computes the shadow wage (Or the opprotunity cost of leisure)
    Inputs:
    Gen-
    Year-
    Yclass
    Country-

    Outputs:
    Changes 

    '''

    IK=None
    KK=None

    Damp_Shw=.6

    for j in range (WM.gens-Gen): #(Gen,WM.gens)
        IK=GetYearBecomingJ(Year,Gen,j)
        if j>= WM.RetirementAge[IK,Country]:
            break
        if WM.Leisure[j,IK,Yclass,Country] > WM.Hours:
            break
        WM.Shadow_Wage[j,IK,Yclass,Country]*=Damp_Shw

    if j >= WM.RetirementAge[IK,Country] or WM.Leisure[j,IK,Yclass,Country] > WM.Hours:
        KK=j

        for j in range(WM.gens-KK): #(KK,WM.gens)
            IK=GetYearBecomingJ(Year,Gen,j)
            Oldsh=WM.Shadow_Wage[j,IK,Yclass,Country]
            WM.Shadow_Wage[j,IK,Yclass,Country]=WM.Consump_Price[IK,Country]*WM.Alp/\
                    (WM.Marg_Wage_Tax_Rate[j,IK,Yclass,Country]-WM.Marg_Pension_Tax_Rate[j,IK,Yclass,Country]-\
                    WM.Marg_Health_Tax_Rate[j,IK,Yclass,Country]-WM.Marg_Disability_Tax_Rate[j,IK,Yclass,Country])*\
                    (WM.Hours/WM.Consump[j,IK,Yclass,Country])**(-1./WM.Rho)-GetWage(j,IK,Yclass,Country)

            WM.Shadow_Wage[j,IK,Yclass,Country]= WM.Shadow_Wage[j,IK,Yclass,Country]*Damp_Shw+(1.-Damp_Shw)*Oldsh
           
    return None

def AvgIndexedEarnings(Year,Country): #Fortran 2921-2942
    '''
    Inputs:
    Year-The year you wish to average on
    Country-Choosing a nation

    Outputs:
    Fills in the array.

    '''

    for i in xrange(WM.Yclasses):
        WM.Avg_Indexed_Earnings[Year,i,Country]=0.
        for j in xrange(WM.startworkingage, WM.RetirementAge[Year,Country]-1):
            IK=GetYearBecomingJ(Year,WM.RetirementAge[Year,Country],j)
            if Year >= 0 and IK<0 :
                IK=0
            WM.Avg_Indexed_Earnings[Year,i,Country]+=WM.Hours-WM.Leisure[j,IK,i,Country]*Sum_Wage(j,IK,i,Country)

        WM.Avg_Indexed_Earnings[Year,i,Country]=WM.Avg_Indexed_Earnings[Year,i,Country]/\
                (WM.RetirementAge[Year,Country]-WM.startworkingage)

    return None


def GoodsMarket(Year): #Fortran 2979-3047
    '''
    Inputs:
    Year-The Current year

    Outputs:

    '''
    Goods_Market=0.
    #DD_World[year]=0. #Not super necessary since the global arrays are declared to have entires as zero.
    #YY_World[year]=0.
    #Net_Endow_Invest=0.
    IPL1=GetYearBecomingJ(Year,0,1)

    for i in xrange(WM.countries):
        #Aggregate Supply
        #WM.YY[Year,i]=0.
        #WM.National_Income[Year,i]=0.
        #WM.GDP[Year,i]=0.

        for j in xrange(WM.Yclasses):
            WM.YY[Year,i]=WM.YY[Year,i]+WM.Labor[Year,j,i]*WM.Wage_Index_Class[Year,j,i]
            #print "Wage_Index Class in year: ", Year, ":",WM.Wage_Index_Class[Year,j,i]

        WM.YY[Year,i]=WM.YY[Year,i]+WM.Capital[Year,i]*WM.R[Year,i]+WM.Agg_Assets_Migrants[Year,WM.Yclasses,i]

        WM.GDP[Year,i]=WM.YY[Year,i]+WM.Country_Endow_Share[Year,i]*WM.Endowment[Year]

        if (WM.IRUN==0):
            WM.YY_0[Year,i]=WM.YY[Year,i]

        #Calculating the net change in investment

        WM.Invest[Year,i]=(1.+WM.Tech)*(1.+WM.NPOP[IPL1,WM.firstcountry])*WM.Capital[IPL1,i]-(1-WM.Del)*WM.Capital[Year,i]

        #C,I,G, altogether gets to Y, keep the total growing, watch the Economy Fly! (Aggregate demand)
        WM.DD[Year,i]=WM.CC[Year,i]+WM.Invest[Year,i]+WM.Government_Expenditures[Year,i]
        WM.DD_World[Year]+=WM.DD[Year,i]
        WM.YY_World[Year]+=WM.YY[Year,i]

    WM.Net_Endow_Invest[Year]=(1.+WM.Tech)*(1.+WM.NPOP[IPL1,WM.firstcountry])*WM.PVEndowment[Year]

    WM.Tot_Gov_Endow_Revenue[Year]=WM.Tot_Gov_Endow_Share[Year]*WM.Endowment[Year]

    #Now Getting the Markets to Balance
    #Add net investment in the asset to world demand

    WM.YY_World[Year]+=WM.PVEndowment[Year]*WM.RG[Year]+WM.Tot_Gov_Endow_Revenue[Year]

    #Add Net investment in the asset to world demand
    
    WM.DD_World[Year]+=WM.Net_Endow_Invest[Year]

    #CHECK HERE FOR THE Market to Clear

    if (abs(WM.YY_World[Year]-WM.DD_World[Year])<WM.tol):
        Goods_Market=1
        print "Good Year!", Year
    else :
        print "Bad Year, this much excess supply:", (WM.YY_World[Year]-WM.DD_World[Year])

    return Goods_Market


def GetEfficientPopulation(Gen,Year,Yclass,Country):
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''
    EfficientPopulation=0

    EfficientPopulation=(1.+WM.Tech)**(WM.startworkingage-Gen)*WM.POP[Gen,Year,Yclass,Country]/WM.POP[WM.startworkingage,Year,-1,0] 
    #0 is put in to indicate the first country, namely, the United States of America. MURICA
            
    return EfficientPopulation

def Initialize_Scenario():  #1086-1078
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''

#Read DATA function

def Population_Development():   #1085-1301
    #FORTRAN PURPOSE: This subroutine calls READ_DATA, and then gets generates the remaining demographic data needed
    #PURPOSE: This subroutine calls READ_DATA, and then gets generates the remaining demographic data needed

    def plotstuff(name, index, years):
        for y in range(len(years)):
          yeartograph = years[y]
          if yeartograph - 2008 < WM.Years+1:
            #IMPORTANT! Plots the sum of the skill classes for each year

            plt.plot(range(WM.gens+1), WM.POP[:-1, yeartograph-2008, -1 , index])
          else:
            print "\nERROR: WE HAVE ONLY SIMULATED UP TO THE YEAR", WM.Years+2008,\
                    "AND YOU REQUESTED TO PLOT THE YEAR", yeartograph
            print"*SEE plot_population_distribution IN class Region(object)*\n"
            time.sleep(15)

        plt.title(str(name + " Population Distribution"))
        plt.legend(years)
        plt.show()
        plt.clf()


        plt.title(str(name+ " Total Population Beginning in 2008"))
        plt.plot(range(1,yeartograph-2008), WM.POP[-1, 1:yeartograph-2008, -1 , index])
        plt.show()

    #Assigning Values to Survival_Probability
    for r in range(WM.firstcountry, WM.lastcountry+1):
        for t in xrange(WM.Years+1):
            for g in range(WM.gens+1):
                if t == 1:
                    WM.Survival_Probability[g, t, r] = 1.
                else:
                    if g == 0:
                        WM.Survival_Probability[0, t, r] = 1-WM.Mortality_rates[0, t, r]
                    else:
                        WM.Survival_Probability[g, t, r] = WM.Survival_Probability[g-1, t-1, r] * (1-WM.Mortality_rates[g, t, r])


    #Child Mortality in India
    WM.Mortality_rates[:WM.startworkingage, 0, 5] = 0.01
    for t in range(1,51):
        WM.Mortality_rates[:WM.startworkingage, t, 5] = WM.Mortality_rates[:WM.startworkingage, t-1, 5] - 0.007/50
    WM.Mortality_rates[:WM.startworkingage, 51:-1, 5] = 0.003
    

    #Child Mortality in Russia
    WM.Mortality_rates[:WM.startworkingage-2, 0, 6] = 0.0022
    for t in range(1,43):
        WM.Mortality_rates[:WM.startworkingage-2, t-1, 6] = WM.Mortality_rates[:WM.startworkingage-2, t-1, 6] - 0.00095/42.
    WM.Mortality_rates[:WM.startworkingage-2, 43:-1, 6] = 0.0011
    

    #ASSIGNING VALUES TO POP MATRIX
    #NOTE: THIS CODE OPERATES UNDER THE ASSUMPTION OF A POPULATION MATIRX OF SIZE (GENS+2, YEARS+2, YCLASSES+2, REGIONS), with the extra rows in the first 3 dimensions containing sums. CHANGE THIS once fortran code works
    for t in xrange(1, WM.Years+1):
        #If the current year t is less than or equal to the greatest year for which we have data (Transyear), use the current year t as the index i for pulling data
        if t <= WM.TransYear:
            ExogenousGrowth = 1
                
        #If the current year t is beyond the greatest year for which we have data (Transyear), use the most recent year as the index i for pulling data
        elif t > WM.TransYear:
            ExogenousGrowth = (1+WM.Exog_Npop)**(t-WM.TransYear)

        #Gets the number of immigrants for each region and each class
        immigrants = np.einsum('kr,gr -> gkr', WM.C_Share_Migrants[t, :, :], WM.Migrants[1:,:] * WM.Migration_Scale[1:, t, :])

        #Gets the new population by adding immigrants and the previous native population (without its oldest generation)
        #NOTE: The matrix newpopulation is shifted by one year. Thus it dosn't have any entries for newborns, so its lowest index contains 1 year-olds
        newpopulation = WM.POP[:-2,t-1,:-1, :] + immigrants*ExogenousGrowth

        #Gets the surviving fraction of the population and stores it in the population timepath for the current year t
        WM.POP[1:-1,t,:-1, :] = np.einsum('gr,gkr -> gkr', (1-WM.Mortality_rates[1:,t, :]), newpopulation)

        #If we aren't past the transition year yet...
        if t <= WM.TransYear:
            #Gets the number of newborns for each region and class using fertility rates and the fertile population, taking into account infant mortality
            #Syntax for C = np.einsum('fr,fkr -> kr', A, B): A is fxr dimensions, B is fxkxr dimensions, and C is kxr dimensions. f=fertile years, k=number of skill classes, r=number of regions/countries
            newborns = np.einsum('fr,fkr -> kr', WM.Fertility[:,t, :], WM.POP[WM.startfertilityage:WM.lastfertilityage+1,t,:WM.Yclasses, :] * (1-WM.Mortality_rates[0,t, :]))

            #Adds the newborns to any immigrant newborns that came this year
            WM.POP[0,t,:-1, :] += newborns

        #If we are past the transition year...
        else:
            #Population growth determined exogenously
            WM.POP[0,t,:, :] = (1+WM.Exog_Npop)*WM.POP[0,t-1,:, :]


    #Getting the sum of skill classes and generations respectively
    WM.POP[:-1,:-1,-1,:] = np.sum(WM.POP[:-1,:-1,:,:], axis=2)
    WM.POP[-1,:-1,-1,:] = np.sum(WM.POP[:-1,1:,-1,:], axis=0)


    #Calculating population growth rates
    WM.NPOP[-1,:] = WM.Exog_Npop
    WM.NPOP[0,:] = WM.POP[22, 0, -1, :]/WM.POP[WM.startworkingage, 0, -1, :] - 1
    for t in xrange(1, WM.Years+1):
        WM.NPOP[t, :] = WM.POP[WM.startworkingage, t, -1, :]/WM.POP[WM.startworkingage, t-1, -1, :] - 1


    #Calculating fertility rates endogenously from trans_year+1 to 300 
    for t in xrange(WM.TransYear+1, WM.Years+1):
        WM.Fertility[:, t, :] = (1+WM.Exog_Npop) * WM.Fertility[:, t-1, :] * WM.POP[WM.startfertilityage:WM.lastfertilityage+1, t-1, -1, :]/WM.POP[WM.startfertilityage:WM.lastfertilityage+1, t-1, -1, :]
    WM.Total_Fertility[WM.TransYear+1:WM.Years+1,:] = np.sum(WM.Fertility[:,WM.TransYear+1:WM.Years+1,:], axis=0)


    #CALL POPULATION PRINT SUMMARY


    #Adopt population in -1 from YEARS
    if WM.IRun == 0:
        WM.POP[:, -1, :, :] = WM.POP[:, WM.Years+1, :, :]        
        WM.Mortality_rates[:, -1, :] = WM.Mortality_rates[:, WM.Years+1, :]
        WM.Survival_Probability[:, -1, :] = WM.Survival_Probability[:, WM.Years+1, :]
        WM.Fertility_ss = WM.Fertility[:, WM.Years+1, :]
        WM.TotalFertilityss = WM.Total_Fertility[WM.Years+1, :]


    #Getting Efficient Population
    for k in range(WM.Yclasses):
        for g in range(WM.startworkingage, WM.gens+1):
            WM.POP_Efficient += WM.POP[g, :, k, :]*(1+WM.Tech)**(WM.startworkingage-g)


    WM.POP[-1, :, :, :] = np.sum(WM.POP, axis=0)
    WM.POP[:,-1, :, :] = WM.POP[:,-2, :, :]

    #plotstuff("USA", 0, np.array([2008, 2058, 2308]))

def Get_Aggregate_Variables(Year): #1947-2067                                                               
    '''
    Damps Capital and labor for iterating in transition.
    Inputs:
    Year-
    Outputs:
    Aggregate Variables
    '''
    OldLabor=np.zeros((WM.Yclasses,WM.countries))
    IPL1=GetYearBecomingJ(Year,0,1)
    WM.Agg_Assets_World[Year]=0.
    WM.Trade_Balance_World[Year]=0.
    HHB=0. #HHB is the government debt for the whole world except the US
    HHK=0. #HHK is the sum of capital deman of the open economies

    OldCap=WM.Capital[Year,WM.firstcountry]
    #Total amounts of consumption/labor etc. is found by summing over individuals

    for i in xrange(WM.countries):
        WM.CC[Year,i]=Sum_Individual_Variables(Year,1,i)

        for j in xrange(WM.Yclasses):
            OldLabor[j,i]=WM.Labor[Year,j,i]

        WM.Labor[Year,WM.Yclasses,i]=Sum_Individual_Variables(Year,2,i)

        for j in xrange(WM.Yclasses):
            WM.Labor[Year,j,i]=OldLabor[j,i]+WM.Damp*(WM.Labor[Year,j,i]-OldLabor[j,i])

        WM.Agg_Assets[Year,i]=Sum_Individual_Variables(Year,3,i)
        WM.Agg_Assets_World[Year]+=WM.Agg_Assets[Year,i]


    #Note that in the steady state, there is no endowment

    WM.PVEndowment[-1]=0

    #Then, we calculate the prsent value of the endowment flow to households
    #Using the no arbitrage condition between owning an oil well and owning
    #a unit of capital, we can find the value of oil assets owned by
    #the private sector

    if Year<WM.Years+1:
        WM.PVEndowment[Year]=0
        for q in xrange(0,WM.Years-Year):
            WM.PVEndowment[Year]=((1/(1+WM.RG[WM.Years-q]))*(WM.PVEndowment[Year]))+\
                    (1.0-WM.Tot_Gov_Endow_Share[(WM.Years-q)])*WM.Endowment[WM.Years-q]

    if (WM.firstcountry != WM.lastcountry):
        for i in xrange(1,WM.lastcountry+1):
            HHK+=WM.Capital[Year,i]
            HHB+=WM.Debt[Year,i]

    #Capital in the US iS NOT calculated using the interest FOC like everybody else.
    #Rather, it is calculated by taking total world assets and subtracting off
    #capital in other countries, debt in all countries, and privately held oil assets (PVEndowment)
    WM.Capital[Year,WM.firstcountry]=WM.Agg_Assets_World[Year]-HHK-HHB-\
            WM.Debt[Year,WM.firstcountry]-WM.PVEndowment[Year]
    #This damps the guess of capital in the US
    WM.Capital[Year,WM.firstcountry]=OldCap+WM.Dampk*(WM.Capital[Year,WM.firstcountry]-OldCap)

    #Sum Transfers to Foreign Regions
    for i in xrange(WM.countries):
        WM.TRF[Year,i]=0.
        for j in xrange(WM.countries):
            if i != j:
                for k in xrange(WM.Yclasses):
                    for l in xrange(WM.WorkToDeath):
                        WM.TRF[Year,i]+=WM.H_Transfer[l,Year,k,j,i]*GetEfficientPopulation(l,Year,k,j)-\
                                WM.H_Transfer[l,Year,k,i,j]*GetEfficientPopulation(l,Year,k,i)

    #Gets foreign assets and trade Balance
    if WM.firstcountry != WM.lastcountry:
        for i in xrange(WM.countries):
            WM.Foreign_Assets[Year,i]=WM.Agg_Assets[Year,i]-WM.Debt[Year,i]-WM.Capital[Year,i]

            if Year==301 or Year==WM.Years:
                WM.TradeBalance[Year,Country]=((1.+WM.Tech)*(1.+WM.NPOP[Year,WM.firstcountry])-(1.+WM.RG[Year]))*WM.Foreign_Assets[Year,i]+\
                        WM.TRF[Year-1,i]
            if (Year<300 or Year>WM.First_Solution_Year):
                WM.Trade_Balance[Year-1,i]=(1.+WM.Tech)*(1.+WM.NPOP[Year,WM.firstcountry])*\
                        WM.Foreign_Assets[Year,i]-(1.+WM.RG[Year-1])*WM.Foreign_Assets[Year-1,i]+WM.TRF[Year-1,i]
            WM.Trade_Balance_World[Year]+=WM.Trade_Balance[Year,i]
            WM.Foreign_Assets_World[Year]+=WM.Foreign_Assets[Year,i]

    WM.Foreign_Assets_World-=WM.PVEndowment[Year]

    WM.Trade_Balance_World[Year]-=((1.+WM.Tech)*(1.+WM.NPOP[IPL1,WM.firstcountry])*WM.PVEndowment[IPL1]-WM.PVEndowment[Year])+WM.PVEndowment[Year]*WM.RG[Year]



    return None

def Sum_Individual_Variables(Year, JD, Country): #2074-2158                                                          
    '''
    PURPOSE: General routine to compute aggregate consumption (JD = 1), labor (JD = 2), assets (JD = 3), and wages (JD = 4)
    Inputs:
    Year-
    JD-
    Country-

    Outputs:
    Fills in the array
    '''
    IndividualVariables=0.

    if JD==3:
        WM.Agg_Assets_Migrants[Year,WM.Yclasses,Country]=0.

    for i in xrange(WM.Yclasses):
        
        if JD==2:
            WM.Labor[Year,i,Country]=0.
        if JD==3:
            WM.Agg_Assets_Migrants[Year,i,Country]=0.
            WM.Agg_Assets_For_Bequests[Year,i,Country]=0.

        for j in xrange(WM.WorkToDeath):

            IK=GetYearBecomingJ(Year,j,j+1)

            if JD==1:
                IndividualVariables+=(WM.Consump[j,Year,i,Country]+KIDS(j,Year,i,Country)*WM.Consump_Kids[j,Year,i,Country])*\
                        GetEfficientPopulation(j,Year,i,Country)

            if JD==2:
                if WM.Leisure[j,Year,i,Country]<=WM.Hours or j >= WM.RetirementAge[Year,Country]:
                    WM.Labor[Year,i,Country]+=WM.Age_Efficiency[j,Year,i,Country]*(WM.Hours-WM.Leisure[j,Year,i,Country])*\
                            GetEfficientPopulation(j,Year,i,Country)

            if JD==3:
                WM.Agg_Assets_For_Bequests[Year,i,Country]+=WM.Assets[j,Year,i,Country]*(1.+WM.Tech)**(WM.startworkingage-j)*\
                        WM.POP[j,Year,i,Country]*WM.Mortality_rates[j,Year,Country]/(1.-WM.Mortality_rates[j,Year,Country])

                IndividualVariables+=WM.Assets[j,Year,i,Country]*GetEfficientPopulation(j,Year,i,Country)/\
                        (1.-WM.Mortality_rates[j,Year,Country])

                if Year==WM.Years+1:
                    WM.Agg_Assets_Migrants[Year,i,Country]+=WM.Assets[j,Year,i,Country]*WM.Migrants[j,Country]*\
                            WM.C_Share_Migrants[Year,i,Country]*WM.Migration_Scale[j,Year,Country]*\
                            (1.-WM.Mortality_rates[j,Year,Country])*\
                            (1.+WM.Exog_Npop)**(WM.Years-(WM.Trans_Year-WM.First_Year+2000)+1)*\
                            (1.+WM.Tech)**(22-j)/WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]

                elif Year<WM.Years+1 and Year < WM.TransYear-WM.First_Year+2000:
                    WM.Agg_Assets_Migrants[Year,i,Country]+=WM.Assets[j,IK,i,Country]*WM.C_Share_Migrants[IK,i,Country]*\
                            WM.Migration_Scale[j,IK,Country]*WM.Migrants[j,Country]*(1.-WM.Mortality_rates[j,Year,Country])*\
                            (1.+WM.Tech)**(22-j)/WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]

                elif Year>= (WM.TransYear-WM.First_Year+2000) and Year<= WM.Years:
                    WM.Agg_Assets_Migrants[Year,i,Country]+=WM.Assets[j,IK,i,Country]*WM.C_Share_Migrants[IK,i,Country]*\
                            WM.Migration_Scale[j,IK,Country]*WM.Migrants[j,Country]*(1.-(j,IK,Country))*\
                            (1.+WM.Exog_Npop)**(IK-(WM.Trans_Year-WM.First_Year+2000))*(1.+WM.Tech)**(22-j)/\
                            WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]


            if JD==4:
                if WM.Leisure[j,Year,i,Country] <= WM.Hours or j >= WM.RetirementAge[Year,Country]:
                        IndividualVariables+=((WM.Hours-WM.Leisure[j,Year,i,Country])*GetWage(j,Year,i,Country))**2 *\
                                GetEfficientPopulation(j,Year,i,Country)

        if JD==2:
            IndividualVariables+=WM.Labor[Year,i,Country]
        if JD==3:
            WM.Agg_Assets_Migrants[Year,WM.Yclasses,Country]+=WM.Agg_Assets_Migrants[Year,i,Country]


    return IndividualVariables

def Get_Marginal_Products(Year, Country): #2165-2191                                                                   
    '''
    PURPOSE: Compute wages for all countries, capital productivity only for the US from production function.
    Inputs:
    Year-
    Country-

    Outputs:
    Nothing direct, it just changes some of the global arrays    
    '''
    for i in xrange(WM.Yclasses):
        HLabor=1.
        for j in xrange(WM.Yclasses):
            if i != j:
                HLabor=HLabor*WM.Labor[Year,j,Country]**WM.Beta[j]

        WM.Wage_Index_Class[Year,i,Country]=WM.Beta[i]*WM.AP[Year]*\
                WM.Capital[Year,Country]**WM.Alpha*WM.Labor[Year,i,Country]**(WM.Beta[i]-1.)*HLabor

    if Country==WM.firstcountry:
        HLabor=1.
        for i in xrange(WM.Yclasses):
            HLabor=HLabor*WM.Labor[Year,i,Country]**WM.Beta[i]

        WM.R[Year,Country]=WM.Alpha*WM.AP[Year]*WM.Capital[Year,Country]**(WM.Alpha-1.)*HLabor


    return None

def Health_System(Year, Country): #2444-2639                                                                               
    #PURPOSE: Calculate components of the country's health system.
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''
    TDAMP=1.
    TDAM=1.

    WM.Health_Benefits[Year,Country]=0.
    Earnings_Base_Payroll_Ind=np.zeros((WM.WorkToDeath,WM.Yclasses))
    Earnings_Base_Payroll=np.zeros(WM.Yclasses+1)

    Get_Avg_Labor_Earnings(Year,Country)

    Max_Taxable_Earnings=WM.Contribution_Ceiling[Year,Country]*WM.Avg_Labor_Earnings[Year,Country]

    for a in xrange(WM.gens):
        if Country <= 2 or Country==5 or Country==6:
            if Year>=0 and Year<=25:
                WM.Health_Benefits_Ind[a,Year,Country]=WM.Health_Benefits_Ind[a,0,Country]*(WM.YY_0[Year,Country]/\
                        WM.YY[0,Country]*WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]/\
                        WM.POP[WM.startworkingage,0,WM.Yclasses,WM.firstcountry]*\
                        WM.POP_Efficient[0,Country]/WM.POP_Efficient[Year,Country])*1.02**Year
            elif Year>=26 and Year<=35:
                WM.Health_Benefits_Ind[a,Year,Country]=WM.Health_Benefits_Ind[a,25,Country]*\
                        (WM.YY_0[Year,Country]/WM.YY_0[25,Country]*\
                        WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]/\
                        WM.POP[WM.startworkingage,25,WM.Yclasses,WM.firstcountry]*WM.POP_Efficient[25,Country]/\
                        WM.POP_Efficient[Year,Country])*1.01**(Year-25)

            elif Year >=36:
                WM.Health_Benefits_Ind[a,Year,Country]=WM.Health_Benefits_Ind[a,35,Country]*\
                        (WM.YY_0[Year,Country]/WM.YY_0[35,Country]*WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]/\
                        WM.POP[WM.startworkingage,35,WM.Yclasses,WM.firstcountry]*WM.POP_Efficient[35,Country]/\
                        WM.POP_Efficient[Year,Country])
        elif Country>=3 and Country <=4:
            if Year>=0 and Year <=40:
                WM.Health_Benefits_Ind[a,Year,Country]=WM.Health_Benefits_Ind[a,0]*\
                        (WM.YY_0[Year,Country]/WM.YY[0,Country]*WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]/\
                        WM.POP[WM.startworkingage,0,WM.Yclasses,WM.firstcountry]*WM.Pop_Efficient[0,Country]/\
                        WM.POP_Efficient[Year,Country])*1.04**Year
            elif Year >=41:
                WM.Health_Benefits_Ind[a,Year,Country]=WM.Health_Benefits_Ind[a,41,Country]*\
                        (WM.YY_[Year,Country]/WM.YY_0[41,Country]*\
                        WM.POP[WM.startworkingage,Year,WM.Yclasses,WM.firstcountry]/\
                        WM.POP[WM.startworkingage,41,WM.Yclasses,WM.firstcountry]*WM.POP_Efficient[41,Country]/\
                        WM.POP_Efficient[Year,Country])

        if Year==WM.Years+1:
            WM.Health_Benefits_Ind[a,Year,Country]=WM.Health_Benefits_Ind[a,Year,Country]

    for i in xrange(WM.Yclasses):
        for j in xrange(WM.WorkToDeath):
            WM.Health_Benefits[Year,Country]+=(WM.Health_Benefits_Ind[j,Year,Country]+Kids_Health_Benefits(j,Year,i,Country))*\
                    GetEfficientPopulation(j,Year,i,Country)

    for i in xrange(WM.Yclasses):
        Earnings_Base_Payroll[i]=0.
        for j in xrange(WM.startworkingage,WM.RetirementAge[Year,Country]):
            Earnings_Base_Payroll_Ind[j,i]=Sum_Wage(j,Year,i,Country)*(WM.Hours-WM.Leisure[j,Year,i,Country])

            if Earnings_Base_Payroll_Ind[j,i]<0.:
                Earnings_Bas_Payroll_Ind[j,i]=0.
            if Country==0:
                Earnings_Base_Payroll[i]+=Earnings_Base_Payroll_Ind[j,i]*GetEfficientPopulation(j,Year,i,Country)
            else:
                if Earnings_Base_Payroll_Ind[j,i] <= Max_Taxable_Earnings:
                    Earnings_Base_Payroll[i]+=Earnings_Base_Payroll_Ind[j,i]*GetEfficientPopulation(j,Year,i,Country)
                else:
                    Earnings_Base_Payroll[i]+=Max_Taxable_Earnings*GetEfficientPopulation(j,Year,i,Country)

        Earnings_Base_Payroll[WM.Yclasses]+=Earnings_Base_Payroll[i]

    if Earnings_Base_Payroll[WM.Yclasses] != 0:
        WM.Agg_Health_Tax_Rate[Year,Country]=(1.-WM.Mu2_gov[Year,Country]+WM.Mu2_gov[Year,Country])*\
                WM.Health_Benefits[Year,Country]/Earnings_Base_Payroll[WM.Yclasses]

    else:
        print "PY=0, Country=",Country,", Year= ", Year

    WM.Agg_Marg_Health_Tax_Rate[Year,WM.Yclasses,Country]=0.
    WM.Agg_Avg_Health_Tax_Rate[Year,WM.Yclasses,Country]=0.
    for i in xrange(WM.Yclasses):
        WM.Agg_Marg_Health_Tax_Rate[Year,i,Country]=0.
        WM.Agg_Avg_Health_Tax_Rate[Year,i,Country]=0.
        for j in xrange(WM.RetirementAge[Year,Country]-WM.startworkingage):
            if Country==0:
                WM.Marg_Health_Tax_Rate[j,Year,i,Country]=WM.Marg_Health_Tax_Rate[j,Year,i,Country]*(1-TDAMP)+\
                        TDAMP*WM.Agg_Health_Tax_Rate[Year,Country]
                WM.Avg_Health_Tax_Rate[j,Year,i,Country]=WM.Avg_Health_Tax_Rate[j,Year,i,Country]*(1.-TDAMP)+\
                        TDAMP*WM.Agg_Health_Tax_Rate[Year,Country]
                WM.Agg_Marg_Health_Tax_Rate[Year,i,Country]+=WM.Marg_Health_Tax_Rate[j,Year,i,Country]*\
                        Earnings_Base_Payroll_Ind[j,i]*GetEfficientPopulation(j,Year,i,Country)
                WM.Agg_Avg_Health_Tax_Rate[Year,i,Country]+=WM.Marg_Health_Tax_Rate[j,Year,i,Country]*\
                        Earnings_Base_Payroll_Ind[j,i]*GetEfficientPopulation(j,Year,i,Country)
            else:
                if Earnings_Base_Payroll_Ind[j,i]<=Max_Taxable_Earnings:
                    if i < WM.Yclasses:
                        WM.Marg_Health_Tax_Rate[j,Year,i,Country]=WM.Marg_Health_Tax_Rate[j,Year,i,Country]*(1.-TDAMP)+\
                                TDAMP*WM.Agg_Health_Tax_Rate[Year,Country]
                    else:
                        WM.Marg_Health_Tax_Rate[j,Year,i,Country]=0.
                    WM.Avg_Health_Tax_Rate[j,Year,i,Country]=WM.Avg_Health_Tax_Rate[j,Year,i,Country]*\
                            (1.-TDAMP)+TDAMP*WM.Agg_Health_Tax_Rate[Year,Country]
                    WM.Agg_Marg_Health_Tax_Rate[Year,i,Country]+=WM.Marg_Health_Tax_Rate[j,Year,i,Country]*\
                            Earnings_Base_Payroll_Ind[j,i]*GetEfficientPopulation(j,Year,i,Country)
                    WM.Agg_Avg_Health_Tax_Rate[Year,i,Country]+=WM.Avg_Health_Tax_Rate[j,Year,i,Country]*\
                            Earnings_Base_Payroll_Ind[j,i]*GetEfficientPopulation(j,Year,i,Country)
                else:
                    WM.Marg_Health_Tax_Rate[j,Year,i,Country]=0.
                    if Earnings_Base_Payroll_Ind[j,i]>0:
                        WM.Avg_Health_Tax_Rate[j,Year,i,Country]=WM.Avg_Health_Tax_Rate[j,Year,i,Country]*(1.-TDAM)+\
                                TDAM*WM.Agg_Health_Tax_Rate[Year,Country]*Max_Taxable_Earnings/Earnings_Base_Payroll_Ind[j,i]
                    else:
                        WM.Avg_Health_Tax_Rate[j,Year,i,Country]=0.
                        WM.Agg_Marg_Health_Tax_Rate[Year,i,Country]+=WM.Marg_Health_Tax_Rate[j,Year,i,Country]*\
                            Max_Taxable_Earnings*GetEfficientPopulation(j,Year,i,Country)
                        WM.Agg_Avg_Health_Tax_Rate[Year,i,Country]+=WM.Avg_Health_Tax_Rate[j,Year,i,Country]*\
                                Max_Taxable_Earnings*GetEfficientPopulation(j,Year,i,Country)
            
            if WM.Agg_Health_Tax_Rate[Year,Country]==0.:
                WM.Marg_Health_Tax_Rate[j,Year,i,Country]=0.
                WM.Avg_Health_Tax_Rate[j,Year,i,Country]=0.
                WM.Agg_Marg_Health_Tax_Rate[Year,i,Country]=0.
                WM.Agg_Avg_Health_Tax_Rate[Year,i,Country]=0.

        WM.Agg_Marg_Health_Tax_Rate[Year,WM.Yclasses,Country]+=WM.Agg_Marg_Health_Tax_Rate[Year,i,Country]
        WM.Agg_Avg_Health_Tax_Rate[Year,WM.Yclasses,Country]+=WM.Agg_Avg_Health_Tax_Rate[Year,i,Country]
        WM.Agg_Marg_Health_Tax_Rate[Year,i,Country]=WM.Agg_Marg_Health_Tax_Rate[Year,i,Country]/Earnings_Base_Payroll[i]
        WM.Agg_Avg_Health_Tax_Rate[Year,i,Country]=WM.Agg_Avg_Health_Tax_Rate[Year,i,Country]/Earnings_Base_Payroll[i]

    WM.Agg_Marg_Health_Tax_Rate[Year,WM.Yclasses,Country]=WM.Agg_Marg_Health_Tax_Rate[Year,WM.Yclasses,Country]/\
            Earnings_Base_Payroll[WM.Yclasses]
    WM.Agg_Avg_Health_Tax_Rate[Year,WM.Yclasses,Country]=WM.Agg_Avg_Health_Tax_Rate[Year,WM.Yclasses,Country]/\
            Earnings_Base_Payroll[WM.Yclasses]

                




    return None


def Disability_Insurance(YEAR, COUNTRY): #2265-2771                                                                    
    #PURPOSE: Calculate components of the country's disability insurance system
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''
    return None
                                                        

def KIDS(gen, year, yclass, country):                                                                   
    #PURPOSE: Compute number of kids of parents age GEN in YEAR 
    Kids = 0

    if gen <= 46:
        UU = 1
    else:
        UU = gen - 45

    if gen - WM.startfertilityage < 20:
        MM = gen - WM.startfertilityage
    else:
        MM = 20

    if gen >= WM.startfertilityage and gen <= 65:
        if gen <= 45 and year > -1: Kids = WM.Fertility[gen-WM.startfertilityage, year, country]
        if gen <= 45 and year == -1: Kids = WM.Fertility_ss[gen-WM.startfertilityage, country]
        for j in range(UU, MM+1):
            if year > -1: 
                Kids += WM.POP[j, year, yclass, country] * \
                    WM.Fertility[gen-j-WM.startfertilityage, year-j, country] / WM.Total_Fertility[year-j, country] \
                    / WM.POP[gen, year, yclass, country]
            if year == -1: 
                Kids += WM.POP[j, year, yclass, country] * \
                        WM.Fertility_ss[gen-j-WM.startfertilityage, country] / WM.TotalFertilityss[country]\
                        / WM.POP[gen, year, yclass, country]

    return Kids


def Population_Summary():  #3547-3671                                                                           
    #PURPOSE: Print demograpic info   
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''
    return None


def Year_of_Retirement(Gen, Year, Country):  #3721-3750                                                                 
    #PURPOSE: GET_YEAR_OF_RETIREMENT is the year in which someone age GEN in year YEAR has retired.  
    #As in the other GET_YEAR routines, this returns year -1 or YEARS if a year beyond this is called for 
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''
    retire=0.
    IK=0.
    IM=0.
    IPL=0.
    IKL=0.

    if Gen >= WM.RetirementAge[Year,Country]:
        IM=Gen-1
        IK1=Year-1
        for k in xrange(IM, WM.startworkingage,-1):
            if IK<-1:
                IK=-1
            if WM.RetirementAge[IK1,Country]<=k:
                IK=IK-1
        retire=IK1+1
    else:
        IPL=Gen+1
        IK1=Year+1
        for k in xrange(IPL,WM.gens):
            retire=IKL-1
    if retire>WM.Years:
        retire=WM.Years
    if retire<-1:
        retire=-1
    if Year==-1:
        retire=-1
    if Year==WM.Years:
        retire=WM.Years

    return retire


def Output(YEAR, I1, I2, COUNTRY): #3757-3918                                                                              
    #PURPOSE: Write output.    
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''
    return None


def Initial_Equilibrium(): #3925-4042                                                                                   
    #PURPOSE: Write data of the initial equilibrium
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''
    return None


def Results():     #4049-4395                                                                                               
    #PURPOSE: Print Results 
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''
    return None


def Hicksian_EV(): #4400-4455 DO WE ALREADY HAVE THIS???
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''
    return None


def Intialize_Variables():
    #PURPOSE: Make initial guesses
    '''
    Inputs:
    Gen-
    Year-
    YClass-
    Country-

    Outputs:
    Fills in the array
    '''

    return None

def Get_Avg_Labor_Earnings(Year,Country): #Fortran 2950-2973
    WorkPop=0.0

    for j in xrange(WM.startworkingage,WM.RetirementAge[Year,Country]-1):
        for k in xrange(WM.Yclasses):
            WM.Avg_Labor_Earnings[Year,Country]+=(WM.Hours-WM.Leisure[j,Year,k,Country])*Sum_Wage(j,Year,k,Country)*\
                    GetEfficientPopulation(j,Year,k,Country)
        WorkPop+=WM.POP[j,Year,k+1,Country]/WM.POP[WM.startworkingage,Year,k+1,WM.firstcountry]

    WM.Avg_Labor_Earnings[Year,Country]=WM.Avg_Labor_Earnings[Year,Country]/WorkPop

    return WM.Avg_Labor_Earnings

def Kids_Health_Benefits(Gen,Year,YClass,Country):
    '''
    Inputs:
        Gen-
        Year-
        YClass-
        Country-
    Outputs:
        Returns Health Benefits for Kids

    '''
    KidsHealthBenefits=0.
    UU=0

    if Gen<=46:
        UU=1
    else:
        UU=Gen-45

    if Gen-WM.startfertilityage<20:
        MM=Gen-WM.startfertilityage

    else:
        MM=20

    if Gen<=22 or Gen>=65:
        if Gen<=45 and Year<WM.Years+1:
            KidsHealthBenefits=WM.Fertility[Gen,Year,Country]*WM.Health_Benefits_Ind[0,Year,Country]
        if Gen<=45 and Year==WM.Years+1:
            KidsHealthBenefits=WM.Fertility_ss[Gen,Country]*WM.Health_Benefits_Ind[0,Year,Country]
        for j in xrange(UU,MM):
            if Year<WM.Years+1:
                KidsHealthBenefits+=WM.POP[j,Year,YClass]*WM.Fertility[Gen-j,Year-j,Country]/WM.Total_Fertility[Year-j,Country]\
                        /WM.POP[Gen,Year,YClass,Country]*WM>Health_Benefits_Ind[j,Year,Country]
            if Year==WM.Years+1:
                KidsHealthBenefits+=WM.POP[j,Year,YClass,Country]*WM.Fertility_ss[Gen-j,Country]/\
                        WM.TotalFertility_ss[Country]/WM.POP[Gen,Year,YClass,Country]*WM.Health_Benefits_Ind[j,Year,Country]

    return KidsHealthBenefits

def Read_Data():
    return None


def PopDevelopment():

    #First bits of fortran code are redundant, since the data is already read in and the appropriate arrays are already zeros.

    if WM.IRUN==0:
        Read_Data()

    #First, we calculate the rates
    for i in xrange(WM.countries):
        for j in xrange(WM.Years+1):
            for k in xrange(WM.gens+1):
                if j==0:
                    WM.Survival_Probability[k,j,i]=1.0
                else:
                    if k==0:
                        WM.Survival_Probability[0,j,i]=(1.-WM.Mortality_rates[0,j,i])
                    else:
                        WM.Survival_Probability[k,j,i]=WM.Survival_Probability[k-1,j-1,i]*(1.-WM.Mortality_rates[k,j,i])

    #Children in India
    for i in xrange(WM.Years+1):
        for j in xrange(0,WM.lasteducation+1): #These go through the childhood years
            if i==0:
                WM.Mortality_rates[j,i,4]=.01
            if i>=1 and i<=50:
                WM.Mortality_rates[j,i,4]=WM.Mortality_rates[j,i-1,4]-(.007/50.)
            if i>50:
                WM.Mortality_rates[j,i,4]=0.003

    #Children in mother Russia
    for i in xrange(WM.Years+1):
        for j in xrange(0,18+1):
            if i==0:
                WM.Mortality_rates[j,i,5]=0.0022
            if i>=1 and i<=42:
                WM.Mortality_rates[j,i,5]=WM.Mortality_rates[j,i-1,5]-(.00095/42.)
            if i>42:
                WM.Mortality_rates[j,i,5]=.0011

    #Calculation population in each year.
    #This is done by given what was read in (year zero population, future, migration and fertility)
    for i in xrange(WM.countries):
        for j in xrange(1,WM.Years+1):
            for k in xrange(WM.gens,0,-1):
                if k>0:
                    for l in xrange(WM.Yclasses):
                        if j <= WM.TransYear:
                            WM.POP[k,j,l,i]=(WM.POP[k-1,j-1,l,i]+WM.Migration_Scale[k,j,i]*\
                                    WM.C_Share_Migrants[j,l,i]*WM.Migrants[k,i])*(1.-WM.Mortality_rates[k,j,i])
                        else:
                            WM.POP[k,j,l,i]=(WM.POP[k-1,j-1,l,i]+WM.Migration_Scale[k,j,i]*\
                                    WM.C_Share_Migrants[j,l,i]*WM.Migrants[k,i]*(1.+WM.Exog_Npop)**(j-WM.TransYear))*\
                                    (1.-WM.Mortality_rates[k,j,i])

                        WM.POP[k,j,WM.Yclasses,i]+=WM.POP[k,j,l,i]
                else:
                    for l in xrange(WM.Yclasses):
                        for m in xrange(WM.startfertilityage,WM.lastfertilityage):
                            if j<=WM.TransYear:
                                WM.POP[0,j,l,i]+=WM.POP[m,j,l,i]*WM.Fertility[m,j,i]*(1.-WM.Mortality_rates[0,j,i])
                            else:
                                WM.POP[0,j,l,i]=(1.+WM.Exog_Npop)*WM.POP[0,j-1,i,l]

                        WM.POP[0,j,WM.Yclasses,i]+=WM.POP[0,j,l,i]

                WM.POP[91,j,WM.Yclasses,i]+=WM.POP[k,j,WM.Yclasses,i]

    #Print_Population_Summary()
    
    #Calculate the population growth rates
    for i in xrange(WM.countries):
        WM.NPOP[-1,i]=WM.Exog_Npop
        WM.NPOP[0,i]=WM.POP[22,0,WM.Yclasses,i]/WM.POP[WM.startworkingage,0,WM.Yclasses,i]-1.
        for j in xrange(1,WM.Years+1):
            WM.NPOP[j,i]=WM.POP[WM.startworkingage,j,WM.Yclasses,i]/WM.POP[WM.startworkingage,j-1,WM.Yclasses,i]-1.

    #Calculate fertility rates endogenously from year WM.TransYear+1 onwards
    for i in xrange(WM.countries):
        for j in xrange(WM.TransYear+1,WM.Years+1):
            for k in xrange(WM.startfertilityage-WM.lastfertilityage):
                WM.Fertility[k,j,i]=(1.+WM.Exog_Npop)*WM.Fertility[k,j-1,i]*WM.POP[k,j-1,WM.Yclasses,i]/\
                        WM.POP[k,j,WM.Yclasses,Country]
                WM.Total_Fertility[j,i]+=WM.Fertility[k,j,i]


    #Now, we need to shift population acounts accoring to the base year in transition. Russia is just a little different due to
    #their documentation.

    if WM.First_Year != 2000:
        for i in xrange(WM.firstcountry,5):
            for j in xrange(0,WM.Years-WM.First_Year+2000):
                WM.POP[91,j,WM.Yclasses,i]=WM.POP[91,j+WM.First_Year-2000,WM.Yclasses,i]
                for k in xrange(WM.gens+1):
                    WM.POP[k,j,WM.Yclasses,i]=WM.POP[k,j+WM.First_Year-2000,WM.Yclasses,i]
                    WM.Mortality_rates[k,j,i]=WM.Mortality_rates[k,j+WM.First_Year-2000,i]
                    WM.Survival_Probability[k,j,i]=WM.Survival_Probability[k,j+WM.First_Year-2000,i]
                    for l in xrange(WM.Yclasses):
                        WM.POP[k,j,l,i]=WM.POP[k,j+WM.First_Year-2000,l,i]
            
            for j in xrange(WM.gens+j-WM.First_Year+2000):
                WM.Total_Fertility[j,i]=WM.Total_Fertility[j+WM.First_Year-2000,i]
                for k in xrange(WM.startfertilityage-WM.lastfertilityage):
                    WM.Fertility[k,j,i]=WM.Fertility[k,j+WM.First_Year-2000,i]
            #Now We'll fill the remaining years with the original year-300 values

            for j in xrange(WM.Years-WM.First_Year+2001,WM.Years):
                WM.POP[91,j,WM.Yclasses,i]=WM.POP[91,j,WM.Yclasses,i]*(1.+WM.Exog_Npop)**(j-(WM.Years-WM.First_Year+2000))
                for k in xrange(WM.gens+1):
                    WM.POP[k,j,WM.Yclasses,i]=WM.POP[k,j,WM.Yclasses,i]*(1.+WM.Exog_Npop)**(j-(WM.Years-WM.First_Year+2000))
                    WM.Mortality_rates[k,j,i]=WM.Mortality_rates[k,j,i]
                    WM.Survival_Probability[k,j,i]=WM.Survival_Probability[k,j,i]
                    for l in xrange(WM.Yclasses):
                        WM.POP[k,j,l,i]=WM.POP[k,j,l,i]*(1.+WM.Exog_Npop)**(j-(WM.Years-WM.First_Year+2000))

                WM.Total_Fertility[j,i]=WM.Total_Fertility[WM.Years,i]

                for k in xrange(WM.startfertilityage-WM.lastfertilityage):
                    WM.Fertility[k,j,i]=WM.Fertility[k,WM.Years,i]

    #Print_Population_Summary()

    #Adopt population in -1 from Years
    if WM.IRUN==0:
        for i in xrange(WM.countries):
            WM.POP[91,-1,WM.Yclasses,i]=WM.POP[91,WM.Years,WM.Yclasses,i]
            for k in xrange(WM.gens+1):
                WM.POP[k,-1,WM.Yclasses,i]=WM.POP[k,WM.Years,WM.Yclasses,i]
                WM.Mortality_rates[k,-1,i]=WM.Mortality_rates[k,WM.Years,i]
                WM.Survival_Probability[k,-1,i]=WM.Survival_Probability[k,WM.Years,i]
                for l in xrange(WM.Yclasses):
                    WM.POP[k,-1,l,i]=WM.POP[k,WM.Years,l,i]
            
            WM.TotalFertilityss[i]=WM.Total_Fertility[WM.Years,i]

            for k in xrange(WM.startfertilityage-WM.lastfertilityage):
                WM.Fertility_ss[k,i]=WM.Fertility[k,WM.Years,i]
    #Calculate growth rates, again
    for i in xrange(WM.countries):
        WM.NPOP[-1,i]=WM.Exog_Npop
        WM.NPOP[0,i]=WM.POP[WM.startworkingage,0,WM.Yclasses,i]/WM.POP[22,0,WM.Yclasses,i]-1.
        for j in xrange(1,WM.Years+1):
            WM.NPOP[j,i]=WM.POP[WM.startworkingage,j,WM.Yclasses,i]/WM.POP[WM.startworkingage,j-1,WM.Yclasses,i]-1.

    #Efficient Population in each year and country:
    for i in xrange(WM.countries):
        for j in xrange(0,WM.Years+1):
            for k in xrange(WM.Yclasses):
                for l in xrange(WM.startworkingage, WM.gens+1):
                    WM.POP_Efficient[j,i]+=WM.POP[l,j,k,i]*(1.+WM.Tech)**(WM.startworkingage-l)


def Initialize_Variables():
    WM.RG[:]=.1
    WM.AP[:]=1.

    WM.R[:,:]=.1
    WM.Agg_Assets[:,:]=1.

    WM.YY[:,:]=1.
    WM.YY_0[:,:]=1.
    WM.DD[:,:]=1.
    WM.Capital_Tax_Rate[:,:]=1.

    WM.Capital[:,:]=100.
    WM.Consump_Price[:,:]=1.

    WM.Labor[:,:,:]=10.

    WM.Shadow_Wage[0:69,:,:,:]=.0000001
    WM.Avg_Wage_Tax_Rate[0:69,:,:,:]=1.
    WM.Marg_Wage_Tax_Rate[0:69,:,:,:]=1.

    WM.Transfer[21:90,:,:,:]=.05
    WM.H_Transfer[21:90,:,:,:,:]=.05 
