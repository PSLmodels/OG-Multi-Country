import numpy as np
import scipy as sp

##############################################################################################################
#DATA FILES ZONE
#For more information about the individual data files, please refer to DATA_README.txt
input1="Data_Files/china_fertility.csv"
input2="Data_Files/china_mortality.csv"
input3="Data_Files/contributionceiling.csv"
input4="Data_Files/corporatetax.csv"
input5="Data_Files/debtlevel.csv"
input6="Data_Files/eu_fertility.csv"
input7="Data_Files/eu_mortality.csv"
input8="Data_Files/govpy.csv"
input9="Data_Files/healthcare.csv"
input10="Data_Files/india_fertility.csv"
input11="Data_Files/india_mortality.csv"
input12="Data_Files/inheritancetax.csv"
input13="Data_Files/japan_fertility.csv"
input14="Data_Files/japan_mortality.csv"
input15="Data_Files/korea_fertility.csv"
input16="Data_Files/korea_mortality.csv"
input17="Data_Files/mu2gov.csv"
input18="Data_Files/mu2tax.csv"
input19="Data_Files/mu3.csv"
input20="Data_Files/mu4.csv"
input21="Data_Files/net_migration.csv"
input22="Data_Files/popscale.csv"
input23="Data_Files/population.csv"
input24="Data_Files/retirementage.csv"
input25="Data_Files/russia_fertility.csv"
input26="Data_Files/russia_mortality.csv"
input27="Data_Files/skillclasses.csv"
input28="Data_Files/usa_fertility.csv"
input29="Data_Files/usa_mortality.csv"
input30="Data_Files/healthexpenditurescale.csv"
input31="Data_Files/endowment.csv"

#print "Beginning World Module"

########## GLOBAL CONSTANTS ############ This section corresponds to The Global Module part of the Fortran Code


#INTEGER CONSTANTS#####
Beta=np.array([.4,.25]) #This array is the income-class specific labor input in the production function
countries=7 #Indicates the number of countries, can be changed
firstcountry=0
IRUN=0
IRun=0
lastcountry=6
gens=90 #Total number of age groups
lasteducation=20 #Last year that people are educated
lastfertilityage=45 #People stop having children at this age
MaxIter=500 #Caps the number of iterations that model can use
startfertilityage=23 #People start having children at this age
startworkingage=21 #People become independent from their parents at this age
StartYear=2008 #The starting year for our data
Yclasses=2 #Number of skill classes
Years=300 #Total number of years that we will give the model to converge
tol=.0000001 #This is the tolerance for considering the markets cleared, adjust this is the level, corresponds to SigFig
TransYear=50 #The number of years until the model reaches a steady-state
WorkToDeath=gens-startworkingage+1
First_Year=2008

#REAL CONSTANTS###
Alp=1.5 #Leisure Preference Parameter
Alpha=.35 #Share of production by Kapitol
Damp=.7 #These three govern how much to weigh the old and new guesses when iterating.
Dampk=.05
Dampr=.2
Del=.075 #Depreciation Rate
Exog_Npop=0.0 #After the transition year demographic parameters, especially fertility, adjust to set population growth at this level.
Gamma=.25 #Intertemporal Elasticity of Substitution
Hours=1.0 #Time Endowment
Rho=.4 #Intratemporal Elasticity of Subsitution
Tech=.01 #The rate of technological improvement
Theta=1.0 #Weighs the the utility of children
epsilon=.75
Iter=0



#These will change depending on what the functions want, so I just left them as empty
FirstSolutionYear=None
KIter=None
Niter=None
Debt_Level_Refixed=0
Mod_Fert=None


########## NUMPY ARRAYS ############# 
#These are arranged by dimension and alphabetically.

######ONE DIMENSIONAL##############
Years=300
def GetYearBecomingJ(Year,Gen,JJ): 
    '''
    Description:
        This gives the year in which someone age Gen will turn age JJ. This returns year -1 or Years if a year is beyond
        this is called for
    Inputs:
        -Year:The Current Year
        -Gen:Current age of the cohort
        -JJ: The age we want them to turn
    Outputs:
        -GetYearJ: The year in which they will turn JJ
        
    #Fortran 3704-3716, the year in which someone age GEN will turn
    #age JJ. As in the other get year routines, this returns year -1 or YEARS if a year beyond this is
    #is called for
    '''
    GetYearJ=Year+JJ-Gen
    if GetYearJ > Years:
        GetYearJ=Years
    if GetYearJ < 0:
        GetYearJ=0
    if Year==0:
        GetYearJ=0
    if Year==Years:
        GetYearJ=Years

    return GetYearJ



'''
Array Name:Productivity C
Purpose: This array gives how productivity for works of every class and country vary of time.
Fortran:224-314
'''
Productivity_C=np.zeros((Years+2,Yclasses,countries))
for i in xrange(Years+2):
    for j in xrange(3):
        Productivity_C[i,0,j]=1.
        Productivity_C[i,1,j]=1.
    Productivity_C[i,0,5]=1.
    Productivity_C[i,1,5]=1.


Productivity_C[0,0,1]=.6
Productivity_C[0,1,1]=.6
Productivity_C[0,0,2]=.47
Productivity_C[0,1,2]=.47
Productivity_C[0,0,3]=.06
Productivity_C[0,1,3]=.06
Productivity_C[0,0,4]=.035
Productivity_C[0,1,4]=.035
Productivity_C[0,0,5]=.36
Productivity_C[0,1,5]=.36
Productivity_C[0,0,6]=.6
Productivity_C[0,1,6]=.6

Productivity_C[1,0,1]=.6
Productivity_C[1,1,1]=.6
Productivity_C[1,0,2]=.47
Productivity_C[1,1,2]=.47
Productivity_C[1,0,3]=.06
Productivity_C[1,1,3]=.06
Productivity_C[1,0,4]=.035
Productivity_C[1,1,4]=.035
Productivity_C[1,0,5]=.36
Productivity_C[1,1,5]=.36
Productivity_C[1,0,6]=.6
Productivity_C[1,1,6]=.6


Productivity_C[11,0,1]=1.
Productivity_C[11,1,1]=1.
Productivity_C[11,0,2]=1.
Productivity_C[11,1,2]=1.
Productivity_C[41,0,5]=1.
Productivity_C[41,1,5]=1.
Productivity_C[16,0,3]=1.
Productivity_C[16,1,3]=1.
Productivity_C[76,0,4]=1.
Productivity_C[76,1,4]=1.
Productivity_C[11,0,6]=1.
Productivity_C[11,1,6]=1.


for i in xrange(2,Years+2):
    if (i<=10) :
        Productivity_C[i,:,1]=Productivity_C[i-2,:,1]+(Productivity_C[11,:,1]-Productivity_C[1,:,1])/11.

    if (i >= 12) :
        Productivity_C[i,:,1]=Productivity_C[11,:,1]

    if (i <= 10):
        Productivity_C[i,:,2]=Productivity_C[i-2,:,2]+(Productivity_C[11,:,2]-Productivity_C[1,:,2])/12.

    if (i >= 12) :
        Productivity_C[i,:,2]=Productivity_C[11,:,2]

    if (i <= 15) :
        Productivity_C[i,:,3]=Productivity_C[i-2,:,3]+(Productivity_C[16,:,3]-Productivity_C[1,:,3])/16.

    if (i >= 17) :
        Productivity_C[i,:,3]=Productivity_C[16,:,3]

    if (i <= 75) :
        Productivity_C[i,:,4]=Productivity_C[i-2,:,4]+(Productivity_C[76,:,4]-Productivity_C[1,:,4])/76.

    if (i >= 77) :
        Productivity_C[i,:,4]=Productivity_C[76,:,4]

    if (i <= 40) :
        Productivity_C[i,:,5]=Productivity_C[i-2,:,5]+(Productivity_C[41,:,5]-Productivity_C[1,:,5])/41.

    if (i >= 42) :
        Productivity_C[i,:,5]=Productivity_C[41,:,5]

    if (i <= 10) :
        Productivity_C[i,:,6]=Productivity_C[i-2,:,6]+(Productivity_C[11,:,6]-Productivity_C[1,:,6])/12.

    if (i >= 12) :
        Productivity_C[i,:,6]=Productivity_C[11,:,6]


# Number of countries X 1
'''
Array Name:Governments
Purpose:This is the ecogenous government spending, it stands for other expenditures. For example, military spending. This grows with GDP and population.
Fortran:682-690
'''
Govs= np.zeros(countries) #Note that the indeces are one different from the fortran indeces
Govs[0]=.034*Productivity_C[1,0,0]
Govs[1]=.075*Productivity_C[1,0,1]
Govs[2]=.075*Productivity_C[1,0,2]
Govs[3]=.34*Productivity_C[1,0,3]
Govs[4]=.2*Productivity_C[1,0,4]
Govs[5]=.228*Productivity_C[1,0,5]
Govs[6]=.075*Productivity_C[1,0,6]

'''
Array Name:Education Expenditures Profile
Purpose: Lists how much is relatively spent on students of different ages (<21)
Fortran:190-199
'''
Education_Expenditure_Profile=np.array([1.,1.,1.,1.,1.,1.,24.351512,47.25656,49.479592,53.260933,\
            53.593294,54.594752,53.816327,52.482507,51.14723,50.925656,\
            50.814869,49.537901,45.600583,33.516035,20.985432])

'''
Array Name:Population Scale
Purpose:This allows for weighting the populations. This will be multiplied by the population read in from demographics.
Fortran:190-199
'''
Pop_Scale=np.loadtxt((input22),delimiter=',')


'''
Array Name:Total Fertility
Purpose:Used to calculate the fertility rate endogenously
Fortran:1201-1210
'''
TotalFertilityss=np.zeros(countries)

'''
Array Name:Health Expenditure Scale
Purpose: This array scales the health expenditure profiles for each country
Fortran:730-739
'''
Health_Expenditure_Scale=np.zeros(countries)
Health_Expenditure_Scale[0]=.032*Productivity_C[1,0,0]
Health_Expenditure_Scale[1]=.145*Productivity_C[1,0,1]
Health_Expenditure_Scale[2]=.027*Productivity_C[1,0,2]
Health_Expenditure_Scale[3]=.017*Productivity_C[1,0,3]
Health_Expenditure_Scale[4]=.01*Productivity_C[1,0,4]
Health_Expenditure_Scale[5]=.059*Productivity_C[1,0,5]
Health_Expenditure_Scale[6]=.027*Productivity_C[1,0,6]


'''
Array Name:Education Expeditures Scale
Purpose:Scales how much is spend on education
Fortran:776-784
'''
Education_Expenditure_Scale=np.zeros(countries)
Education_Expenditure_Scale[0]=.0033*Productivity_C[1,0,0]
Education_Expenditure_Scale[1]=.0047*Productivity_C[1,0,1]
Education_Expenditure_Scale[2]=.0038*Productivity_C[1,0,2]
Education_Expenditure_Scale[3]=.0042*Productivity_C[1,0,3]
Education_Expenditure_Scale[4]=.0023*Productivity_C[1,0,4]
Education_Expenditure_Scale[5]=.0064*Productivity_C[1,0,5]
Education_Expenditure_Scale[6]=.0064*Productivity_C[1,0,6]


'''
Array Name:GBA (Not really sure what it stands for, it isn't explained)
Purpose: Used in calculating the initial asset distribution
Fortran:142,816-822,878-890,948-958
'''
GBA=np.zeros(14)

'''
Array Name:Initial Asset Share
Purpose:This represents the initial distribution of world assets
Fortran:1540-1547
'''
Initial_Asset_Share=np.array([.319,.217,.148,.097,.035,.009,.011])

### Years+2 X 1 ###
'''
Array Name: Aggregate Assets World
Purpose:Aggregates the world's total assets
Fortran:
'''
Agg_Assets_World=np.ones(Years+2) #Years+2 is used because the orignal fortran code was indexed from -1 to Years and in fortran the upper bound is included. So these will be indexed from 0 to years+2.

'''
Array Name:AP
Purpose:A total factor productivity term
Fortran:1335-1336
'''
AP=np.zeros(Years+2)

'''
Array Name:DD_World
Purpose:Tracks the demand when we're trying to get the markets to clear
Fortran:1407, 2989-3040
'''
DD_World=np.zeros(Years+2)

'''
Array Name:Endowment
Purpose: I think it contains the oil endowment, but I'm not entirely sure what it does yet
Fortran: 3502-3518
'''
Endowment = np.zeros((Years+2))
Endowment= np.loadtxt((input31),delimiter=',')*5.2

'''
Array Name:Foreign Assets World
Purpose:Tracks foriegn assets and trade balance
Fortran:2044-2069
'''
Foreign_Assets_World=np.zeros(Years+2)

'''
Array Name:
Purpose:
Fortran:
'''
Net_Endow_Invest=np.zeros(Years+2)

'''
Array Name:World Interest Rate
Purpose:Tracks the overall interest rate across the globe
Fortran:1340-1341,
'''
RG=np.zeros(Years+2) 

'''
Array Name:
Purpose:
Fortran:
'''
Trade_Balance_World=np.zeros(Years+2)

'''
Array Name:
Purpose:
Fortran:
'''
Tot_Gov_Endow_Share=np.zeros(Years+2)

'''
Array Name:
Purpose:
Fortran:
'''
Tot_Gov_Endow_Revenue=np.zeros(Years+2)

'''
Array Name:
Purpose:
Fortran:
'''
YY_World=np.zeros(Years+2)

# Years+3 X 1
'''
Array Name:
Purpose:
Fortran:
'''
PVEndowment=np.zeros(Years+3) #Years+3 is used because the Fortran code indexes from -1 to years +1

#####################TWO DIMENSIONAL#####################

# Yclasses X Number of Countries
'''
Array Name: Consumption share 
Purpose:Give the population of high and low earners in every country in year 0. Low earners are type 0 and high is type 1
Fortran:201-216
'''
C_Share=np.loadtxt((input27),delimiter=',')

#Years+2 X Number of Countries
'''
Array Name:
Purpose:
Fortran:
'''
Agg_Assets=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Disability_Tax_Rate=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Health_Tax_Rate=np.zeros((Years+2,countries))

'''
Array Name:Average Labor Earnings
Purpose: Calcualtes the average labor earnings of working generations in period YEAR in country LA.
Fortran:2950-2973
'''
Avg_Labor_Earnings=np.zeros((Years+2,countries))


'''
Array Name:
Purpose:
Fortran:
'''
Agg_Pension_Tax_Rate=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Capital=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Capital_Tax_Rate=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
CC=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Consump_Price=np.zeros((Years+2,countries))

'''
Array Name:Contribution Ceiling
Purpose: Each of the countries have a contribution ceiling for Social Security, this is determined here
Fortran:720-727
'''
Contribution_Ceiling=np.loadtxt((input3),delimiter=',')


'''
Array Name:Corporate Taxes
Purpose:Sets the corporate tax rates by year and country
Fortran:613-620
'''
Corp_Tax=np.loadtxt((input4),delimiter=',')

'''
Array Name:
Purpose:
Fortran:
'''
DD=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Debt=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Debt_Endogenous=np.zeros((Years+2,countries))

'''
Array Name: Debt Level
Purpose: This is where the Debt level is set. At least one of the above revenues has to be endogenous in order for countries to say at the target debt level in every year (Needed for convergence)
Fortran:672-680
'''
Debt_Level=np.loadtxt((input5),delimiter=',')

'''
Array Name:
Purpose:
Fortran:
'''


Defecit=np.zeros((Years+2,countries))
'''
Array Name:
Purpose:
Fortran:
'''

Disability_Benefits=np.zeros((Years+2,countries))

'''
Array Name:Disability Benefits Index
Purpose:Tracks the portion of benefits paid for disability insurance
Fortran:763-774
'''
Disability_Benefits_Ind=np.zeros((Years+2,countries))
Disability_Benefits_Ind[:,0]=.015*Productivity_C[1,0,0]
Disability_Benefits_Ind[:,1]=.031*Productivity_C[1,0,1]
Disability_Benefits_Ind[:,2]=.011*Productivity_C[1,0,2]
Disability_Benefits_Ind[:,3]=0.
Disability_Benefits_Ind[:,4]=0.
Disability_Benefits_Ind[:,5]=0.
Disability_Benefits_Ind[:,6]=0.

'''
Array Name:
Purpose:
Fortran:
'''
Education_Expenditures=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Foreign_Assets=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Government_Discretionary_Spending=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Government_Expenditures=np.zeros((Years+2,countries))

'''
Array Name:Government Endowment Share
Purpose:Keeps track of the share of assets that belong to the government.
Fortran:486-496
'''
Gov_Endow_Share=np.zeros((Years+2,countries))
Gov_Endow_Share[:,0]=.0195
Gov_Endow_Share[:,1]=.0315
Gov_Endow_Share[:,2]=0.
Gov_Endow_Share[:,3]=.0335
Gov_Endow_Share[:,4]=.0125
Gov_Endow_Share[:,5]=.29
Gov_Endow_Share[:,6]=0.

'''
Array Name:
Purpose:
Fortran:
'''
Health_Benefits=np.zeros((Years+2,countries))

'''
Array Name:Inheritance Tax Rate
Purpose:Keeps track of the tax rate paid on bequests
Fortran:604-611
'''
Inheritance_Tax_Rate=np.loadtxt((input12),delimiter=',') #Could have just left it at all zeros, but it could be useful to model inheritance taxes later, so I left it in.

'''
Array Name:
Purpose:
Fortran:
'''

Invest=np.zeros((Years+2,countries))

'''
Array Name: GOVPY 
Purpose: This is the part of the pension beneifts contributed by the government (the other part is financed by taxes)
Fortran:623-630
'''
Mu1=np.loadtxt((input8),delimiter=',')

'''
Array Name: Mu2 Government
Purpose:The part of the Health benefits treated as general government consumption
Fortran:632-639
'''
Mu2_gov=np.loadtxt((input17),delimiter=',')

'''
Array Name:Mu2 Tax
Purpose:This is the part of Health benefits financed by general taxes
Fortran:641-647
'''
Mu2_tax=np.loadtxt((input18),delimiter=',')

'''
Array Name: Mu3 (Disability)
Purpose:The part of disability system paid for through general taxes
Fortran:650-657
'''
Mu3=np.loadtxt((input19),delimiter=',')

'''
Array Name:Mu 4
Purpose:The share of corporate tax revenues which is lump-sum transferred to households
Fortran:659-666
'''
Mu4=np.loadtxt((input20),delimiter=',')


'''
Array Name: Population Growth
Purpose: I think it stores population growth rates, but I'm not sure
Fortran: 1192-1199, Assigned values in Population_Development function
'''
NPOP=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Pension_Benefits=np.zeros((Years+2,countries))

'''
Array Name: Efficient Population
Purpose: Not sure yet, but it seems to be a rate that measures productivity increases
Fortran: Assigned values in 1291-1301 in Population_Development function 
'''
POP_Efficient=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
R=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
RetirementAge=np.loadtxt((input24),delimiter=',').astype(int)

'''
Array Name:
Purpose:
Fortran:
'''
Total_Expenditures=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Trade_Balance=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
TRF=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
YY=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
YY_0=np.zeros((Years+2,countries))

#Last year of Education +1 X Countries
'''
Array Name:Education Expenditures Index
Purpose: Multiplies by the expenditure profile, previously, to ge the spending in each year on each age
Fortran:787-792
'''
Education_Expenditures_Ind=np.zeros((lasteducation+1,countries))
for i in xrange(countries):
    for j in xrange(lasteducation):
        Education_Expenditures_Ind[j,i]=Education_Expenditure_Profile[j]*Education_Expenditure_Scale[i]


#gens+1 X countries
'''
Array Name: Health_Expenditures_Profile
Purpose: Data from csv 'healthcare'. Not quite sure what it's purpose is yet. Contains all 0s for China and India
Fortran: 3439-3445
'''
Health_Expenditures_Profile=np.zeros((gens+1,countries))
Health_Expenditures_Profile[:,:] = np.loadtxt((input9),delimiter=',')
#Sets Korea's equal to Japan's
Health_Expenditures_Profile[:,6] = Health_Expenditures_Profile[:,2]

'''
Array Name:
Purpose:
Fortran:
'''
Migrants = np.zeros((gens+1, countries))
Migrants[1:66,:] = np.loadtxt((input21), delimiter=',')*100


#Fertility Years X Countries
'''
Array Name:
Purpose:
Fortran:
'''
Fertility_ss=np.zeros((lastfertilityage+1-startfertilityage,countries))

#Years +2 X Countries
'''
Array Name:Country Endowment Share
Purpose:These must sum to 1 each year, it is used to calculate GDP figures
Fortran:499-506
'''
Country_Endow_Share=np.zeros((Years+2,countries))
Country_Endow_Share[:,0]=.1938
Country_Endow_Share[:,1]=.0793
Country_Endow_Share[:,2]=0.
Country_Endow_Share[:,3]=.0838
Country_Endow_Share[:,4]=.031
Country_Endow_Share[:,5]=.612
Country_Endow_Share[:,6]=0.

'''
Array Name:
Purpose:
Fortran:
'''
GDP=np.zeros((Years+2,countries))

'''
Array Name:
Purpose:
Fortran:
'''
National_Income=np.zeros((Years+2,countries))

#THREE DIMENSIONAL######

#Years+2 X YClasses X Countries
'''
Array Name:
Purpose:
Fortran:
'''
Agg_Assets_For_Bequests=np.zeros((Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Avg_Indexed_Earnings=np.zeros((Years+2,Yclasses,countries))

'''
Array Name:Consumption Share of Migrants
Purpose:Tracks the share of consumption from incoming migrants
Fortran:334-340
'''
C_Share_Migrants=np.zeros((Years+2,Yclasses,countries))
C_Share_Migrants = np.tile(C_Share.reshape(1, Yclasses, countries), (Years+2, 1, 1))

'''
Array Name:
Purpose:
Fortran:
'''
Wage_Index_Class=np.ones((Years+2,Yclasses,countries))

#Years+2 X YClasses+1 X Countries
'''
Array Name:
Purpose:
Fortran:
'''
Agg_Assets_Migrants=np.zeros((Years+2,Yclasses+1,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Avg_Disability_Tax_Rate=np.zeros((Years+2,Yclasses+1,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Avg_Health_Tax_Rate=np.zeros((Years+2,Yclasses+1,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Avg_Pension_Tax_Rate=np.zeros((Years+2,Yclasses+1,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Avg_Wage_Tax_Rate=np.zeros((Years+2,Yclasses+1,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Marg_Disability_Tax_Rate=np.zeros((Years+2,Yclasses+1,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Marg_Health_Tax_Rate=np.zeros((Years+2,Yclasses+1,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Marg_Pension_Tax_Rate=np.zeros((Years+2,Yclasses+1,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Agg_Marg_Wage_Tax_Rate=np.zeros((Years+2,Yclasses+1,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Labor=np.zeros((Years+2,Yclasses+1,countries))

#gens+1 X Years+2 X Countries
'''
Array Name:Health Benefits Ind
Purpose:Calculates the government health benefits for people of every age, use the scale parameters. For where there is an expenditure missing, another country's distribution is used.
Fortran:
'''
Health_Benefits_Ind=np.zeros((gens+1,Years+2,countries))
for i in xrange(countries):
    for j in xrange(Years+2):
        if i <= 2 or i==6:
            Health_Benefits_Ind[:,j,i]=Health_Expenditures_Profile[:,i]*Health_Expenditure_Scale[i]
        elif i==5:
            Health_Benefits_Ind[:,j,i]=Health_Expenditures_Profile[:,i]*Health_Expenditure_Scale[i]
        elif i==3:
            Health_Benefits_Ind[:,j,i]=Health_Expenditures_Profile[:,2]*Health_Expenditure_Scale[i]
        else:
            Health_Benefits_Ind[:,j,i]=Health_Expenditures_Profile[:,2]*Health_Expenditure_Scale[i]

'''
Array Name:Survival Probability
Purpose: Contains the probability of generation g surviving until year t
Fortran: 1092-1103
'''
Survival_Probability=np.ones((gens+1,Years+2,countries))

'''
Array Name: Mortality rates
Purpose: Constains the probability of dying
Fortran: 3348-3437
'''
Mortality_rates = np.zeros((gens+1, Years+2, countries))
Mortality_rates[68:91,:51,:] = np.dstack((np.transpose(np.loadtxt((input29),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input7),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input14),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input2),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input11),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input26),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input16),delimiter=',')) \
                                ))
Mortality_rates[68:91,51:,:] = np.tile(Mortality_rates[68:91,50,:].reshape(23, 1, 7), (1, 251, 1))

'''
Array Name:
Purpose:
Fortran:
'''

#gens+1 X Years+2 X Countries
'''
Array Name:Migration Scale
Purpose:A multiplier for the amount of migration read in from the demographic file.
Fortran:
'''
Migration_Scale=np.ones((gens+1,Years+2,countries))
Migration_Scale[:,:,6] = Migration_Scale[:,:,6]*.5


#3 X Years+2 X Countries
'''
Array Name:Aleph
Purpose:This array records the proportional term
Fortran:2353
'''
Aleph=np.zeros((3,Years+2,countries))

'''
Array Name:Beth
Purpose:This is the progressive Term
Fortran:2353
'''
Beth=np.zeros((3,Years+2,countries))

#6 X Years+2 X Countries
'''
Array Name:
Purpose:
Fortran:
'''
Tax_Revenues=np.zeros((6,Years+2,countries))

# Fertility Years X gens+years+1 X countries
'''
Array Name:Fertility
Purpose:
Fortran: Readin:3250-3340, Assigning Values: 3453-3469
'''
Fertility=np.zeros((lastfertilityage+1 - startfertilityage,gens+Years+1,countries))

Fertdatacsv = np.zeros((lastfertilityage+1 - startfertilityage, 98, countries))#This is just a dummy variable to assist in reading in and storing Fertility data in python
Fertdatacsv = np.dstack((np.transpose(np.loadtxt((input28),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input6),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input13),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input1),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input10),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input25),delimiter=',')) \
                                ,np.transpose(np.loadtxt((input15),delimiter=',')) \
                                ))
Fertility[:,0:51,:] = Fertdatacsv[:,48:,:] #Fertility values from 2008-2058. Stored in the first 51 entries of Fertility
Fertility[:,343:,:] = Fertdatacsv[:,:48,:] #Fertility values from 1960-2007. All the negative years in the csv files. Stored in Fertility 343-390 so to directly relate to the Fortran code that called these values from negative numbers 
Fertility[:,301:343,:] = np.tile(Fertdatacsv[:,0,:].reshape(23, 1, 7), (1, 42, 1)) #All Fertility from 1918-1959 is equal to the Fertility rate in 1960
Fertility[:,51:301,:] = np.tile(Fertdatacsv[:,-1,:].reshape(23, 1, 7), (1, 250, 1)) #All Fertility from 2051-2038 is equal to the Fertility rate in 2051. This changes in Population_Development function
Fertility[:,:,6] = Fertility[:,:,6]*.9 #Scales Russia Fertility
#Fertility[:,:,7] = Fertility[:,:,7]*.1.5

'''
Array Name:
Purpose:
Fortran:
'''
# gens+Years+1 X Number of Countries
'''
Array Name:Total_Fertility
Purpose:Stores the total fertility rate for a given year and region (sum of all generations)
Fortran:3472-3482
NOTE:Has same strange indexing as Fertility
'''
Total_Fertility=np.zeros((gens+Years+1,countries))
Total_Fertility=np.sum(Fertility, axis=0)

'''
Array Name:
Purpose: Contains the average age for any given year and region
Fortran:3472-3482
NOTE:Has same strange indexing as Fertility
'''
Avg_Birth_Age = np.zeros((gens+Years+1, countries))
Avg_Birth_Age[:,:] = np.tensordot(Fertility, xrange(startfertilityage, lastfertilityage+1), axes=([0],[0]))/Total_Fertility

#Fertility Years X Years +2 X Countries
'''
Array Name:Delta
Purpose:Keeps track of the individual deltas used in the maximization problem
Fortran:352-393
'''
Delta=np.zeros((WorkToDeath,Years+2,countries))

for j in xrange(WorkToDeath):
    Delta[j,:,0]=.01
    Delta[j,:,1]=-.02
    Delta[j,:,2]=-.02
    Delta[j,:,3]=-.07
    Delta[j,:,4]=-.026
    Delta[j,:,5]=-.07
    Delta[j,:,6]=-.01
        #Check for delta changing over time 366-375
'''
for i in xrange(3,3):
    if i==0:
        Delta[0,0,i]=-.07
        Delta[0,26,i]=.01
    for j in xrange(2,Years+1):
        if i==3:
            if j <=25:
                Delta[0,j,i]=Delta[0,j-1,i]+(Delta[0,26,i]-Delta[0,1,i])/26.*1.
            if j >=27:
                Delta[0,j,i]=Delta[0,26,i]*1.
    for k in xrange(22,gens):
        for l in xrange(gens):
            IK=GetYearBecomingJ(0,k,l)
            Delta[l,IK,i]=Delta[l,0,i]*1.
    for k in xrange(Years+1):
        for l in xrange(WorkToDeath):
            IK=GetYearBecomingJ(k,0,l)
            Delta[l,IK,i]=Delta[0,k,i]*1.
'''

#gens +1 X Years+2 X Countries
'''
Array Name:
Purpose:
Fortran:
'''
Pension_Replacement_Rate=np.zeros((gens+1,Years+2,Yclasses,countries))

#3 X Years + 2 X Countries
'''
Array Name:Endogenous Tax Ratio
Purpose:This is the sets the actual rate
Fortran:520-620
'''
Endogenous_Tax_Ratio=np.zeros((3,Years+2,countries))


    #LowerSkill Class Endogenous Tax Ratio
Endogenous_Tax_Ratio[0,0:,0]=.5
Endogenous_Tax_Ratio[0,0:,1]=.65
Endogenous_Tax_Ratio[0,0:,2]=.58
Endogenous_Tax_Ratio[0,0:,3]=.83
Endogenous_Tax_Ratio[0,0:,4]=.83
Endogenous_Tax_Ratio[0,0:,5]=0.
Endogenous_Tax_Ratio[0,0:,6]=0.

    #HigherSkill Class Endogenous Tax Ratio
Endogenous_Tax_Ratio[1,0:,0]=.5
Endogenous_Tax_Ratio[1,0:,1]=.35
Endogenous_Tax_Ratio[1,0:,2]=.42
Endogenous_Tax_Ratio[1,0:,3]=.17
Endogenous_Tax_Ratio[1,0:,4]=.17
Endogenous_Tax_Ratio[1,0:,5]=1.
Endogenous_Tax_Ratio[1,0:,6]=1.


#Working Years+1X Yclasses X Countries REMEMBER THIS CHANGE!!!!!!!!!!!!!!!!!
'''
Array Name: Assets Initial Year
Purpose: This suts up the initial asset distribution by class and age. This is just an initial guess, the initial stead state function in ClawsonOlmstead.py will get the actualy stead state level and distribution.
Fortran:797-1001
'''
Assets_Initial_Year=np.zeros((WorkToDeath,Yclasses,countries))
Assets_Initial_Year[0,:,0]=0.
Assets_Initial_Year[7,:,0]=65.9*.01025*Productivity_C[1,:,0]*.826*.4
Assets_Initial_Year[14,:,0]=196.2*0.01025*Productivity_C[1,:,0]*.826*.4
Assets_Initial_Year[26,:,0]=362.7*0.01025*Productivity_C[1,:,0]*.826*.4
Assets_Initial_Year[43,:,0]=500.2*0.01025*Productivity_C[1,:,0]*.826*.4
Assets_Initial_Year[53,:,0]=465.5*0.01025*Productivity_C[1,:,0]*.826*.4
Assets_Initial_Year[61,:,0]=310.2*0.01025*Productivity_C[1,:,0]*.826*.4

for i in xrange(Yclasses):
    GBA[0]=(Assets_Initial_Year[7,i,0]-Assets_Initial_Year[0,i,0])/7.
    GBA[1]=(Assets_Initial_Year[14,i,0]-Assets_Initial_Year[7,i,0])/7.
    GBA[2]=(Assets_Initial_Year[26,i,0]-Assets_Initial_Year[14,i,0])/12.
    GBA[3]=(Assets_Initial_Year[43,i,0]-Assets_Initial_Year[26,i,0])/17.
    GBA[4]=(Assets_Initial_Year[53,i,0]-Assets_Initial_Year[43,i,0])/10.
    GBA[5]=(Assets_Initial_Year[61,i,0]-Assets_Initial_Year[53,i,0])/8.
    GBA[6]=(0.-Assets_Initial_Year[61,i,0])/9.

    for j in xrange(1,8):
        Assets_Initial_Year[j,i,0]=Assets_Initial_Year[j-1,i,0]+GBA[0]
    for j in xrange(8,14):
        Assets_Initial_Year[j,i,0]=Assets_Initial_Year[j-1,i,0]+GBA[1]
    for j in xrange(15,26):
        Assets_Initial_Year[j,i,0]=Assets_Initial_Year[j-1,i,0]+GBA[2]
    for j in xrange(27,43):
        Assets_Initial_Year[j,i,0]=Assets_Initial_Year[j-1,i,0]+GBA[3]
    for j in xrange(44,53):
        Assets_Initial_Year[j,i,0]=Assets_Initial_Year[j-1,i,0]+GBA[4]
    for j in xrange(54,61):
        Assets_Initial_Year[j,i,0]=Assets_Initial_Year[j-1,i,0]+GBA[5]
    for j in xrange(62,70):
        Assets_Initial_Year[j,i,0]=Assets_Initial_Year[j-1,i,0]+GBA[6]

    #European Union

Assets_Initial_Year[0,:,1]=0.
Assets_Initial_Year[1,:,1]=51.2730*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[4,:,1]=65.3330*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[9,:,1]=89.013*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[14,:,1]=136.144*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[19,:,1]=187.2660*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[24,:,1]=233.0530*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[29,:,1]=279.581*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[34,:,1]=318.0230*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[39,:,1]=289.089*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[44,:,1]=251.6210*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[49,:,1]=200.5960*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[54,:,1]=186.5430*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[59,:,1]=150.860*0.0150*Productivity_C[1,:,1]*.84090*.4
Assets_Initial_Year[69,:,1]=50.290*0.0150*Productivity_C[1,:,1]*.84090*.4

for i in xrange(Yclasses):
    GBA[0]=(Assets_Initial_Year[4,i,1]-Assets_Initial_Year[1,i,1])/3.
    GBA[1]=(Assets_Initial_Year[9,i,1]-Assets_Initial_Year[4,i,1])/5.
    GBA[2]=(Assets_Initial_Year[14,i,1]-Assets_Initial_Year[9,i,1])/5.
    GBA[3]=(Assets_Initial_Year[19,i,1]-Assets_Initial_Year[14,i,1])/5.
    GBA[4]=(Assets_Initial_Year[24,i,1]-Assets_Initial_Year[19,i,1])/5.
    GBA[5]=(Assets_Initial_Year[29,i,1]-Assets_Initial_Year[24,i,1])/5.
    GBA[6]=(Assets_Initial_Year[34,i,1]-Assets_Initial_Year[29,i,1])/5.
    GBA[7]=(Assets_Initial_Year[39,i,1]-Assets_Initial_Year[34,i,1])/5.
    GBA[8]=(Assets_Initial_Year[44,i,1]-Assets_Initial_Year[39,i,1])/5.
    GBA[9]=(Assets_Initial_Year[49,i,1]-Assets_Initial_Year[44,i,1])/5.
    GBA[10]=(Assets_Initial_Year[54,i,1]-Assets_Initial_Year[49,i,1])/5.
    GBA[11]=(Assets_Initial_Year[59,i,1]-Assets_Initial_Year[54,i,1])/5.
    GBA[12]=(Assets_Initial_Year[69,i,1]-Assets_Initial_Year[59,i,1])/10.


    for j in xrange(2,4):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[0]
    for j in xrange(5,9):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[1]
    for j in xrange(10,14):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[2]
    for j in xrange(15,19):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[3]
    for j in xrange(20,24):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[4]
    for j in xrange(25,29):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[5]
    for j in xrange(30,34):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[6]
    for j in xrange(35,39):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[7]
    for j in xrange(40,44):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[8]
    for j in xrange(45,49):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[9]
    for j in xrange(50,54):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[10]
    for j in xrange(55,59):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[11]
    for j in xrange(60,69):
        Assets_Initial_Year[j,i,1]=Assets_Initial_Year[j-1,i,1]+GBA[12]

#Japan

Assets_Initial_Year[0,:,2]=0.
Assets_Initial_Year[1,:,2]=29.020*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[5,:,2]=58.030*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[11,:,2]=107.170*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[16,:,2]=156.320*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[21,:,2]=218.280*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[26,:,2]=280.240*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[31,:,2]=337.760*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[36,:,2]=395.270*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[41,:,2]=448.930*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[46,:,2]=502.590*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[49,:,2]=527.010*0.0097*Productivity_C[1,:,2]*.92010*.4
Assets_Initial_Year[69,:,2]=105.40*0.0097*Productivity_C[1,:,2]*.92010*.4

for i in xrange(Yclasses):
    GBA[0]=(Assets_Initial_Year[5,i,2]-Assets_Initial_Year[1,i,2])/4.
    GBA[1]=(Assets_Initial_Year[11,i,2]-Assets_Initial_Year[5,i,2])/6.
    GBA[2]=(Assets_Initial_Year[16,i,2]-Assets_Initial_Year[11,i,2])/5.
    GBA[3]=(Assets_Initial_Year[21,i,2]-Assets_Initial_Year[16,i,2])/5.
    GBA[4]=(Assets_Initial_Year[26,i,2]-Assets_Initial_Year[21,i,2])/5.
    GBA[5]=(Assets_Initial_Year[31,i,2]-Assets_Initial_Year[26,i,2])/5.
    GBA[6]=(Assets_Initial_Year[36,i,2]-Assets_Initial_Year[31,i,2])/5.
    GBA[7]=(Assets_Initial_Year[41,i,2]-Assets_Initial_Year[36,i,2])/5.
    GBA[8]=(Assets_Initial_Year[46,i,2]-Assets_Initial_Year[41,i,2])/5.
    GBA[9]=(Assets_Initial_Year[49,i,2]-Assets_Initial_Year[46,i,2])/3.
    GBA[10]=(Assets_Initial_Year[69,i,2]-Assets_Initial_Year[49,i,2])/20.


    for j in xrange(2,5):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[0]
    for j in xrange(6,11):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[1]
    for j in xrange(12,16):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[2]
    for j in xrange(17,21):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[3]
    for j in xrange(22,26):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[4]
    for j in xrange(27,31):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[5]
    for j in xrange(32,36):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[6]
    for j in xrange(37,41):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[7]
    for j in xrange(42,46):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[8]
    for j in xrange(47,49):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[9]
    for j in xrange(50,69):
        Assets_Initial_Year[j,i,2]=Assets_Initial_Year[j-1,i,2]+GBA[10]

#The Rest: China, India, Russia and Korea

for i in xrange(Yclasses):
    for j in xrange(WorkToDeath):
        Assets_Initial_Year[j,i,3]=Assets_Initial_Year[j,i,2]*.1
        Assets_Initial_Year[j,i,4]=Assets_Initial_Year[j,i,2]*.1
        Assets_Initial_Year[j,i,5]=Assets_Initial_Year[j,i,2]*.1
        Assets_Initial_Year[j,i,6]=Assets_Initial_Year[j,i,2]*.1 

#print Assets_Initial_Year[7,:,0]


'''
Array Name:
Purpose:
Fortran:
'''
Base_Assets=np.zeros((WorkToDeath,Yclasses,countries))

'''
Array Name:ZXI
Purpose:Used to calculate consumption in future years, it's a factor to multiply the utility function by.
Fortran:1850-1852
'''
ZXI=np.zeros((WorkToDeath,Yclasses,countries))


'''
Array Name:OAP
Purpose:This code allows workers of different ages to have different productivies. Those who are past retirement ages have a productivity of 0.
Fortran:407-450
'''
OAP=np.zeros((WorkToDeath,Years+2,countries))
for i in xrange(countries):
    for j in xrange(Years+2):
        for k in xrange(WorkToDeath):
            if i==0 and j<=41:
                OAP[k,j,i]=1.0
            elif i==0 and j >= 42:
                OAP[k,j,i]=0.
            elif i==1 and j <= 38:
                OAP[k,j,i]=1.
            elif i==1 and j >= 39:
                OAP[k,j,i]=0.
            elif i==2 and j <= 38:
                OAP[k,j,i]=1.
            elif i==2 and j >= 39:
                OAP[k,j,i]=0.
            elif i==3 and j <= 38:
                OAP[k,j,i]=1.
            elif i==3 and j >= 39:
                OAP[k,j,i]=0.
            elif i==4 and j <= 38:
                OAP[k,j,i]=1.
            elif i==4 and j >= 39:
                OAP[k,j,i]=0.
            elif i==5 and j <= 36:
                OAP[k,j,i]=1.
            elif i==5 and j >= 38:
                OAP[k,j,i]=0.
            elif i==6 and j <= 43:
                OAP[k,j,i]=1.
            elif i==6 and j >= 44:
                OAP[k,j,i]=0.

#FOUR DIMENSIONAL#######

#gens +1 X Years +2 X Yclasses +1 X Countries
'''
Array Name:POP
Purpose: Stores the population data. The 92 generation represents the sum of all other generations, and the 3 Yclass is the sum of the population in each skill class
Fortran: 1149-1187
'''
POP=np.zeros((gens+2,Years+2,Yclasses+1,countries))
#POP[:-1,0,:-1,:] = np.multiply(np.tile(np.loadtxt((input23), delimiter=',').reshape((gens+1, 1, countries)),(1, Yclasses, 1)),np.tile(C_Share.reshape((1, Yclasses, countries)), (gens+1, 1, 1)))*1000

POP[:-1,0,:-1,:] = np.einsum('gr,kr -> gkr', np.loadtxt((input23), delimiter=','), C_Share)*Pop_Scale*1000



#4 X 2 X Years+1 X Countries
'''
Array Name:Tax Rate
Purpose:This tracks several of the different taxes charged (Type of parameter, proportional/progressive,year,country)
Fortran:520-667
'''
#The first four are the different types of taxes, the next two are proportional/progressive, year and country
Tax_Rate=np.zeros((4,2,Years+2,countries))
   
    #Consumption Tax
Tax_Rate[0,0,:,0]=-1.
Tax_Rate[0,1,:,0]=0.
Tax_Rate[0,0,:,1]=-1.
Tax_Rate[0,1,:,1]=0.
Tax_Rate[0,0,:,2]=-1.
Tax_Rate[0,1,:,2]=0.
Tax_Rate[0,0,:,3]=-1.
Tax_Rate[0,1,:,3]=0.
Tax_Rate[0,0,:,4]=-1.
Tax_Rate[0,1,:,4]=0.
Tax_Rate[0,0,:,5]=.19
Tax_Rate[0,1,:,5]=0.
Tax_Rate[0,0,:,6]=.1
Tax_Rate[0,1,:,6]=0.

    #Wage Tax
Tax_Rate[1,0,:,0]=-1.
Tax_Rate[1,1,:,0]=0.009
Tax_Rate[1,0,:,1]=-1.
Tax_Rate[1,1,:,1]=0.009
Tax_Rate[1,0,:,2]=-1.
Tax_Rate[1,1,:,2]=0.009
Tax_Rate[1,0,:,3]=-1.
Tax_Rate[1,1,:,3]=0.009
Tax_Rate[1,0,:,4]=-1.
Tax_Rate[1,1,:,4]=0.009
Tax_Rate[1,0,:,5]=.19
Tax_Rate[1,1,:,5]=0.009
Tax_Rate[1,0,:,6]=.1
Tax_Rate[1,1,:,6]=0.

    #Capital Tax
Tax_Rate[2,0,:,0]=.11
Tax_Rate[2,1,:,0]=0.
Tax_Rate[2,0,:,1]=.14
Tax_Rate[2,1,:,1]=0.
Tax_Rate[2,0,:,2]=.08
Tax_Rate[2,1,:,2]=0.
Tax_Rate[2,0,:,3]=.03
Tax_Rate[2,1,:,3]=0.
Tax_Rate[2,0,:,4]=.018
Tax_Rate[2,1,:,4]=0.
Tax_Rate[2,0,:,5]=.24
Tax_Rate[2,1,:,5]=0.005
Tax_Rate[2,0,:,6]=.15
Tax_Rate[2,1,:,6]=0.

    #Parameter Replacement Rate    
Tax_Rate[3,0,:,0]=.65
Tax_Rate[3,1,:,0]=0.1
Tax_Rate[3,0,:,1]=.6
Tax_Rate[3,1,:,1]=0.15
Tax_Rate[3,0,:,2]=.4
Tax_Rate[3,1,:,2]=0.15
Tax_Rate[3,0,:,3]=.3
Tax_Rate[3,1,:,3]=0.
Tax_Rate[3,0,:,4]=.3
Tax_Rate[3,1,:,4]=0.
Tax_Rate[3,0,:,5]=.63
Tax_Rate[3,1,:,5]=0.02
Tax_Rate[3,0,:,6]=.5
Tax_Rate[3,1,:,6]=0.1


#Working Years X Years+2 X Yclasses X Countries
'''
Array Name:Age Efficiency
Purpose:This array uses the productivity parameters to create actual efficiences by age. Uses the exponential function to match an empirical wage distribution.
Fortran:452-473
'''
Age_Efficiency=np.zeros((WorkToDeath,Years+2,Yclasses,countries))
for i in xrange(countries):
    for j in xrange(Years+2):
        for k in xrange(Yclasses):
            for l in xrange(WorkToDeath): #WorkToDeath
                if j<=0 or (j>0 and l==0):
                    Age_Efficiency[l,j,k,i]=Productivity_C[j,k,i]*np.exp(4.46+0.033*(l+1)\
                    -.00067*(l+1)**2)/np.exp(4.50233)*(1.+Tech)**(l)*OAP[l,j,i]
                elif j>0 and l>0:
                    IM21=GetYearBecomingJ(j,l,0)
                    Age_Efficiency[l,j,k,i]=Productivity_C[j,k,i]*np.exp(4.46+0.033*(l+1)\
                    -.00067*(l+1)**2)/np.exp(4.50233)*(1.+Tech)**(l)*OAP[l,j,i]


'''
Array Name:
Purpose:
Fortran:
'''
Assets=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Avg_Health_Tax_Rate=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Avg_Pension_Tax_Rate=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Avg_Wage_Tax_Rate=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Avg_Disability_Tax_Rate=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Consump=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Consump_Kids=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Leisure=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Marg_Disability_Tax_Rate=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Marg_Health_Tax_Rate=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Marg_Pension_Tax_Rate=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Marg_Wage_Tax_Rate=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Pension_Benefits_Ind=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Shadow_Wage=np.ones((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Transfer=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

'''
Array Name:
Purpose:
Fortran:
'''
Utility=np.zeros((WorkToDeath,Years+2,Yclasses,countries))

#Working Years X Years+1 X Yclasses X Countries
'''
Array Name:
Purpose:
Fortran:
'''
Utility_0=np.zeros((WorkToDeath,Years+1,Yclasses,countries))

#Gens+1 X Years+2, 46-23, Yclasses+1, Countries
'''
Array Name: Kid Population
Purpose:Tracks the population of children
Fortran: Russia 75
'''
Kid_Pop=np.zeros((gens+2,Years+2,46-23+1,Yclasses+1,countries))

'''
Array Name: Inheritance
Purpose: Keeps track of generational inheritance
Fortran: Russia 76
'''
Inheritance_By_Age_Parent=np.zeros((gens+2,Years+2,46-23+1,Yclasses+1,countries))


#FIVE DIMENSIONAL######

#Working Years X Years+2 X Yclasses X Countries X Countries
'''
Array Name:
Purpose:
Fortran:
'''
H_Transfer=np.zeros((WorkToDeath,Years+2,Yclasses,countries,countries))



#############################################################################################

#FILLING IN DATA ZONE STARTING WITH:
#Productivity_C

#Shr
shr=np.zeros(gens+1)
shr[21]=.137
shr[22]=.139
shr[23]=.256
shr[24]=.451
shr[25]=.757
shr[26]=1.207
shr[27]=1.834
shr[28]=2.652
shr[29]=3.65
shr[30]=4.784
shr[31]=5.969
shr[32]=7.091
shr[33]=8.018
shr[34]=8.632
shr[35]=8.847
shr[36]=8.632
shr[37]=8.018
shr[38]=7.091
shr[39]=5.969
shr[40]=4.784
shr[41]=3.650
shr[42]=2.652
shr[43]=1.834
shr[44]=1.207
shr[45]=.757
shr[46]=.451
shr[47]=.256
shr[48]=.139
shr[49]=.136

#Filling in the Total governmetn endowment share
for i in xrange(Years+2):
    for j in xrange(countries):
        Tot_Gov_Endow_Share[i]+=Gov_Endow_Share[i,j]



'''

shr=np.array([.137,.139,.256,.451,.757,1.207,1.834,2.652,3.65,4.784,\
            5.969,7.091,8.018,8.632,8.847,8.632,8.018,7.091,\
            5.969,4.784,3.650,2.652,1.834,1.207,.757,.451,.256,.139,.136])
'''

