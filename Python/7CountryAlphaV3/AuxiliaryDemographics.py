from __future__ import division
import csv
import numpy as np
import scipy as sp
import scipy.optimize as opt
import time as time
from matplotlib import pyplot as plt


#DEMOGRAPHICS FUNCTIONS
def getkeyages(S, PrintAges):
    """
    Description:
        -Gets the key ages for calculating demographics based on S

    Inputs:
        S                = Int in [10,80], Number of cohorts
        PrintAges        = Boolean, Prints key ages if true

    Functions called:
        -None

    Objects in Function:
        agestopull        = (S) List, Contains the ages from the data files to use for the demographic data in the program
        LeaveHouseAge     = Int, Age where children become adults and are no longer dependent on parents for consumption
        FirstFertilityAge = Int in (0,S), First age when agents bear children
        LastFertilityAge  = Int in (0,S), Last age when agents bear children
        FirstDyingAge     = Int in (0,S), First age when agents die
        MaxImmigrantAge   = Int in (0,S), Age of the oldest immigrants
       
    Returns: LeaveHouseAge, FirstFertilityAge, LastFertilityAge, MaxImmigrantAge, FirstDyingAge, agestopull
    """

    agestopull = np.arange(80)[np.where(np.in1d(np.arange(80), np.round(np.linspace(0, 80, num=S, endpoint=False))))]
    LeaveHouseAge = np.abs(agestopull - 21).argmin()#The age when agents become independent households and don't rely on parents consumption
    FirstFertilityAge = np.abs(agestopull - 23).argmin()#The age when agents have their first children
    LastFertilityAge = np.abs(agestopull - 45).argmin()#The age when agents have their last children
    MaxImmigrantAge = np.abs(agestopull - 65).argmin()#All immigrants are between ages 0 and MaxImmigrantAge
    FirstDyingAge = np.abs(agestopull - 68).argmin()#The first age agents can begin to die

    #Insures that parents don't die with kids in the home
    if FirstDyingAge - (LeaveHouseAge + LastFertilityAge) <= 0:
        FirstDyingAge+=1

    #Makes sure we aren't pulling anything out of bounds from the data files
    if agestopull[FirstDyingAge] < 68:
        agestopull[FirstDyingAge] = 45
    if agestopull[FirstFertilityAge] < 23:
        agestopull[FirstFertilityAge] = 23
    if agestopull[LastFertilityAge] > 45:
        agestopull[LastFertilityAge] = 45

    if PrintAges:
        print "\t\t\t\t\tModel\tData"
        print "LeaveHouseAge:\t\t", LeaveHouseAge, "\t\t", agestopull[LeaveHouseAge]
        print "FirstFertilityAge:\t", FirstFertilityAge, "\t\t", agestopull[FirstFertilityAge]
        print "LastFertilityAge:\t", LastFertilityAge, "\t\t", agestopull[LastFertilityAge]
        print "MaxImmigrantAge:\t", MaxImmigrantAge, "\t\t", agestopull[MaxImmigrantAge]
        print "FirstDyingAge:\t\t", FirstDyingAge, "\t\t", agestopull[FirstDyingAge]
        print "All ages used from data:", agestopull

    return LeaveHouseAge, FirstFertilityAge, LastFertilityAge, MaxImmigrantAge, FirstDyingAge, agestopull

def plotDemographics(ages, datasets, I, S, T, I_touse, T_touse = None, compare_across = "T", data_year = 0):
    """
        Description:
            - Displays two plots that show the following:
                plot1: For each country: Mortality, fertility, and immigration rates, 
                                         initial and steady state population shares, 
                                         and the transition path of the total population

                plot2: Will show 2 different plots depending on the input value 'compare_across'.
                         If compare_across == "T": plot2 will display a plot of Nhat for each year in the function input T_touse with each country 
                         If compare_across == "I": plot2 will display a plot for each country of Nhat in whatever year in input value 'data_year' is

        Inputs:
            - ages                  = tuple: Contains FirstFertilityAge, LastFertilityAge, FirstDyingAge, and MaxImmigrantAge from the OLG class
            - datasets              = tuple: Contains the arrays FertilityRates, MortalityRates, ImmigrationRates, and Nhat from the OLG class
            - I                     = Int: Number of Countries
            - S                     = Int: Number of Cohorts
            - T                     = Int: Number of the total amount of time periods
            - I_touse               = List: [I], Roster of countries that are being used
            - T_touse               = List: [Unknown], List of years in plot2 given from user input
            - compare_across        = String: (Either "T" or "I"), changes the output of plot2 (see function description)
            - data_year             = Int: The year plot1 will display for each countries demographic data

        Variables Called from Object:
            - FirstFertilityAge     = Int: First age where agents give birth
            - LastFertilityAge      = Int: Last age where agents give birth
            - FirstDyingAge         = Int: First age where mortality rates effect agents
            - MaxImmigrantAge       = Int: No immigration takes place for cohorts older than this age
            - FertilityRates        = Array: [I,S,T], Fertility rates for each country, cohort, and year
            - MortalityRates        = Array: [I,S,T], Mortality rates for each country, cohort, and year
            - ImmigrationRates      = Array: [I,S,T], Immigration rates for each country, cohort, and year
            - Nhat                  = Array: [I,S,T], Population shares for each country, cohort, and year

        Other Functions Called:
            -None

        Objects in Function:
            - subplotdim_dict       = Dictionary: [6], Contains keys for each int 2-7 that maps to a dimensionality of the subplots of plot2
            - magic_int             = Int: Value from subplotdim_dict that indicates the number and dimensionality of subplots of plot2

        Outputs:
            - None
    """

    LeaveHouseAge, FirstFertilityAge, LastFertilityAge, FirstDyingAge, MaxImmigrantAge = ages
    FertilityRates, MortalityRates, ImmigrationRates, Nhat, Kids = datasets

    if T_touse is None or T_touse == "default":
        T_touse = [0, S//4, S//2, S, T]

    #The firstPlot is a subplot of key demographic data and the population dynamics for each country
    def firstPlot():
        plt.subplot(231)
        for i in range(I):
            plt.plot(range(FirstDyingAge, S-1), MortalityRates[i,FirstDyingAge:-1,data_year])
        plt.title("Mortality Rates", fontsize=14)
        plt.xlabel('Age')
        plt.ylabel('Mortality Rate')


        plt.subplot(232)
        for i in range(I):
            plt.plot(range(FirstFertilityAge,LastFertilityAge+1), FertilityRates[i,FirstFertilityAge:LastFertilityAge+1,data_year])
        plt.legend(I_touse, prop={'size':11}, loc="upper right")
        plt.title("Fertility Rates", fontsize=14)
        plt.xlabel('Age')
        plt.ylabel('Fertility Rate')


        plt.subplot(233)
        for i in range(I):
            plt.plot(range(MaxImmigrantAge), ImmigrationRates[i,:MaxImmigrantAge,data_year])
        plt.title("Immigration Rates", fontsize=14)
        plt.xlabel('Age')
        plt.ylabel('Immigration Rate')

        #Kids
        plt.subplot(234)
        for i in range(I):
            plt.plot(range(FirstFertilityAge,LastFertilityAge+LeaveHouseAge), Kids[i,FirstFertilityAge:LastFertilityAge+LeaveHouseAge,data_year])
        plt.xlabel('Age')
        plt.ylabel('Kids')
        plt.title("Kids", fontsize=14)

        plt.subplot(235)
        for i in range(I):
            plt.plot(range(S), Nhat[i,:,data_year])
        plt.xlabel('Age')
        plt.ylabel('Population Share')
        plt.title("Initial Population Shares", fontsize=14)

        #Transition path for total population of each country from the initial shares to the steady-state
        plt.subplot(236)           
        for i in range(I):
            plt.plot(range(T), np.sum(Nhat[i,:,:T], axis=0))
        plt.title("Total Pop Shares Transition Path", fontsize=14)
        plt.xlabel('Year')
        plt.ylabel('Total Population Share')

        plt.show()

    #The secondPlot compares population shares across years or across countries depending on the function input 'compare_across'
    def secondPlot():

        #Dictionary that contains ideal dimensions of the subplot depending on how many countries are being used
        subplotdim_dict = {2:222, 3:222, 4:222, 5:232, 6:232, 7:242}

        #If we want to compare each country in a given year...
        if compare_across == "T":

            if len(T_touse) == 1:
                for i in range(I):
                    plt.plot(range(S), Nhat[i,:,T_touse[0]])
                if T_touse[0] < 0:
                    T_touse[0] += T+1
                plt.title("Time t =" + str(T_touse[0]), fontsize=14)
                plt.xlabel('Age')
                plt.ylabel('Population Share')
                plt.legend(I_touse, loc="upper right")

            else:

                if len(T_touse) > 8: raise ValueError("Too many years to plot")
                magic_int = subplotdim_dict[len(T_touse)]

                plt.subplot(magic_int-1)
                for i in range(I): plt.plot(0,0)
                plt.legend(I_touse, loc="center")
                plt.gca().axes.get_xaxis().set_visible(False)
                plt.gca().axes.get_yaxis().set_visible(False)

                for count, t in enumerate(T_touse):

                    plt.subplot(magic_int+count)
                    for i in range(I):
                        plt.plot(range(S), Nhat[i,:,t])
                        if t < 0:
                            t += T+1
                        plt.title("Time t =" + str(t), fontsize=14)
                        plt.xlabel('Age')
                        plt.ylabel('Population Share')

        #If we want to compare years for each country...
        elif compare_across == "I":
            magic_int = subplotdim_dict[I]

            plt.subplot(magic_int-1)
            for t in range(T): plt.plot(0,0)
            legend = ["t = " + str(T_touse[index]) if t >= 0 else "t = " + str(T+1+T_touse[index]) for index, t in enumerate(T_touse)]
            plt.legend(legend, loc="center")
            plt.gca().axes.get_xaxis().set_visible(False)
            plt.gca().axes.get_yaxis().set_visible(False)

            for i in range(I):
                plt.subplot(magic_int+i)

                for t in T_touse:
                    plt.plot(range(S), Nhat[i,:,t])
                    plt.title(I_touse[i], fontsize=14)
                    plt.xlabel('Age')
                    plt.ylabel('Population Share')
        
        else: raise TypeError(compare_across + " is not a valid name for 'compare_across'. Choose either 'T' or 'I'")

        plt.tight_layout()

        plt.show()

    firstPlot()
    secondPlot()
