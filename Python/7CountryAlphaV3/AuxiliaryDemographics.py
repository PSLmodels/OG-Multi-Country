from __future__ import division
import csv
import numpy as np
import scipy as sp
import scipy.optimize as opt
import time as time
from matplotlib import pyplot as plt


#DEMOGRAPHICS FUNCTIONS
def getkeyages(S, PrintAges, UseStaggeredAges):
    """
    Description:
        -Gets the key ages for calculating demographics based on S

    Inputs:
        S                = Int in [10,80], Number of cohorts
        PrintAges        = Bool, Prints key ages if true
        UseStaggeredAges = Bool, Gets ages to use from a linspace rather than a sequence of integers

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
    if UseStaggeredAges:
        agestopull = np.arange(80)[np.where(np.in1d(np.arange(80), np.round(np.linspace(0, 80, num=S, endpoint=False))))]
        LeaveHouseAge = np.abs(agestopull - 21).argmin()#The age when agents become independent households and don't rely on parents consumption
        FirstFertilityAge = np.abs(agestopull - 23).argmin()#The age when agents have their first children
        LastFertilityAge = np.abs(agestopull - 45).argmin()#The age when agents have their last children
        MaxImmigrantAge = np.abs(agestopull - 65).argmin()#All immigrants are between ages 0 and MaxImmigrantAge
        FirstDyingAge = np.abs(agestopull - 68).argmin()#The first age agents can begin to die
    else:
        agestopull = range(S)
        LeaveHouseAge = int(np.round(S/80.*21))#The age when agents become independent households and don't rely on parents consumption
        FirstFertilityAge = int(np.round(S/80.*23))#The age when agents have their first children
        LastFertilityAge = int(np.round(S/80.*45))#The age when agents have their last children
        MaxImmigrantAge = int(np.round(S/80.*65))#All immigrants are between ages 0 and MaxImmigrantAge
        FirstDyingAge = int(np.round(S/80.*68))#The first age agents can begin to die

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
        print "All ages to pull:", agestopull
        print "Adjusted to S =", S, ":    LeaveHouseAge =", LeaveHouseAge, "FirstFertilityAge =", FirstFertilityAge, \
        "LastFertilityAge =", LastFertilityAge, "FirstDyingAge =", FirstDyingAge, "MaxImmigrantAge =", MaxImmigrantAge
        if UseStaggeredAges:
            print "To pull from the data:  LeaveHouseAge =", agestopull[LeaveHouseAge], "FirstFertilityAge =", agestopull[FirstFertilityAge], \
            "LastFertilityAge =", agestopull[LastFertilityAge], "FirstDyingAge =", agestopull[FirstDyingAge], "MaxImmigrantAge =", agestopull[MaxImmigrantAge],"\n"

    return LeaveHouseAge, FirstFertilityAge, LastFertilityAge, MaxImmigrantAge, FirstDyingAge, agestopull

def plotDemographics(ages, datasets, I, S, T, I_touse, T_touse = None, compare_across = "T", data_year = 0):

    FirstFertilityAge, LastFertilityAge, FirstDyingAge, MaxImmigrantAge = ages
    FertilityRates, MortalityRates, ImmigrationRates, Nhat = datasets


    if T_touse is None or T_touse == "default":
        T_touse = [0, S//4, S//2, S, T]

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

        plt.subplot(234)
        for i in range(I):
            plt.plot(range(S), Nhat[i,:,data_year])
        plt.xlabel('Age')
        plt.ylabel('Population Share')
        plt.title("Initial Population Shares", fontsize=14)

        plt.subplot(235)
        for i in range(I):
            plt.plot(range(S), Nhat[i,:,-1])
        plt.xlabel('Age')
        plt.ylabel('Population Share')
        plt.title("Steady State Population Shares", fontsize=14)

        plt.subplot(236)           
        for i in range(I):
            plt.plot(range(T), np.sum(Nhat[i,:,:T], axis=0))
        plt.title("Total Pop Shares Transition Path", fontsize=14)
        plt.xlabel('Year')
        plt.ylabel('Total Population Share')

        plt.show()

    def secondPlot():

        subplotdim_dict = {2: (122, False), 3:(222, False), 4:(222, False), 5:(232, True), 6:(232, True), 7:(242, True), 8:(242, True)}

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
                full_screen = False

            else:

                if len(T_touse) > 8: raise ValueError("Too many years to plot")
                magic_int = subplotdim_dict[len(T_touse)][0]
                full_screen = subplotdim_dict[len(T_touse)][1]

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

        elif compare_across == "I":
            magic_int = subplotdim_dict[I][0]
            full_screen = subplotdim_dict[I][1]

            plt.subplot(magic_int-1)
            for i in range(I): plt.plot(0,0)
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
        
        else: raise TypeError(compare_across + " is not a valid name for 'compare_across'")

        plt.tight_layout()
        #if full_screen: plt.get_current_fig_manager().window.showMaximized()
        plt.show()

    firstPlot()
    secondPlot()