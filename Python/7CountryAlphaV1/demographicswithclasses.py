import numpy as np
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

GENS = 90
FIRST_DEATH_AGE = 68
OLDEST_IMMIGRANTS = 65
FIRST_FERTILITY_YEAR = 23
LAST_FERTILITY_YEAR = 45
FIRST_WORK_YEAR = 21
SKILL_GROUPS = 2

"""
TODO:
#net_migration only goes to age 65?
#Constant migration each year?
#Low skill babies become low skill parents?
#Make some paramater like Maxyears = 300 that is the furthest number of years the model can go
#Look at paper to compare projections past 2050
#Make delta change for each Region, or import it from somewhere
#Make sure all the dimentions of data in here correspond to WorldModel (fertility is different for example)
#Update data readin files format
"""


def Region_information():
		"""
		********************************************************************************************************************
		REGION CLASS DATA TYPES:
		self.name: string
			-name of the region, eg. "Russia", "EU"
			-Used when generically printing a Region and to gather data from the folder Data_files

		self.index: integer
			-USA = 1, EU = 2, Japan = 3, China = 4, India = 5, Russia = 6, Korea = 7
			-Used also to gather data from outside .csv files, specifically from population.csv and net_migration.csv

		self.initial_population
			-Vector of length GENS+1 for ages 0-GENS+1 that contains exact number of people in each
			-Comes from population.csv for specific country
			-Note: population.csv is in thousands of people, while self.initial_population is in exact number of people

		self.net_migration
			-Vector of length 65 for ages 1-65 that has the total number of immigrants for per year for each age group
			-Comes from net_migration.csv for specific country
			-Note: net_migration.csv is in hundreds of people, while self.net_migration is in exact number of people

		self.fertility_rates
			-Has 23 rows for ages 23-45 (when agents give birth) and 51 columns for years 2008-2058
			-Essentially the transpose of the data in self.name_fertility.csv for nonnegative years

		self.mortality_rates
			-Has 23 rows for ages 68-90 and 51 columns for years 2008-2058
			-Essentially the transpose of the data in self.name_mortality.csv

		self.skill_distributions
			-Vector of length 2 that has the percentage of total population in each skill group

		self.population_timepath: np.array((GENS, number of simulated years, SKILL_GROUPS))
			-91 rows for ages 0-90, Has 1 column for year 2008, and 2 depth for skill classes 0-1. 
			-Will add more columns for additional years simulated
			-Stores the population size for each generation, year, and skill group

		self.KID_mat: np.array((GENS, number of simulated years, SKILL_GROUPS))
			-91 rows for ages 0-90, Has 1 column for year 2008, and 2 depth for skill classes 0-1. 
			-Will add more columns for additional years simulated
			-Stores the number of kids each agent of that generation has. This implies the following:
				-Ages 0-22 = 0
				-Ages 45-68 = A constant because they have kids from their fertile years, but aren't having new kids
				-Ages 68-90 = 0

		self.timeendowment: np.array(number of simulated years)
			-Stores the time endowment for each year. This is h(a, i) in the paper

		self.delta: double
			-Technilogical progress

		REGION CLASS FUNCTIONS:
		readindata(): 

		newKids(year): 
			-Creates a vector of kids length 91x2 that is added to self.KID_mat for each year

		simulate_demographics(Years): 
			-Simulates population changes for "Years" years
			-Kills old people, makes new babies, adds immigration, gets a new time endowment, and adds new kids and generation to their respective matrices

		plot_demographics(year): 
			-Plots Population distribution, Fertility, Mortality, and Net Migration for the inputed year across ages

		plot_population_distribution(years): 
			-Takes in a list of years and plots the Population Distribution of those years on the same graph

		plot_total_population(year): 
			-Plots the change in total population over time from 2008 to the inputed year

		get_total_population(year):
			-Returns the changes in totalpopulation up until the given year

		get_population_distribution(year, returnclasses = False):
			-Returns the population distribution for a given year. If returnclasses = True it returns the population with the extra dimension for classes
		
		get_fertility_rate(year):
			-Returns the fertility rate for a given year

		get_mortality_rate(year, returnall = False):
			-Returns mortality rate for a given year
			-If returnall == True, it returns the whole array of mortality rates

		get_total_netmigration_rate():
			-Returns the net migration data across age
			-As of now it is the same for each year

		get_kids_mat(year, returnall = False):
			-Returns the distribution of kids taken from KID_mat for the given year
			-If returnall == True, it returns the entire KID_mat

			
		******************************************************************************************************************
		"""
class Region(object):

	def __init__(self, name, index):

		self.name = name
		self.index = index
		self.delta = .01

		def readindata():

			with open('Data_files/Backup/population.csv','r') as csv_file:
				csv_reader=csv.reader(csv_file)
				popdata = []
				for age in csv_reader:
					for region in csv_reader:
						popdata.append(region[self.index])
				popdata = np.array(popdata).astype(np.float)*1000

			with open('Data_files/Backup/skillclasses.csv','r') as csv_file:
				csv_reader=csv.reader(csv_file)
				skilldata = []
				for skill in csv_reader:
					for region in csv_reader:
						skilldata.append(region[self.index])

				skilldata = np.array(skilldata).astype(np.float)

			with open("Data_files/Backup/"+(self.name).lower()+"_fertility.csv","r") as csv_file:
				csv_reader=csv.reader(csv_file)
				fertdata = []
				for age in csv_reader:
					for year in csv_reader:
						fertdata.append(year[1:])
				fertdata = np.transpose(np.array(fertdata).astype(np.float)[48:,])

			with open("Data_files/Backup/"+str(self.name).lower()+"_mortality.csv","r") as csv_file:
				csv_reader=csv.reader(csv_file)
				mortdata = []
				for age in csv_reader:
					for year in csv_reader:
						mortdata.append(year[1:])
				mortdata = np.transpose(np.array(mortdata).astype(np.float))

			with open('Data_files/Backup/net_migration.csv','r') as csv_file:
				csv_reader=csv.reader(csv_file)
				migdata = []
				for age in csv_reader:
					for region in csv_reader:
						migdata.append(region[self.index])
				migdata = np.array(migdata).astype(np.float)*100 
			
			#Takes initial population and migration and give them an 2nd dimention of length 2 for each skill group

			popskilldata = np.zeros((GENS+1, SKILL_GROUPS))
			migskilldata = np.zeros((OLDEST_IMMIGRANTS, SKILL_GROUPS))

			for k in range(SKILL_GROUPS):
				popskilldata[:,k] = popdata*skilldata[k]
				migskilldata[:,k] = migdata*skilldata[k]

			return popskilldata, migskilldata, fertdata, mortdata, skilldata

		self.initial_population, self.net_migration, self.fertility_rates, self.mortality_rates, self.skill_distributions = readindata()
		self.population_timepath = np.zeros((91,1,2))
		self.population_timepath[:,0,:] = self.initial_population
		self.KID_mat = self.newKids(2008).reshape((91, 1, 2))
		self.timeendowment = np.ones(1)

	def __repr__(self):
		return self.name

	def newKids(self, year):

		if year > 2050:
			#We only have fertility rates after the year 2050
			year = 2050

		#Column that contains the number of kids each fertile age group has (ages 23-45 in this case)
		fertilenumkids = np.cumsum(self.fertility_rates[0:FIRST_FERTILITY_YEAR+1, year-2008])
		
		#Combines the number of kids for age groups 0-22, 23-45, 45-65, and 66-90 into one vector of length 91
		numkids = np.hstack((np.zeros(23) , fertilenumkids , np.ones((20))*fertilenumkids[-1] , np.zeros(25)))

		#Adding a column of numkids for each skill group
		kidsvec = np.tile(numkids.reshape(GENS+1, 1), (1, SKILL_GROUPS))

		return kidsvec

	def simulate_demographics(self, Years):

		#Years: Number of years to simulate the population changes

		population_vec = self.population_timepath[:,-1,:]#A 91x2 array for each age and skill class for the most recent year
		newpopulation_vec = np.zeros((91,2))
		Transyear = 50#Number of years for which we have data. Anything beyond this will using data for the 50th year

		for t in xrange(1, Years):
			if t <= Transyear:
				i = t
			elif t > Transyear:
				i = Transyear

			#Old People dying in each skill class
			for k in range(2):
				newpopulation_vec[68:,k] = population_vec[68:,k] - population_vec[68:,k]*self.mortality_rates[:,i]

			#Getting the total number of babies born this year
			num_babies = np.dot(np.transpose(population_vec[23:46,:]),self.fertility_rates[:,i])

			#Adding newborns as youngest generation and killing the oldest generation by only stacking until the last generation of the population_vec
			newpopulation_vec = np.vstack((np.transpose(num_babies),population_vec[:-1,:]))
			
			#Adding immigration
			newpopulation_vec[1:66,] += self.net_migration

			#Storing the population distribution in a timepath matrix
			self.population_timepath = np.column_stack((self.population_timepath,np.reshape(newpopulation_vec, (91, 1, 2))))

			#Re-initializing the previous year as the current year
			population_vec = newpopulation_vec

			#Applying the KID function to the current year and storing it in the timepath for KIDS
			self.KID_mat = np.column_stack((self.KID_mat, np.reshape(self.newKids(t+2008), (91, 1, 2))))

			#Increase time endowment in line with technical progress
			self.timeendowment = np.append(self.timeendowment, (1+self.delta)*self.timeendowment[-1])
			

	def plot_demographics(self, year):

		#IMPORTANT!! PLOTS THE SUM OF THE SKILL CLASSES. THATS' WHAT THE axis=1 NONSENSE IS

		num_simulated_years = self.population_timepath.shape[1]
		if year - 2008 >= num_simulated_years:
			print "\nERROR: WE HAVE ONLY SIMULATED UP TO THE YEAR", num_simulated_years+2008, "AND YOU REQUESTED TO PLOT DATA FOR THE YEAR", year
			print"*SEE plot_demog_distribution IN class Region(object)*\n"
			time.sleep(10)
			return None

		year_index = year - 2008
		plt.clf()
		plt.suptitle(str(self.name+" Data for "+str(year)))

		plt.subplot(2, 2, 1)
		plt.plot(range(91),self.population_timepath[:,year_index,:].sum(axis=1))
		plt.title("Population Distribution")
		plt.grid()

		plt.subplot(2, 2, 2)
		plt.plot(range(23,46),self.get_fertility_rate(year))
		plt.xlim(23, 46)
		plt.title("Fertility Rates")
		plt.grid()

		plt.subplot(2, 2, 3)
		plt.plot(range(68,91), self.get_mortality_rate(year))
		plt.xlim(68, 89)
		plt.ylim(0, 1)
		plt.title("Mortality Rates")
		plt.grid()

		plt.subplot(2, 2, 4)
		plt.plot(range(65), self.get_total_netmigration_rate())
		plt.title("Total Net Migration")
		plt.grid()

		plt.show()
		plt.clf()

	def plot_population_distribution(self, years):
		years = np.array(years)

		for y in range(len(years)):
			yeartograph = years[y]
			num_simulated_years = self.population_timepath.shape[1]

			if yeartograph - 2008 < num_simulated_years:
				#IMPORTANT! Plots the sum of the skill classes for each year
				plt.plot(range(91), self.population_timepath[:,yeartograph-2008,:].sum(axis=1))
			else:
				print "\nERROR: WE HAVE ONLY SIMULATED UP TO THE YEAR", num_simulated_years+2008-1, "AND YOU REQUESTED TO PLOT THE YEAR", yeartograph
				print"*SEE plot_population_distribution IN class Region(object)*\n"
				time.sleep(15)
				return None

		plt.title(str(self.name + " Population Distribution"))
		plt.legend(years)
		plt.show()

	def plot_total_population(self, year):

		totalpopulation = self.population_timepath.sum(axis=0).sum(axis=1)[:year-2008+1]

		plt.plot(range(2008, year+1), totalpopulation/1000000)
		plt.title(self.name+" Population Change from 2008-"+ str(year))
		plt.xlim(2008, year)
		plt.xlabel('Year')
		plt.ylabel('Population (Millions)')
		plt.show()

	def get_total_population(self, year):
		totalpopulation = self.population_timepath.sum(axis=0).sum(axis=1)[:year-2008+1]
		return totalpopulation

	def get_population_distribution(self, year, returnclasses = False):
		if returnclasses == True:
			if year-2008 < self.population_timepath.shape[1]:
				return self.population_timepath[:,year-2008,:]
			else:
				print "We have only calculated up till year", self.population_timepath.shape[0], "so we are returning data for that year"
				return self.population_timepath[-1,:]
		else: #if returnclasses == False
			if year-2008 < self.population_timepath.shape[1]:
				return self.population_timepath[:,year-2008,:].sum(axis=1)
			else:
				print "We have only calculated up till year", self.population_timepath.shape[0], "so we are returning data for that year"
				return self.population_timepath[-1,:].sum(axis=1)

	def get_fertility_rate(self,year, returnall = False):

		if returnall == True:
			return self.fertility_rates

		if year-2008 < self.fertility_rates.shape[1]:
			return self.fertility_rates[:,year-2008]
		else:
			print "\nThis data is too far out, so we are returning the steady-state value\n"
			return self.fertility_rates[:,-1]

	def get_mortality_rate(self, year, returnall = False):

		if returnall == True:
			return self.mortality_rates

		if year-2008 < self.mortality_rates.shape[1]:
			return self.mortality_rates[:,year-2008]
		else:
			print "\nThis data is too far out, so we are returning the steady-state value\n"
			return self.mortality_rates[:,-1]

	def get_total_netmigration_rate(self):
		return self.net_migration.sum(axis=1)

	def get_kids_mat(self, year, returnall = False):
		if returnall == True:
			return self.KID_mat
		return self.KID_mat[:, year-2008, :]


def compare_countries(countries, year):

	plt.clf()
	plt.suptitle("Demographics "+str(year))
	
	plt.subplot(1, 2, 1)
	for r, region in enumerate(countries):
		plt.plot(range(91),region.get_population_distribution(year))
	plt.title("Age Distribution")
	plt.xlabel('Age')
	plt.ylabel('Population')
	plt.grid()

	plt.subplot(1,2,2)
	for r, region in enumerate(countries):
		plt.plot(range(2008, year+1), region.get_total_population(year)/1000000)
	plt.title("Total Population (Millions)")
	plt.xlim(2008, year)
	plt.xlabel('Year')
	plt.legend(countries, loc = "upper left", prop={'size':11})
	plt.grid()
	plt.show()
	plt.clf()

	plt.suptitle(" Data for the Year "+str(year))
	plt.legend(regionlist, loc = "lower right")

	plt.subplot(2, 2, 1)
	for r, region in enumerate(countries):
		plt.plot(range(23,46),region.get_fertility_rate(year))
	plt.title("Fertility")
	plt.xlim(23, 46)
	plt.grid()
		
	plt.subplot(2, 2, 2)
	for r, region in enumerate(countries):
		plt.plot(range(68,91),region.get_mortality_rate(year))
	plt.title("Mortality")
	plt.xlim(68, 89)
	plt.ylim(0, 1)
	plt.legend(countries, loc = "upper left", prop={'size':11})
	plt.grid()

	plt.subplot(2, 2, 3)
	for r, region in enumerate(countries):
		plt.plot(range(65),region.get_total_netmigration_rate())
	plt.title("Net Migration")
	plt.grid()

	plt.show()
	plt.clf()


USA = Region("USA", 1)
EU = Region("EU", 2)
Japan = Region("Japan", 3)
China = Region("China", 4)
India = Region("India", 5)
Russia = Region("Russia", 6)
Korea = Region("Korea", 7)

regionlist = [USA, EU, Japan, China, India, Russia, Korea]

for index, region in enumerate(regionlist):
	region.simulate_demographics(300)


#Russia.plot_demographics(2058)


USA.plot_population_distribution([2008, 2013])




print "Done"