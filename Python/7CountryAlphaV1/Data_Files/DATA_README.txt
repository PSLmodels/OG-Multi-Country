This file contains documentation about each of the data files (saved as .csv). Note that from the original fortran version, I decided to clip out the labels and instead placed them here, this way I can take advantage of the np.loadtxt command that python has to make the process of loading in data speedier. The whole point of this is to cutdown on for loops, since Python isn’t as efficient at for loops as fortran. The files are listed in alphabetical order.

******************************************
FORMAT:
file_name.csv
FortranCode: Some of the csv files correspond to lines of fortran code. This is indicated here.
COLUMNS: Here is where the labels for the columns are listed
ROWS: Here is where the labels for the rows are listed.
******************************************


china_fertility.csv
FortranCode:
COLUMNS: Range of fertility ages from 23-45
ROWS:Years, indexed from -48 to 50


-----------------------------------------------------------


china_mortality.csv
FortranCode:
COLUMNS: Range of death years from 68 to 90
ROWS: Years indexed from 0 to 50

-----------------------------------------------------------


contributionceiling.csv, this inputs the the contribution ceilings for pension systems
FortranCode: 720-726
COLUMNS: Countries: USA, EU, Japan, China, India, Russia, Korea
ROWS: -1 to 300, years.

-----------------------------------------------------------


coporatetax.csv
FortranCode: 612-619
COLUMNS: Countries, USA, EU, Japan, China, India, Russia, Korea
ROWS: Years, with the initial index is -1 and the last one is 300

-----------------------------------------------------------


debtlevel.csv, this file is where the target debt level is set. At least one of the above revenues has to be endogenous in order for countries to stay at the target debt level in every year (needed for convergence)
FortranCode: 671-679
COLUMNS: USA, EU, Japan, China, India, Russia, Korea
ROWS: Years, starting at -1 and going to 300

-----------------------------------------------------------


eu_fertility.csv
FortranCode:
COLUMNS:Range of fertility ages from 23-45
ROWS:Years, indexed from -48 to 50

-----------------------------------------------------------


eu_mortality.csv
FortranCode:
COLUMNS: Range of death years from 68 to 90
ROWS: Years indexed from 0 to 50

-----------------------------------------------------------


govpy.csv, part of the pension benefits contributed by the government
FortranCode:623-629
COLUMNS:Countries, USA, EU, Japan, China, India, Russia, Korea
ROWS: Years, with the index beginning at -1 and going to 300

-----------------------------------------------------------


healthcare.csv
FortranCode:
COLUMNS:
ROWS:

-----------------------------------------------------------


india_fertility.csv
FortranCode:
COLUMNS:Range of fertility ages from 23-45
ROWS:Years, indexed from -48 to 50

-----------------------------------------------------------


india_mortality.csv
FortranCode:
COLUMNS: Range of death years from 68 to 90
ROWS: Years indexed from 0 to 50

-----------------------------------------------------------


inheritancetax.csv
FortranCode:604-610
COLUMNS:Countries, USA, EU, Japan, China, India, Russia, Korea
ROWS:Years
NOTE THAT THIS SEEMS REDUNDANT,SINCE IT’S SET UP IN PYTHON AS ALL ZEROS ANYWAY.
-----------------------------------------------------------


japan_fertility.csv
FortranCode:
COLUMNS:Range of fertility ages from 23-45
ROWS:Years, indexed from -48 to 50

-----------------------------------------------------------


japan_mortality.csv
FortranCode:
COLUMNS: Range of death years from 68 to 90
ROWS: Years indexed from 0 to 50

-----------------------------------------------------------


korea_fertility.csv
FortranCode:
COLUMNS:Range of fertility ages from 23-45
ROWS:Years, indexed from -48 to 50

-----------------------------------------------------------


korea_mortality.csv
FortranCode:
COLUMNS: Range of death years from 68 to 90
ROWS: Years indexed from 0 to 50

-----------------------------------------------------------


mu2gov.csv, part of the health benefits treated as government consumption
FortranCode:632-638
COLUMNS:Countries, USA, EU, Japan, China, India, Russia, Korea
ROWS: Years, starting index at -1 and it goes to 300

-----------------------------------------------------------


mu2tax.csv, part of the health benefits financed by general taxes
FortranCode:640-647
COLUMNS:Countries, USA EU, Japan, China, India, Russia, Korea
ROWS: Years, starting index at -1 and it goes to 300

-----------------------------------------------------------


mu3.csv, part of the disability system paid for through general taxes
FortranCode: 650-656
COLUMNS:Countries: USA, EU, Japan, China, India, Russia, Korea
ROWS: Years: index starting at -1 and goes to 300

-----------------------------------------------------------


mu4.csv, records the share of corporate tax revenues which is lump-sum transferred to HH.
FortranCode: 659-665
COLUMNS: Countries: USA, EU, Japan, China, India, Russia, Korea
ROWS: Years: index starting at -1 and goes to 300

-----------------------------------------------------------


net_migration.csv
FortranCode:
COLUMNS: Countries: USA, EU, Japan, China, India, Russia, Korea
ROWS:age from 1 to 65

-----------------------------------------------------------


popscale.csv, This section of the code allows for weighting populations. It will be later multiplied by the population read in from the demographic file
FortranCode: 190-199
Rows:Countries: USA, EU, Japan, China, India, Russia, Korea

-----------------------------------------------------------


population.csv
FortranCode:
COLUMNS: Countries: USA, EU, Japan, China, India, Russia, Korea
ROWS: Ages from 0 to 90

-----------------------------------------------------------


retirementage.csv
FortranCode: 692-727
COLUMNS=Countries: USA, EU, Japan China, India, Russia Korea
ROWS=The Year, with the index beginning at -1 and going to 300

-----------------------------------------------------------


russia_fertility.csv
FortranCode:
COLUMNS:Range of fertility ages from 23-45
ROWS:Years, indexed from -48 to 50

-----------------------------------------------------------


russia_mortality.csv
FortranCode:
COLUMNS: Range of death years from 68 to 90
ROWS: Years indexed from 0 to 50

-----------------------------------------------------------


skillclasses.csv, gives the portions of the population in each skill class
FortranCode: 203-216
COLUMNS:Countries, USA, EU, Japan, China, India, Russia, Korea
ROWS:First Row: Low Skill, Second Row: High Skill

-----------------------------------------------------------


usa_fertility.csv
FortranCode:
COLUMNS:Range of fertility ages from 23-45
ROWS:Years, indexed from -48 to 50

-----------------------------------------------------------


usa_mortality.csv
FortranCode:
COLUMNS: Range of death years from 68 to 90
ROWS: Years indexed from 0 to 50

-----------------------------------------------------------

