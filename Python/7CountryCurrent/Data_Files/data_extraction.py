import pandas as pd 

# ****************************Extracting the Population******************************
toskip=range(0,1)+range(94,3520)
data = pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
data['JAPAN']=data['JAPAN']-data['KOREA']
data.to_csv("population.csv", index=False)

# ****************************Extracting the Net Migration****************************
toskip=range(0,95)+range(164, 3520)
netm=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="      ", skiprows=toskip)
netm.to_csv("net_migration.csv", index=False)

# *******************************Extract Fertility USA********************************
toskip=range(0,165)+range(266, 3520)
USfert=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,266)+range(367,3520)
USfert1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,367)+range(468,3520)
USfert2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
USfert1=USfert1.drop('YEAR',1)
USfert2=USfert2.drop('YEAR',1)

#concat
totfert=pd.concat([USfert,USfert1,USfert2],join='outer',axis=1)
totfert.to_csv("usa_fertility.csv", index=False)


# *******************************Extract Fertility EU********************************
toskip=range(0,469)+range(570, 3520)
fert=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,570)+range(671,3520)
fert1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,671)+range(772,3520)
fert2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
fert1=fert1.drop('YEAR',1)
fert2=fert2.drop('YEAR',1)

#concat
totfert=pd.concat([fert,fert1,fert2],join='outer',axis=1)
totfert.to_csv("eu_fertility.csv", index=False)

# *****************************Japan Extract Fertility *****************************
toskip=range(0,774)+range(875, 3520)
fert=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,875)+range(976,3520)
fert1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,976)+range(1077,3520)
fert2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
fert1=fert1.drop('YEAR',1)
fert2=fert2.drop('YEAR',1)

#concat
totfert=pd.concat([fert,fert1,fert2],join='outer',axis=1)
totfert.to_csv("japan_fertility.csv", index=False)

# *****************************China Extract Fertility *****************************
toskip=range(0,1079)+range(1180, 3520)
fert=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,1180)+range(1281,3520)
fert1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,1281)+range(1382,3520)
fert2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
fert1=fert1.drop('YEAR',1)
fert2=fert2.drop('YEAR',1)

#concat
totfert=pd.concat([fert,fert1,fert2],join='outer',axis=1)
totfert.to_csv("china_fertility.csv", index=False)


# *****************************India Extract Fertility *****************************
toskip=range(0,1384)+range(1485, 3520)
fert=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,1485)+range(1586,3520)
fert1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,1586)+range(1687,3520)
fert2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
fert1=fert1.drop('YEAR',1)
fert2=fert2.drop('YEAR',1)

#concat
totfert=pd.concat([fert,fert1,fert2],join='outer',axis=1)
totfert.to_csv("india_fertility.csv", index=False)

# *****************************Russia Extract Fertility *****************************
toskip=range(0,1688)+range(1789, 3520)
fert=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,1789)+range(1890,3520)
fert1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,1890)+range(1991,3520)
fert2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
fert1=fert1.drop('YEAR',1)
fert2=fert2.drop('YEAR',1)

#concat
totfert=pd.concat([fert,fert1,fert2],join='outer',axis=1)
totfert.to_csv("russia_fertility.csv", index=False)

# *****************************Korea Extract Fertility *****************************
toskip=range(0,1993)+range(2094, 3520)
fert=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2094)+range(2195,3520)
fert1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2195)+range(2296,3520)
fert2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
fert1=fert1.drop('YEAR',1)
fert2=fert2.drop('YEAR',1)

#concat
totfert=pd.concat([fert,fert1,fert2],join='outer',axis=1)
totfert.to_csv("korea_fertility.csv", index=False)




# *****************************USA Extract Mortality*****************************
toskip=range(0,2299)+range(2351, 3520)
mort=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2352)+range(2404,3520)
mort1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2405)+range(2457,3520)
mort2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
mort1=mort1.drop('YEAR',1)
mort2=mort2.drop('YEAR',1)

#concat
totmort=pd.concat([mort,mort1,mort2],join='outer',axis=1)
totmort.to_csv("usa_mortality.csv", index=False)

# *****************************EU Extract Mortality*****************************
toskip=range(0,2460)+range(2512, 3520)
mort=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2513)+range(2565,3520)
mort1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2566)+range(2618,3520)
mort2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
mort1=mort1.drop('YEAR',1)
mort2=mort2.drop('YEAR',1)

#concat
totmort=pd.concat([mort,mort1,mort2],join='outer',axis=1)
totmort.to_csv("eu_mortality.csv", index=False)

# *****************************Japan Extract Mortality*****************************
toskip=range(0,2621)+range(2673, 3520)
mort=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2674)+range(2726,3520)
mort1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2727)+range(2779,3520)
mort2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
mort1=mort1.drop('YEAR',1)
mort2=mort2.drop('YEAR',1)

#concat
totmort=pd.concat([mort,mort1,mort2],join='outer',axis=1)
totmort.to_csv("japan_mortality.csv", index=False)

# *****************************China Extract Mortality*****************************
toskip=range(0,2782)+range(2834, 3520)
mort=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2835)+range(2887,3520)
mort1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2888)+range(2940,3520)
mort2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
mort1=mort1.drop('YEAR',1)
mort2=mort2.drop('YEAR',1)

#concat
totmort=pd.concat([mort,mort1,mort2],join='outer',axis=1)
totmort.to_csv("china_mortality.csv", index=False)

# *****************************India Extract Mortality*****************************
toskip=range(0,2943)+range(2995, 3520)
mort=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,2996)+range(3048,3520)
mort1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,3049)+range(3101,3520)
mort2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
mort1=mort1.drop('YEAR',1)
mort2=mort2.drop('YEAR',1)

#concat
totmort=pd.concat([mort,mort1,mort2],join='outer',axis=1)
totmort.to_csv("india_mortality.csv", index=False)

# *****************************Russia Extract Mortality*****************************
toskip=range(0,3104)+range(3156, 3520)
mort=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,3157)+range(3209,3520)
mort1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,3210)+range(3262,3520)
mort2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
mort1=mort1.drop('JAHR',1)
mort2=mort2.drop('JAHR',1)

#concat
totmort=pd.concat([mort,mort1,mort2],join='outer',axis=1)
totmort.to_csv("russia_mortality.csv", index=False)

# *****************************Korea Extract Mortality*****************************
toskip=range(0,3265)+range(3317, 3520)
mort=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,3318)+range(3370,3520)
mort1=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)
toskip=range(0,3371)+range(3423,3520)
mort2=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="     ", skiprows=toskip)

#to concat, we drop the year variables
mort1=mort1.drop('YEAR',1)
mort2=mort2.drop('YEAR',1)

#concat
totmort=pd.concat([mort,mort1,mort2],join='outer',axis=1)
totmort.to_csv("korea_mortality.csv", index=False)


# *********************Healthcare Spending******************************************
toskip=range(0,3425)+range(3517,3520)
hc=pd.io.parsers.read_table("POP_Korea_06_10_14_test2.dat",sep="      ", skiprows=toskip)

hc.to_csv("healthcare.csv")