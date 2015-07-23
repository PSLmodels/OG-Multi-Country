import numpy as np
import scipy as sp
import AlphaFunctions as co
import WorldModule as WM
import cython_speedy as se
import time

"""

Governments KEY
0-United States of America
1-European Union
2-Japan
3-People's Republic of China
4-India
5-Russian Federation
6-South Korea
INSERT NEW COUNTRIES HERE
"""
functiontimer=np.zeros(30)

#print Agg_Assets_World.shape
#print Agg_Assets_World

start=time.time()
#print "C version", se.testy(5,5)
functiontime=time.time()-start
#print functiontime
co.Population_Development()
#co.PopDevelopment()
co.Initialize_Variables()


WM.IRUN=0
WM.NITER=1
print co.Sum_Wage(30,20,0,1)

#co.Initialize()
co.SteadyState()
