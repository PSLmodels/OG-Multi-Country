import numpy as np
import Main as mn
import sys
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE
import time

"""
This piece of code is intented to work with the 7 Country OLG model. It cycles through various combinations
of cohort levels and numbers of countries and sees where the code works and where it breaks. It outputs two
sets of sheets: one set is CSV files of booleans and the other set is CSV files of times of how long each 
combination took.

"""

comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()

Sinputs=np.array([20])
S=len(Sinputs)


IUB=8
SigUB=6

Sigstart=1
Istart=2

checker=np.zeros(((SigUB-Sigstart)*(IUB-Istart)*(S)),dtype=bool) if rank==0 else None
timer=np.zeros(((SigUB-Sigstart)*(IUB-Istart)*(S)),dtype=float) if rank==0 else None

if rank!=0:
    checkerpiece=np.zeros(((IUB-Istart)*(S)),dtype=bool)
    timerpiece=np.zeros(((IUB-Istart)*(S)))
else: 
    checkerpiece=np.zeros(0,dtype=bool)
    timerpiece=np.zeros(0)

comm.Scatterv([checker,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.BOOL],checkerpiece)

comm.Scatterv([timer,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.DOUBLE],timerpiece)

if rank!=0:
    sig=rank
    for s in xrange(S):
        for i in xrange(Istart,IUB):
            try:
                start=time.time()
                mn.Multi_Country(Sinputs[s],i,sig)
                tottime=time.time()-start
                checkerpiece=np.reshape(checkerpiece,(IUB-Istart,S))
                checkerpiece[i-Istart,s]=True
                checkerpiece=np.reshape(checkerpiece, (IUB-Istart)*(S))
                timerpiece=np.reshape(timerpiece,(IUB-Istart,S))
                timerpiece[i-Istart,s]=tottime
                timerpiece=np.reshape(timerpiece, (IUB-Istart)*(S))
                
                print "Success at:", Sinputs[s], "cohorts,",i,"countries",sig,"slope and took", tottime, "seconds to complete"
            except:
                checkerpiece=np.reshape(checkerpiece,(IUB-Istart,S))
                checkerpiece[i-Istart,s]=False
                checkerpiece=np.reshape(checkerpiece, (IUB-Istart)*(S))
                timerpiece=np.reshape(timerpiece,(IUB-Istart,S))
                timerpiece[i-Istart,s]=0
                timerpiece=np.reshape(timerpiece, (IUB-Istart)*(S))

                print "Failure at:", Sinputs[s], "cohorts,",i,"countries",sig,"slope"

comm.Barrier()

comm.Gatherv(checkerpiece,[checker,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.BOOL])

comm.Gatherv(timerpiece,[timer,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.DOUBLE])


if rank==0:
    checker=np.reshape(checker,(SigUB-Sigstart,IUB-Istart,S))
    for z in xrange(SigUB-Sigstart):
        label="Curve_" +str(z)+".csv"
        np.savetxt(label,checker[z,:,:],delimiter=",")
if rank==0:
    timer=np.reshape(timer,(SigUB-Sigstart,IUB-Istart,S))
    for z in xrange(SigUB-Sigstart):
        label="Time_"+str(z)+".csv"
        np.savetxt(label,timer[z,:,:],delimiter=",")

