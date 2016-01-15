import numpy as np
import Main as mn
import sys
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE
import time

comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()

Sinputs=np.array([20,25,30,35,40,45,50,55,60,65,70,75,80])
S=13


IUB=8
SigUB=6

Sigstart=1
Istart=2

checker=np.zeros(((SigUB-Sigstart)*(IUB-Istart)*(S)),dtype=bool) if rank==0 else None
timer=np.zeros(((SigUB-Sigstart)*(IUB-Istart)*(S)),dtype=bool) if rank==0 else None

if rank!=0:
    checkerpiece=np.zeros(((IUB-Istart)*(S)),dtype=bool)
    timerpiece=np.zeros(((IUB-Istart)*(S)))
else: 
    checkerpiece=np.zeros(0,dtype=bool)
    timerpiece=np.zeros(0)

comm.Scatterv([checker,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.BOOL],checkerpiece)

comm.Scatterv([timer,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.BOOL],timerpiece)

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
                timerpiece=np.reshape(checkerpiece,(IUB-Istart,S))
                timerpiece[i-Istart,s]=tottime
                timerpiece=np.reshape(checkerpiece, (IUB-Istart)*(S))
                
                print "Success at:", Sinputs[s], "cohorts,",i,"countries",sig,"slope"
                #comm.Send(checkerpiece,dest=0)
            except:
                checkerpiece=np.reshape(checkerpiece,(IUB-Istart,S))
                checkerpiece[i-Istart,s]=False
                checkerpiece=np.reshape(checkerpiece, (IUB-Istart)*(S))
                timerpiece=np.reshape(checkerpiece,(IUB-Istart,S))
                timerpiece[i-Istart,s]="FAILED"
                timerpiece=np.reshape(checkerpiece, (IUB-Istart)*(S))

                print "Failure at:", Sinputs[s], "cohorts,",i,"countries",sig,"slope"

                #comm.Send(checkerpiece,dest=0)


comm.Gatherv(checkerpiece,[checker,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.BOOL])

comm.Gatherv(timerpiece,[timer,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.DOUBLE])



'''
for sig in xrange(Sigstart,SigUB):
    for i in xrange(Istart,IUB):
        for s in xrange(Sstart,SUB):
            
            try:
                mn.TheWholeSmack(s,i,sig)
                checker[sig-Sigstart,i-Istart,s-Sstart]=True

            except:
                checker[sig-Sigstart,i-Istart,s-Sstart]=False
'''
               
                
#print checker[:,:,:]
if rank==0:
    checker=np.reshape(checker,(SigUB-Sigstart,IUB-Istart,S))
    for z in xrange(SigUB-Sigstart):
        label="Curve_" +str(z)+".csv"
        np.savetxt(label,checker[z,:,:],delimiter=",")
if rank==0:
    timer=np.reshape(timer,(SigUB-Sigstart,IUB-Istart,S))
    for z in xrange(SigUB-Sigstart):
        label="Time_"+str(z)+".csv"
        np.savetext(label,timer[z,:,:],delimiter=",")

