import numpy as np
import Main as mn
import sys
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

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

if rank!=0:
    checkerpiece=np.zeros(((IUB-Istart)*(S)),dtype=bool) 
else: 
    checkerpiece=np.zeros(0,dtype=bool)

comm.Scatterv([checker,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.BOOL],checkerpiece)

if rank!=0:
    sig=rank
    for s in xrange(S):
        for i in xrange(Istart,IUB):
            try:
                mn.TheWholeSmack(Sinputs[s],i,sig)
                checkerpiece=np.reshape(checkerpiece,(IUB-Istart,S))
                checkerpiece[i-Istart,s]=True
                checkerpiece=np.reshape(checkerpiece, (IUB-Istart)*(S))
                print "Success at:", Sinputs[s], "cohorts,",i,"countries",sig,"slope"
                #comm.Send(checkerpiece,dest=0)
            except:
                checkerpiece=np.reshape(checkerpiece,(IUB-Istart,S))
                checkerpiece[i-Istart,s]=False
                checkerpiece=np.reshape(checkerpiece, (IUB-Istart)*(S))
                print "Failure at:", Sinputs[s], "cohorts,",i,"countries",sig,"slope"

                #comm.Send(checkerpiece,dest=0)


comm.Gatherv(checkerpiece,[checker,(0,(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S),(IUB-Istart)*(S)),\
        (0,0,(IUB-Istart)*(S),2*(IUB-Istart)*(S),3*(IUB-Istart)*(S),4*(IUB-Istart)*(S)),MPI.BOOL])


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

