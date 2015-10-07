import numpy as np
import Main as mn
import sys
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

comm=MPI.COMM_WORLD
rank=comm.Get_rank()
size=comm.Get_size()

Sinputs=np.array([20,40,60,80]) #Input the levels of numbers of cohorts you want checked
S=Sinputs.size
IUB=8 #Upper bound for the number of countries
Istart=2 #Lower bound for the number of countries, absolute minimum is 2.
SigUB=4 #Upper bound for the sigma parameter (The slope of the curve)
Sigstart=1 #Lower Bound for the sigma parameter
IRange=IUB-Istart #Range of the number of countries
SigRange=SigUB-Sigstart #Range of the slope parameter

#REMEMBER THAT THE NUMBER OF PROCESSORS-1 MUST DIVIDE SIG_RANGE*IRANGE*S.
#THE FIRST PROCESSOR IS USED FOR CSV GENERATION.

checker=np.zeros(((SigRange)*(IRange)*(S)),dtype=bool) if rank==0 else None

if rank!=0:
    checkerpiece=np.zeros(((IRange)*(S)),dtype=bool) 
else: 
    checkerpiece=np.zeros(0,dtype=bool)

def Send(Sigrange,IRange):
    '''
    Description: This automatically creates the tuple for that will tell
    MPI where to send each of the pieces of checkerpiece.

    Inputs:
        -Sigrange: Range of the sigma values
        -IRange: Range of the number of countries

    Outputs:
        -Tuple1

    '''
    tuple1 = (0,)
    for i in xrange(Sigrange):
        tuple1+= (IRange*S,)
    
    return tuple1

def Recieve(Sigrange,IRange):
    '''
    Description: Similar to the Send, it automatically creates the tuple
    that is needed to send

    Inputs:
        -Sigrange:
        -IRange

    Outputs:
        -tuple2

    '''
    tuple2 = (0,0,)
    for i in xrange(1,Sigrange):
        tuple2+=(IRange*S*i,)

    return tuple2

tuple1=Send(SigRange,IRange)
tuple2=Recieve(SigRange,IRange)

comm.Scatterv([checker,tuple1,tuple2,MPI.BOOL],checkerpiece)

if rank!=0:
    sig=rank
    for s in xrange(S):
        for i in xrange(Istart,IUB):
            try:
                mn.Multi_Country(Sinputs[s],i,sig)
                checkerpiece=np.reshape(checkerpiece,(IRange,S))
                checkerpiece[i-Istart,s]=True
                checkerpiece=np.reshape(checkerpiece, (IRange)*(S))
                print "Success at:", Sinputs[s], "cohorts,",i,"countries",sig,"slope"
            except:
                checkerpiece=np.reshape(checkerpiece[:,np.newaxis],(IRange,S))
                checkerpiece[i-Istart,s]=False
                checkerpiece=np.reshape(checkerpiece, (IRange)*(S))
                print "Failure at:", Sinputs[s], "cohorts,",i,"countries",sig,"slope"



comm.Gatherv(checkerpiece,[checker,tuple1,tuple2,MPI.BOOL])


#This creates the csv sheets, one for each slope parameter, that shows the results of
#the robustness check. Each cell has a "True" where the model worked, and a "False"
#Where it did not.
if rank==0:
    checker=np.reshape(checker,(SigRange,IRange,S))
    for z in xrange(SigRange):
        label="Curve_" +str(z+1)+".csv"
        np.savetxt(label,checker[z,:,:],delimiter=",")

