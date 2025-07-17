# Read and plot two dimensional data written in the .dat format
#

import sys
import numpy as np
import matplotlib.pyplot as plt

def readfortranvec(n,infile) :
    vecout=[]
    i=0
    while i<n :
        nextline=list(map(float,next(infile).split()))
        i+=len(nextline)
        vecout+=nextline
    return vecout

FileName=sys.argv[1]
print("Read File",FileName)
infile=open(FileName,"r")
for LabelIn in infile :
# LabelIn=next(infile)
    flag=int(next(infile))
    assert flag == -2, "Not a 2D array"
    N = int(next(infile))
    M = int(next(infile))
    XVal = readfortranvec(N,infile)
    YVal = readfortranvec(M,infile)
    ZValLin = readfortranvec(N*M,infile)
    ZVal=np.reshape(ZValLin,(N,M),order='F')
    plt.plot()
    if LabelIn.split()[0] == "TDL" :
        curves=plt.contour(XVal,YVal,ZVal.transpose(),np.linspace(0,500,11))
        plt.clabel(curves,fmt="%d")
    else :
        curves=plt.contour(XVal,YVal,ZVal.transpose(),11)
        plt.clabel(curves)
    plt.title(LabelIn[0:39]+"\n"+LabelIn[40:],y=.98)
# plt.savefig("testout.pdf")
    plt.show()