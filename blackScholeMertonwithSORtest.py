#import Libraries
from math import *
import numpy
import numpy as np
import matplotlib . pyplot as plt
from decimal import *
import random
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import time
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

#set r and sigma as suggested
r=0.02
sigma=0.3
#T is maturity date
T=100
#strike price is X
X=50
#no timestamps
M=100
#used in equation
k=T/M
# create zero matrix nxn
#Size of matrix

n=150
#initial matrix for data colection during iterations
outputDataSparse= numpy.zeros(shape=(T,n))
#Initialise A Matrix
A= numpy.zeros(shape=(n,n))
#create zero b matrix nx1
coefficientsSparse=numpy.zeros(shape=(n,1))



## FUNCTION DEFINITION
# This of course returning CSR with convention that first element is 0 rather than 1. Python implementation
def csr_store(B):
    # Obtain size of matrix. Cycling through each row    
    val = []
    col = []
    rowStart = []    
    for i in range(B.shape[0]):
        # Set this to zero at start of each new row. Used to populate rowStart vector        
        startCheck = 0
        
        # Iterating through each column
        for j in range(B.shape[0]):
            
            # Ignore zero entries
            if B[i,j] != 0:
                
                # When first entry in a row reached iterate check value
                startCheck += 1
                # Appending relevant values to val and col vectors
                val.append(B[i,j])
                col.append(j)
                
                # Here if startCheck is 1 then new rowStart value needs to be appended
                if startCheck == 1:
                    rowStart.append(len(val)-1)
                    
    # Finally append the last value to row start, namely the length of the val vector
    rowStart.append(len(val))
    # Output is returned as a dictionary that can be referenced by the various different key values    
    return{'val':val, 'col':col, 'rowStart':rowStart}


# For the moment just have maxits check. Integrate section 5 once complete
def sparse_sor(A,b,maxits,e,w):
    
    x = np.zeros(len(A['rowStart'])-1)
    
    for k in range(maxits):
        
        for i in range(len(A['rowStart'])-1):
            
            i_sum = 0
            for j in range(csr['rowStart'][i],csr['rowStart'][i+1]):
                
                i_sum += A['val'][j] * x[A['col'][j]]
               
                if csr['col'][j] == i:
                # Getting diagonal entry            
                    d = A['val'][j]
        
            x[i] += w * (b[i] - i_sum)/d
                        
    return x               
   


#fn,m = f(nh,mk).
#n is intervals 0 - Smax
#S0 = 0, Sn = nh = Sn−1 + h, 1  n  N − 1, SN = Smax
#t0 = 0, tm = mk = tm−1 + k, m = 0, 1, . . . ,M.
#fn,M = max{X − nh, 0}, n = 0, 1, . . . ,N
#The value of the put option when S = 0 is just the strike price X, giving
#f0,m = X, m = 0, 1, . . . ,M.
#As assumed earlier, the value of the put option is zero when S = Smax, giving
#fN,m = 0, m = 0, 1, . . . ,M.
#These define the value of the option along three edges of our grid (where S = 0, S = Smax and
#t = T) and allow us to get started.

                                                                     
#set the size of the results matrix and fill with zeros
#this is the number from M-1 down to 0 of iterarions of solve
#first build the b part of the solution
#code is the same as fnM=max{S-nk}
#just use increments initially and look at 50 to 100 strike price in increments of 1
for i in range(0,n):
       
       if (i==0):
           #b[i]=50+((k/2)*((sigma*sigma)-r)*X)
           coefficientsSparse[i]=max(X-(i*k),0)
           coefficientsSparse[i]=coefficientsSparse[i]+((k/2)*((sigma*sigma)-r)*X)
       else:
           #b[i]=b[i-1]+2
           coefficientsSparse[i]=max(X-(i*k),0)
start=coefficientsSparse
start1=coefficientsSparse[0]
stop1=coefficientsSparse[n-1]
#print(start)

#first build initial version of the regression coefficients 
#perform lin regression for the range of time samples and use these for next ieration    
t1 = time.time()
for p in range(0,T,1):
    outputDataSparse[p,0:n]=coefficientsSparse.T
    #loop for M-1 to 0
    # run the sparse sor
    #only create the matrix A on the first iteration
    if (p==0):
        for i in range(0, n):
              # loop for BSM coefficient creation
              for j in range(0, n):
                  #all other matrix solutions except top right and bottom left
                  if (i==j) and (j!=0) and (j!=n-1):
                      A[i,j]=(1+(k*r)+(k*(sigma*sigma)*((i+1)*(i+1))))
                      A[i-1,j]=((-i*k)/2)*(((i)*sigma*sigma)+r)
                      A[i+1,j]=((-(i+2)*k)/2)*(((i+2)*sigma*sigma)-r)
                      
                  #top left unique set up
                  if (i==0) and (j==0):
                      A[i,j]=(1+(k*r)+(k*(sigma*sigma)*((i+1)*(i+1))))
                      A[i+1,j]=((-(i+2)*k)/2)*(((i+2)*sigma*sigma)-r)
                      
                  #bottom right unique set up    
                  if (i==n-1) and (j==n-1):
                      A[i,j]=(1+(k*r)+(k*(sigma*sigma)*((i+1)*(i+1))))
                      A[i-1,j]=((-i*k)/2)*(((i)*sigma*sigma)+r)
                  
    # run the sparse sor
    #only perform tghe csr on the first iteration
    if (p==0):
        csr = csr_store(A)
    #time the solve
    t1 = time.time()
    #solve the new co-efficients andset the extremities
    coefficientsSparse=sparse_sor(csr,coefficientsSparse,100,0.1,1.25)
    coefficientsSparse[0]=start1
    coefficientsSparse[n-1]=stop1
    t2 = time.time()
    
    #print(t2-t1)

outputDataNew=numpy.zeros(shape=(T,n))    
for p in range(0,M,1):
    #outputDataNew[0:len(outputDataSparse),p]=outputDataSparse[0:len(outputDataSparse),n-p-1]
    outputDataNew[p,0:n]=outputDataSparse[M-p-1,0:n]

# now plot 2D results for the BSM
fig = plt.figure()
fig.suptitle('Implemented BSM', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
plt.plot(outputDataNew)
ax.set_xlabel('number of iterations M-1 to 0')
ax.set_ylabel('Stock Price BSM value using SOR')
plt.show()


#plot the 3D results for BSM for the parameters
fig = plt.figure()
fig.suptitle('Implemented BSM', fontsize=14, fontweight='bold')
ax = fig.gca(projection='3d')
x = np.arange(0, n, 1)
y = np.arange(0, M, 1)
X, Y = np.meshgrid(x, y)
Z=outputDataSparse
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_zlim(0, 50)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_ylabel('Days to expire')
ax.set_xlabel('Stock Price')
ax.set_zlabel('Options')
plt.show()
