from math import *
import numpy
import numpy as np
import matplotlib . pyplot as plt
from decimal import *
import random
from sklearn import linear_model
from mpl_toolkits.mplot3d import Axes3D
from scipy import ndimage
import time


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

                                                                     
# create zero matrix nxn
#small matrix for debugging
n=50
#large matrix
outputData= numpy.zeros(shape=(T,n))
A= numpy.zeros(shape=(n,n))
#create zero b matrix nx1
b=numpy.zeros(shape=(n,1))

#set the size of the results matrix and fill with zeros
#this is the number from M-1 down to 0 of iterarions of solve

#set r and sigma as suggested
r=0.025
sigma=0.035
#T is maturity date
T=90
#strike price is X
X=100
#no the timestamps
M=90
#used in equation
k=T/M


#first build the b part of the solution
#code is the same as fnM=max{S-nk}
#just use increments initially and look at 50 to 100 strike price in increments of 1
for i in range(0,n):
       if (i==0):
           #b[i]=50+((k/2)*((sigma*sigma)-r)*X)
           b[i]=max(X-(i*k),0)
           b[i]=b[i]+((k/2)*((sigma*sigma)-r)*X)
       else:
           #b[i]=b[i-1]+2
           b[i]=max(X-(i*k),0)
start=b
#check b values
#print(b)

#first build initial version of the regression coefficients 
#perform lin regression for the range of time samples and use these for next ieration    
#t1 = time.time()
for p in range(0,T,1):
    outputData[p,0:n]=b.T
    #loop for M-1 to 0
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
                  
    #use linear regression until SOR is completed
    #print(A)
    #print(b)
    t1 = time.time()
    coefficients=np.linalg.lstsq(A,b)
    t2 = time.time()
    coefficients=np.asarray(coefficients)
    #pick the new coefficients for next iteration of the loop
    b=coefficients[0]
    #print(b)
    

    print(t2-t1)

outputDataNew=numpy.zeros(shape=(T,n))    
for p in range(0,M,1):
    #outputDataNew[0:len(outputDataSparse),p]=outputDataSparse[0:len(outputDataSparse),n-p-1]
    outputDataNew[p,0:n]=outputData[M-p-1,0:n]

# now plot the results for the BSM
fig = plt.figure()
fig.suptitle('Implemented BSM', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
plt.plot(outputDataNew)
ax.set_xlabel('number of iterations M-1 to 0')
ax.set_ylabel('Stock Price BSM value using SOR')
plt.show()



fig = plt.figure()
ax = fig.gca(projection='3d')
x = np.arange(0, n, 1)
y = np.arange(0, M, 1)
X, Y = np.meshgrid(x, y)
Z=outputDataNew
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_zlim(0, 100)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_ylabel('number of iterations M-1 to 0')
ax.set_xlabel('Stock Price BSM value using LG')
ax.set_zlabel('Z axis')
plt.show()
