# -*- coding: utf-8 -*-
"""
MIS40530 Programming Assignment 3

nas_Boland_Gorman_ODonoghue_prog2

Student Name: Micheál Boland
Student Number: 15204343
Date: 04/11/2016 

"""

import numpy as np

#### SECTION 1. Code to read in file nas_Sor.in

# Have this to be reading in tab separated files.    
with open('nas_Sor.in','r') as infile:
    n = int(infile.readline())
    A = np.zeros((n,n))
    b = np.zeros(n)
    mat_temp = infile.readlines()
    
    for row in range(n+1):
        if row != n:
            for col in range(n):
                # Now have this reading in one element at a time. Easy to pass into the CSR storage format
                A[row,col] = float(mat_temp[row].strip().split("\t")[col])
        else:
            b[0:] = mat_temp[row].strip().split("\t")
    
    print(A)
    print(b)

#### END SECTION 1

#### SECTION 2. Code to store sparse matix in CSR format

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

## END FUNCTION DEFINITION

# Using example from notes
B = np.matrix([[21,0,0,12,0,0],
               [0,0,49,0,0,0],
                [31,16,0,0,0,23],
                [0,0,0,85,0,0],
                [55,0,0,0,91,0],
                [0,0,0,0,0,41]])

csr = csr_store(B)            
print('val = ',csr['val'])
print('col = ',csr['col'])
print('rowStart = ',csr['rowStart'])

#### END SECTION 2

#### SECTION 3: Check for zero values on diagonal entries

# Have this implemented that we need to get a value of n, for a n x n matrix back in the value diag_check 
# when the col and rowStart are searched based on logic below
diag_check = 0

# The length of (rowStart - 1) is n for the n x n matirx
for i in range(len(csr['rowStart'])-1):
    
    # Interval iteration between the i-th and (i+1)-th entries of rowStart
    for j in range(csr['rowStart'][i],csr['rowStart'][i+1]):
        # If the column entry is equal to the row entry then there is a non-zero diagonal in the matrix
        if csr['col'][j] == i:
            # Check gets incremented            
            diag_check += 1

# FInal check to see that the expected number of non-zero diagonal entries were found            
if diag_check == len(csr['rowStart']) - 1:
    print("Non-zero values present in all diagonal entries")
else:
    print('Error: Zero entries found on diagonal')    

#### END SECTION 3

#### SECTION 4: Function to calculate vector norms

# Just need a function to calculate vector norm at each iteration of x

#### END SECTION 4

#### SECTION 5: Function to check stopping reason

# In here need check for divergence at each step or if maxits reached
# Write as a function that can be called after every iteration

#### END SECTION 5

#### SECTION 6 Sparse SOR Implementation

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
   

c = np.array([1,2,3,4,5,6])
print(sparse_sor(csr,c,5,0.1,1.25))


#### END SECTION 6

### SECTION 7: Writing to nas_Sor.out

# Need to have option to have custom filename otherwise nas_Sor.out
# Need to work on formatting here to have the entries aligned
with open('nas_Sor.out','w') as outfile:
    outfile.write('Stopping reason\tMax num of iterations\tNumber of iterations\tMacine epsilon\tX seq tolerance\tResidual seq tolerance\n')
    outfile.write('TestA\tTestB\tTestC\tTestD\tTestE\tTestF\n')

### END SECTION 7

#-------------------
#-------------------
# 4/11/2016 This section has the SOR working for small matrices---------

def SOR(A,b,w,maxits):

    x = np.zeros_like(b)
    
    for i in range(maxits):
        
        x_k = np.zeros_like(x)
        
        for j in range(A.shape[0]): #Returning length of first row of A
            s1 = np.dot(A[j,:j],x_k[:j])
            s2 = np.dot(A[j,j + 0:],x[j + 0:])
            x_k[j] = x[j] + (b[j] - s1 - s2)*(w/A[j,j])
            
        x = x_k
    
    return x

A = np.matrix([[3.0, -1.0, 1.0],[-1.0, 3.0, -1.0],[1.0, -1.0, 3.0]])
b = np.array([-1.0,7.0,-7.0])
w = 1.25
maxits = 100
    
#print(SOR(A,b,w,maxits))

#-----------------

    



