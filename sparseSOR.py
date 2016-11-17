# -*- coding: utf-8 -*-
"""
MIS40530 Programming Assignment 3

nas_Boland_Gorman_ODonoghue_prog2

Student Name: Micheál Boland
Student Number: 15204343
Date: 04/11/2016 

"""

import numpy as np
import math
import sys
import copy

## FUNCTION DEFINITIONS ##############################
# This takes in filename that has matrix stored in required format and stores it in compressed format
def csr_store(file_name):
    val = []
    col = []
    rowStart = []    
    
    with open(file_name,'r') as infile:
        n = int(infile.readline())
        b = np.zeros((n,1))
        colSum = np.zeros(n)
        mat_temp = infile.readlines()
        
        for i in range(n+1):
            if i != n:
                startCheck = 0
                rowSum = 0
                
                for j in range(n):
                    
#                    if (i==j) and float(mat_temp[i].strip().split("\t")[j]) == 0:
#                        diagZero = 0
#                    elif (i==j) and float(mat_temp[i].strip().split("\t")[j]) != 0:
#                        diagZero = 1
#                        diagEntry = float(mat_temp[i].strip().split("\t")[j])
#                        colSum[j] -= diagEntry
                        
                    if float(mat_temp[i].strip().split("\t")[j]) != 0:
                        
                        startCheck += 1
                        #colSum[j] += float(mat_temp[i].strip().split("\t")[j])
                        #rowSum += float(mat_temp[i].strip().split("\t")[j])
                        
                        val.append(float(mat_temp[i].strip().split("\t")[j]))
                        col.append(j)
                        
                        if startCheck == 1:
                            rowStart.append(len(val)-1)
                
                #print(i , ":",rowSum)
                #print(i, ":", diagEntry - (rowSum - diagEntry))
                
            else:
                # b read in as a column vector            
                b[0:,0] = mat_temp[i].strip().split("\t")
        
        #print(colSum)
        rowStart.append(len(val))
        return {'val':val, 'col':col, 'rowStart':rowStart, 'b':b}  
        
## 
    
def vect_norm(x,y):
    if len(x) == len(y):
        norm = 0
        
        for i in range(len(x)):
            norm += (x[i]-y[i])**2
            
        return float(math.sqrt(norm))
        
    else:
        print("ERROR: length of x and y are not equal")
        
def sparse_sor(A,b,maxits,e,w):
    
    x = np.zeros(len(A['rowStart'])-1)
    x_old = []
    k = 0
    vect_norm_old = 0
        
    while k < maxits:
        
        # This is the check for zero on diagonal. This sum should equal number of rows. If less than there is a zero entry on diagonal
        diagCount = 0                
        #Equating two lists in Python is pass by value, not by reference! Need to use copy         
        x_old = copy.copy(x)
        
        for i in range(len(A['rowStart'])-1):
                        
            i_sum = 0
            for j in range(A['rowStart'][i],A['rowStart'][i+1]):
                
                i_sum += A['val'][j] * x[A['col'][j]]
               
                if A['col'][j] == i:
                # Getting diagonal entry            
                    d = A['val'][j]
                    # Zero diagaonal check. There should be the same number of diagaonal entries as rows
                    diagCount += 1
        
            x[i] += w * (b[i,0] - i_sum)/d
                
        # Divergence and Convergence check
        if k > 0 and (vect_norm(x,x_old) <= (e + (4.0 * sys.float_info.epsilon * np.linalg.norm(x_old)))):
                  
            print("\nSparseSOR Result: X Convergence")            
            return {"x":x,"stopReason":"x Sequence Convergence","numIts":k+1}
        elif k > 0 and (vect_norm_old < vect_norm(x,x_old)):
        
            if diagCount < len(A['rowStart'])-1:
                print("\nSparseSOR Result: Divergence from zero on Diagonal")
                return {"x":x,"stopReason":"Zero on diagonal","numIts":k+1}
            else:
                print("\nSparseSOR Result: X Divergence")
                return {"x":x,"stopReason":"x Sequence Divergence","numIts":k+1}
            
        vect_norm_old = vect_norm(x,x_old)   
        k += 1
    
    print("\nSparseSOR Result: maxIts Reached")
    return {"x":x,"stopReason":"Max Iterations Reached","numIts":k}

## END FUNCTION DEFINITIONS ############################

## This section prompts for user input for desired input or output filename
# User prompted for output filename. If no value entered default nas_Sor.in selected.    
ip_name = input("Please enter full input filename including extension: ")
if ip_name == "":
    ip_name = "nas_Sor.in"
    print("Default ", ip_name, " selected as input filename")
else:
    print(ip_name, " selected as input filename")

# User prompted for output filename. If no value entered default nas_Sor.out selected.    
op_name = input("Please enter full output filename including extension: ")
if op_name == "":
    op_name = "nas_Sor.out"
    print("Default ",op_name, " selected as output filename")
else:
    print(op_name, " selected as output filename")
##

A = csr_store1(ip_name)
print("\nINPUT MATRICES:")
print("A in CSR format")
print("val:",A['val'])
print("col:",A['col'])
print("rowStart:",A['rowStart'])
print("b COLUMN MATRIX:")
b = A['b']
print(b)

#SparseSOR input parameters
maxIts = 100
e = 10**-5
result = sparse_sor(A,b,maxIts,e,1.25)

## Writing to nas_Sor.out
headings = ['Stopping reason','Max num of iterations','Number of iterations','Macine epsilon','X seq tolerance']
content = [result['stopReason'],str(maxIts),str(result['numIts']),str(sys.float_info.epsilon),str(e)]

# This section will append traling spaces where necessary to align in output file
for i in range(len(headings)):
    if len(headings[i])>len(content[i]):
        content[i] = content[i].ljust(len(headings[i]))
    elif len(headings[i])<len(content[i]):
        headings[i] = headings[i].ljust(len(content[i]))

# Writing to output file
with open(op_name,'w') as outfile:
    for entry in headings:    
        outfile.write(str(entry)+"\t")
    outfile.write('\n')
    for entry in content:    
        outfile.write(str(entry)+"\t")
    if (result['stopReason'] == "Max Iterations Reached") or (result['stopReason'] == "x Sequence Convergence"):
        outfile.write('\n')
        for entry in result['x']:
            outfile.write(str(entry)+"\t")

print("Output saved to: ", op_name)

######## TRASH #################################

##-------------------
##-------------------
## 4/11/2016 This section has the SOR working for small matrices---------
#
#def SOR(A,b,w,maxits):
#
#    x = np.zeros_like(b)
#    
#    for i in range(maxits):
#        
#        x_k = np.zeros_like(x)
#        
#        for j in range(A.shape[0]): #Returning length of first row of A
#            s1 = np.dot(A[j,:j],x_k[:j])
#            s2 = np.dot(A[j,j + 0:],x[j + 0:])
#            x_k[j] = x[j] + (b[j] - s1 - s2)*(w/A[j,j])
#            
#        x = x_k
#    
#    return x
#
#A = np.matrix([[3.0, -1.0, 1.0],[-1.0, 3.0, -1.0],[1.0, -1.0, 3.0]])
#b = np.array([-1.0,7.0,-7.0])
#w = 1.25
#maxits = 100
#    
##print(SOR(A,b,w,maxits))

#-----------------

# Using example from notes
#B = np.matrix([[21,0,0,12,0,0],
#               [0,1,49,0,0,0],
#                [31,16,65,0,0,23],
#                [0,0,0,85,0,0],
#                [55,0,0,0,42,0],
#                [0,0,0,0,0,41]])
#
#csr = csr_store(B)            
#print('val = ',csr['val'])
#print('col = ',csr['col'])
#print('rowStart = ',csr['rowStart'])
## This of course returning CSR with convention that first element is 0 rather than 1. Python implementation
#def csr_store(B):
#    # Obtain size of matrix. Cycling through each row    
#    val = []
#    col = []
#    rowStart = []    
#    for i in range(B.shape[0]):
#        # Set this to zero at start of each new row. Used to populate rowStart vector        
#        startCheck = 0
#        
#        # Iterating through each column
#        for j in range(B.shape[0]):
#            
#            # Ignore zero entries
#            if B[i,j] != 0:
#                
#                # When first entry in a row reached iterate check value
#                startCheck += 1
#                # Appending relevant values to val and col vectors
#                val.append(B[i,j])
#                col.append(j)
#                
#                # Here if startCheck is 1 then new rowStart value needs to be appended
#                if startCheck == 1:
#                    rowStart.append(len(val)-1)
#                    
#    # Finally append the last value to row start, namely the length of the val vector
#    rowStart.append(len(val))
#    # Output is returned as a dictionary that can be referenced by the various different key values    
#    return{'val':val, 'col':col, 'rowStart':rowStart}


    



