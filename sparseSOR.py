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
from sparseSORModule import vect_norm, csr_store, sparse_sor
import matplotlib.pyplot as plt


### FUNCTION DEFINITIONS ##############################

def fname_prompt():
    ## This function prompts for user input for desired input or output filename
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
        
    return {'ip_name':ip_name,'op_name':op_name}

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
                        
                    if float(mat_temp[i].strip().split("\t")[j]) != 0:
                        
                        startCheck += 1
                        val.append(float(mat_temp[i].strip().split("\t")[j]))
                        col.append(j)
                        
                        if startCheck == 1:
                            rowStart.append(len(val)-1)
                
            else:
                # b read in as a column vector            
                b[0:,0] = mat_temp[i].strip().split("\t")
        
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

# Prompting user for filename input
fnames = fname_prompt()

# Storage of A in CSR format
A = csr_store(fnames['ip_name'])
print("\nINPUT MATRICES:")
print("A in CSR format")
print("val:",A['val'])
print("col:",A['col'])
print("rowStart:",A['rowStart'])

print("b COLUMN MATRIX:")
b = A['b']
print(b)

#SparseSOR input parameters
maxIts = 10000
e = 10**-5
w = 1.25
result = sparse_sor(A,b,maxIts,e,w)

# Creating lists for headers and content for .out file
headings = ['Stopping reason','Max num of iterations','Number of iterations','Machine epsilon','X seq tolerance']
content = [result['stopReason'],str(maxIts),str(result['numIts']),str(sys.float_info.epsilon),str(e)]

# This section will append traling spaces where necessary to align headers with content in output file
for i in range(len(headings)):
    if len(headings[i])>len(content[i]):
        content[i] = content[i].ljust(len(headings[i]))
    elif len(headings[i])<len(content[i]):
        headings[i] = headings[i].ljust(len(content[i]))

# Writing to output file
with open(fnames['op_name'],'w') as outfile:
    for entry in headings:    
        outfile.write(str(entry)+"\t")
    outfile.write('\n')
    for entry in content:    
        outfile.write(str(entry)+"\t")
    if (result['stopReason'] == "Max Iterations Reached") or (result['stopReason'] == "x Sequence Convergence"):
        outfile.write('\n')
        for entry in result['x']:
            outfile.write(str(entry)+"\t")

print("Output saved to: ", fnames['op_name'])

# Plotting tolerance versus number of iterations. Uncomment when required
#maxIts = 10000
#e = [10**-1,10**-2,10**-3,10**-4,10**-5,10**-6,10**-7,10**-8]
#
## Collecting values for iterations for defined tolerances
#for i in range(len(e)):
#    iterations.append(sparse_sor(A,b,maxIts,e[i],1.25)['numIts'])
#
## Plotting of result
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.set_ylim(-0.01,0.11)
#plt.plot(iterations,e,'r-')
#plt.plot(iterations,e,'ro')
#for i,j in zip(iterations,e):
#    ax.annotate(("numIts="+str(i)+"\n"+"Tolerance="+str(j)),xy=(i,j))
#plt.xlabel("Tolerance")
#plt.ylabel("Number of Iterations")
#plt.title("Tolerance vs Number of Iterations")
#plt.show()


    



