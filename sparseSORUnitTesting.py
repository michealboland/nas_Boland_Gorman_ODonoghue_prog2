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

def fname_prompt():
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
        
    return {'ip_name':ip_name,'op_name':op_name}


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

# Input and Output File Name Handling Test Cases

print("---INPUT AND OUTPUT FILE NAMING TESTING---")
print("TC 1.1: USER DEFINED INPUT FOR BOTH INPUT AND OUTPUT FILENAMES")
print("Please enter a user defined filename for both the input and output filename")

fnames = fname_prompt()

print("TC 1.1 RESULT")
print("Value of input filename: ",fnames['ip_name'])
print("Value of output filename: ",fnames['op_name'])

print("TC 1.2: USER DEFINED INPUT FOR ONLY INPUT FILENAME")
print("Please enter a user defined filename for only the input filename. Please hit enter when prompted for output filename")

fnames = fname_prompt()

print("TC 1.2 RESULT")
print("Value of input filename: ",fnames['ip_name'])
print("Value of output filename: ",fnames['op_name'])

print("TC 1.3: USER DEFINED INPUT FOR ONLY OUTPUT FILENAME")
print("Please enter a user defined filename for only the output filename. Please hit enter when prompted for input filename")

fnames = fname_prompt()

print("TC 1.3 RESULT")
print("Value of input filename: ",fnames['ip_name'])
print("Value of output filename: ",fnames['op_name'])

print("TC 1.4: USER DEFINED INPUT FOR NEITHER INPUT NOR OUTPUT FILENAMES")
print("Please hit enter when prompted for both the input and output filenames")

fnames = fname_prompt()

print("TC 1.4 RESULT")
print("Value of input filename: ",fnames['ip_name'])
print("Value of output filename: ",fnames['op_name'])

# Vector Norm Calculation Test Cases

print("TC 2.1:  CALCULATION OF VECTOR NORM OF TWO ROW VECTORS OF EQUAL LENGTH")
print("The vectors norm of the below vectors are to be calculated")
x = np.array([1,2,3,4,5])
y = np.array([9,8,7,6,5])
print("Vector 1: ",x)
print("Vector 2: ",y)
print("TC 2.1 RESULT")
print("Vector Norm: ", vect_norm(x,y))


print("TC 2.2:  CALCULATION OF VECTOR NORM OF TWO COLUMN VECTORS OF EQUAL LENGTH")
print("The vectors norm of the below vectors are to be calculated")
x.shape = (5,1)
y.shape = (5,1)
print("Vector 1: \n",x)
print("Vector 2: \n",y)
print("TC 2.2 RESULT")
print("Vector Norm: ", vect_norm(x,y))


print("TC 2.3:  CALCULATION OF VECTOR NORM OF A COLUMN AND ROW VECTOR OF EQUAL LENGTH")
print("The vectors norm of the below vectors are to be calculated")
x.shape = (1,5)
print("Vector 1: \n",x)
print("Vector 2: ",y)
print("TC 2.3 RESULT")
print("Vector Norm: ", vect_norm(x,y))


print("TC 2.4:  CALCULATION OF VECTOR NORM OF TWO ROW VECTORS OF UNEQUAL LENGTH")
print("The vectors norm of the below vectors are to be calculated")
x = np.array([1,2,3,4,5])
y = np.array([9,8,7,6])
print("Vector 1: ",x)
print("Vector 2: ",y)
print("TC 2.4 RESULT")
print("Vector Norm: ", vect_norm(x,y))

# Reading content from file and storing directly in CSR format

print("TC 3.1:  TAKING MATRIX A AND b FROM FILE AND STORING A IN CSR FORMAT")
print("The below matrices are stored in the file type1.in")
print("A:")
print("[225.1270 81.5323  16.1611 34.0935 60.1489")
print("44.0649	 316.7961	37.0439 16.5390 93.6523")
print("99.9771	 30.8879	0	    11.7380 37.3152")
print("39.5335	 20.6465	75.6235	318.1132 38.4645")
print("70.6488	 7.7633	   94.5225	60.3797	  383.8700]")
print("b:")
print("[1")	
print("2")	
print("3")	
print("4")	
print("5]")
print("\nPlease enter type1.in when prompted for input file name")

fnames = fname_prompt()
test = csr_store(fnames['ip_name'])

print("TC 3.1 RESULT")
print("\nMatrix A in CSR format and b column vectors as read from file")
print("val: ",test['val'])
print("col: ",test['col'])
print("rowStart: ",test['rowStart'])
print("b: \n",test['b'])