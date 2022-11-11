#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np
from numpy.linalg import norm, qr
import copy

def power(in_mat, size, maxIter, errLim):
    n = size
    mat_a = np.matrix(in_mat, dtype = float)
    x_vec = np.zeros(shape=(size,1), dtype = float)
    x_vec[0,0] = 1.0
    for no in range(maxIter):
        temp = mat_a*x_vec
        e1 = norm(temp)
        new_x = temp/e1
        e2 = norm(mat_a*new_x)
        err = abs(e2-e1)*100/e2
        x_vec = new_x
        if err < errLim:
            break
    
    return x_vec.T.tolist()[0], e2, no

def inversepower(in_mat, size, maxIter, errLim):
    mat_a = np.matrix(in_mat, dtype = float).I
    return power(mat_a, size, maxIter, errLim)

def inversePowerShift(in_mat, size, maxIter, errLim):
    epsE = float(input("Find Eigenvalue closest to:"))
    mat_a = np.matrix(in_mat, dtype = float)
    mat_a = mat_a - epsE*np.identity(size)
    mat_a = mat_a.I
    x_vec, e, no = power(mat_a, size, maxIter, errLim)
    return x_vec, (1/e)+epsE, no

def qrM(in_mat, size, maxIter, errLim):
    mat_a = np.matrix(in_mat, dtype = float)
    eVals = [mat_a[x,x] for x in range(size)]
    for no in range(maxIter):
        q, r = qr(mat_a)
        new_a = r*q
        new_eVals = [new_a[x,x] for x in range(size)]
        errLis = [(new_eVals[x] - eVals[x])*100/eVals[x] for x in range(size)]
        err = max(errLis)
        eVals = new_eVals
        mat_a = new_a
        if err < errLim:
            break
    return eVals,no+1

def printVec(vec):
    for x in vec:
        print(round(x,4))

#with open("input.txt", 'r') as f:
#    n = int(f.readline().strip())
#    mat_a = []
#    for x in range(n):
#        row = map(float, f.readline().strip().split(' '))
#        mat_a.append(row)
#    maxIter = int(f.readline().strip())
#    errLim = float(f.readline().strip())
#    epsE = float(f.readline().strip())
    
n = int(input("Enter the number of rows in the matrix:"))
print("Enter the matrix:")
mat_a=[list(map(float, input().split())) for x in range(n)]
maxIter = int(input("Maximum iterations:"))
errLim = float(input("Maximum relative approximate error (in percentage):")) / 100

display_str = """Enter the number of method you want to use:
1: Power Method
2: Inverse Power Method
3: Inverse Power Method with shift 
4: QR method
Choose option 1-4: """
choice = int(input(display_str).strip())
if choice == 1:
    x_vec, eVal, no = power(mat_a, n, maxIter, errLim)
    print("\nEigenvalue")
    print (round(eVal,4))
    print ("\nEigenvector")
    printVec(x_vec)
    print ("\nIterations")
    print (no)
elif choice == 2:
    x_vec, eVal, no = inversepower(mat_a, n, maxIter, errLim)
    print ("\nEigenvalue")
    print (round(1/eVal,4))
    print ("\nEigenvector")
    printVec(x_vec)
    print ("\nIterations")
    print (no)
elif choice == 3:
    x_vec, eVal, no = inversePowerShift(mat_a, n, maxIter, errLim)
    print ("\nEigenvalue")
    print (round(eVal,4))
    print ("\nEigenvector")
    printVec(x_vec)
    print ("\nIterations")
    print (no)
elif choice == 4:
    eVal, no = qrM(mat_a, n, maxIter, errLim)
    print ("\nEigenvalues")
    printVec(eVal)
    print ("\nIterations")
    print (no)
else:
    print("\nInvalid input")

