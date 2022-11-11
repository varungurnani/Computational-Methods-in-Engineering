#!/usr/bin/env python
# coding: utf-8

# In[10]:


from copy import copy, deepcopy
import numpy
from math import *

def printmatrix(a):
    n=len(a)
    m=len(a[0])
    
    for i in range(n):
        for j in range(m):
            print(round(a[i][j],3),end="  ")
        print("\n")
        

def printvector(x):
    for i in range(len(x)):
        print(round(x[i],3))
    
    
def GE_without_pivoting(mat,n):
    #Forward elimination
    for i in range(0,n-1):
        for j in range(i+1,n):
            alpha=-(mat[j][i])/mat[i][i]
            mat[j][i]=0
            for k in range(i+1,n+1):
                mat[j][k]=mat[j][k]+alpha*mat[i][k]

    #Back substitution
    x=[0]*n
    for i in range(n-1,-1,-1):
        if i==n-1:
            x[n-1]=mat[n-1][n]/mat[n-1][n-1]
            continue
        sumi=0
        for j in range(1+i,n):
            sumi=sumi+mat[i][j]*x[j]
        x[i]=(mat[i][n]-sumi)/mat[i][i]

    print("\nX")
    printvector(x)

def GE_with_pivoting(mat,n):
    I=numpy.eye(n)
    for i in range(0,n-1):
        count=i
        for j in range(i+1,n):
            if mat[count][i]<mat[j][i]:

                count=j
        newpivot=mat[count]
        newpivoti=deepcopy(I[count])
        mat[count]=mat[i]
        I[count]=I[i]
        mat[i]=newpivot
        I[i]=newpivoti

    GE_without_pivoting(mat,n)

def LU_decomposition_by_using_doolittle_without_pivoting(U,n):
    L= [[0 for x in range(n)] for y in range(n)]
    for i in range(n):
        L[i][i]=1
    for i in range(0,n-1):
        for j in range(i+1,n):
            alpha=-(U[j][i])/U[i][i]
            L[j][i]=-1*alpha
            U[j][i]=0
            for k in range(i+1,n):
                U[j][k]=U[j][k]+alpha*U[i][k]

    #Forward substitution
    y=[0]*n
    for i in range(0,n):
        if i==0:
            y[i]=U[0][n]
            continue
        sumi=0
        for j in range(0,i):
            sumi=sumi+L[i][j]*y[j]
        y[i]=(U[i][n]-sumi)

    #Back substitution
    x=[0]*n
    for i in range(n-1,-1,-1):
        if i==n-1:
            x[n-1]=y[n-1]/U[n-1][n-1]
            continue
        sumi=0
        for j in range(1+i,n):
            sumi=sumi+U[i][j]*x[j]
        x[i]=(y[i]-sumi)/U[i][i]
    print('\nX')
    printvector(x)
    print('\nL')
    printmatrix(L)
    
    for l in range(len(U)):
        U[l].pop(-1)
    print('\nU')
    printmatrix(U)


def LU_decomposition_by_using_crout_without_pivoting(U,n):
    L= [[0 for x in range(n)] for y in range(n)]
    for i in range(n):
        L[i][i]=U[i][i]
        for t in range(0,n):
            U[i][t]=U[i][t]/L[i][i]
        for j in range(i+1,n):
            L[j][i]=U[j][i]
            U[j][i]=0
            for k in range(i+1,n):
                U[j][k]=U[j][k]-L[j][i]*U[i][k]
    L1=deepcopy(L)
    U1=deepcopy(U)
    

    #Forward substitution
    y=[0]*n
    for i in range(0,n):
        if i==0:
            y[i]=U[0][n]/L[0][0]
            continue
        sumi=0
        for j in range(0,i):
            sumi=sumi+L[i][j]*y[j]
        y[i]=(U[i][n]-sumi)/L[i][i]

    #Back substitution
    x=[0]*n
    for i in range(n-1,-1,-1):
        if i==n-1:
            x[n-1]=y[n-1]
            continue
        sumi=0
        for j in range(1+i,n):
            sumi=sumi+U[i][j]*x[j]
        x[i]=(y[i]-sumi)
    print("\nX")
    printvector(x)
    print("\nL")
    printmatrix(L1)
    for l in range(len(U1)):
        U1[l].pop(-1)
    print("\nU")
    printmatrix(U1)

def Cholesky_decomposition(U,n):
    
    l=numpy.zeros(shape=(n,n),dtype=float)

    for i in range(n):
        for j in range(n):
            sum = 0
 
            if (j == i):
                for k in range(j):
                    sum += (l[j][k])**2
                l[j][j] = numpy.sqrt(U[j][j]-sum)
            else:
                for k in range(j):
                    sum += (l[i][k] * l[j][k])
                if l[j][j]==0:
                    l[i][j] = (U[i][j] - sum)
                else:
                    l[i][j] = (U[i][j] - sum)/l[j][j]
    for i in range(n):
        for j in range(i+1,n):
            l[i][j]=0

    #solving for x vector
    b=U*numpy.matrix([[0],[0],[0],[1]])
    y=(numpy.linalg.inv(l))*b

    x=(numpy.linalg.inv(l.transpose()))*y
    
    
    print("\nX")
    for i in range(len(x)):
        print(round(float(x[i]),2))
    
    print("\nL")
    printmatrix(l)


n=int(input('Enter the number of equations:'))
print('Enter the augmented matrix:')
a=[list(map(float, input().split())) for x in range(n)]
g=input('Please select your method:\nA. Gauss Elimination method without pivoting\nB. Gauss Elimination method with pivoting\nC. LU decomposition by using Doolittle method without pivoting\nD. LU decomposition by using Crout method without pivoting\nE. Cholesky decomposition for symmetric positive definite matrix\n')


if g=='A':
    GE_without_pivoting(a,n)
elif g=='B':
    GE_with_pivoting(a,n)
elif g=='C':
    LU_decomposition_by_using_doolittle_without_pivoting(a,n)
elif g=='D':
    LU_decomposition_by_using_crout_without_pivoting(a,n)
elif g=='E':
    Cholesky_decomposition(a,n)

