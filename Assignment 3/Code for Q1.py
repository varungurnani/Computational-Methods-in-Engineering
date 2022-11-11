#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import numpy.polynomial.polynomial as P
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def lagrange(n,known_pairs):
    coefficients=np.array([0])
        
    for i in range(n):
        coef=np.array([1])
        for j in range(n):
            if j==i:
                continue
            
            coef=P.polymul(coef,(1/(known_pairs[i][0]-known_pairs[j][0]),-known_pairs[j][0]/(known_pairs[i][0]-known_pairs[j][0])))
        coefficients=P.polyadd(coefficients,P.polymul(coef,(known_pairs[i][1])))
        
    return coefficients


# Two ways to take input, 1st-> from .txt file, 2nd-> from user in runtime    
with open("<write the whole address of input file here>\q1.txt", "r") as f:
    n = int(f.readline())
    known_pairs = np.array([list(map(float, f.readline().split())) for i in range(n)])
    m = int(f.readline())
    to_find = np.array([float(f.readline()) for i in range(m)])
    
# n=int(input("Enter the number of values:"))
# known_pairs=[]
# print("Enter them one pair at a time separated by space:")
# for i in range(n):
#     known_pairs.append(list(map(float,input().split())))
# m=int(input("Enter the number of values to be predicted:"))
# to_find=[]
# print("Enter them one pair at a time separated by space:")
# for i in range(m):
#     to_find.append(float(input()))


print("Interpolated values y at given x*\n")

# Lagrange polynomials
p=lagrange(n,known_pairs)
x=np.arange(min(to_find)-0.2,max(to_find)+0.2,0.001)
f=np.poly1d(p)
ans=[f(i) for i in to_find]
plt.plot(x,np.polyval(p,x),label="Lagrange",c='blue')
t=[known_pairs[i][0] for i in range(n)]
u=[known_pairs[i][1] for i in range(n)]
plt.scatter(t,u,label="Data",c='black')
plt.grid()
print("Lagrange Polynomials")
for i in range(m):
    print(f"{to_find[i]}     {round(np.poly1d(p)(to_find[i]),4)}")

# Normal Cubic Spline
ns = CubicSpline(t,u,bc_type='natural', extrapolate=True)
plt.plot(x,ns(x),label="Spline",c='green')
plt.legend()
print("\nCubic Spline")
for i in range(m):
    print(f"{to_find[i]}     {round(float(ns(to_find[i])),4)}")

