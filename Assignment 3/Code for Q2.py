#!/usr/bin/env python
# coding: utf-8

# In[64]:


import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import polynomial as P

def Fit(Degree,pairs):
    A_matrix=np.zeros((Degree+1,Degree+1))
    for i in range(Degree+1):
        for j in range(Degree+1):
            power=i+j
            sum=0
            for k in range(n):
                sum+=pow(pairs[k][0],power)
            A_matrix[i][j]=sum
    B_matrix=np.zeros((Degree+1,1))
    for i in range(Degree+1):
        power=i
        sum=0
        for k in range(n):
            sum+=(pow(pairs[k][0],power))*pairs[k][1]
        B_matrix[i][0]=sum

    x=np.linalg.inv(A_matrix).dot(B_matrix)

    return x
    

# Two ways to take input, 1st-> from .txt file, 2nd-> from user in runtime
with open("<write the whole address of input file here>\q2.txt", "r") as f:
    n = int(f.readline())
    pairs = np.array([list(map(float, f.readline().split())) for i in range(n)])
    
# n=int(input("Enter the number of values:"))
# pairs=[]
# print("Enter them one pair at a time separated by space:")
# for i in range(n):
#     pairs.append(list(map(float,input().split())))

# plotting the curves
t=[pairs[i][0] for i in range(n)]
u=[pairs[i][1] for i in range(n)]
x=np.arange(0,1,0.001)
plt.scatter(t,u,label="Data",c='black')
plt.grid()


l=Fit(1,pairs)

# finding R2-Score
sumx=0
for i in u:
    sumx+=i
mean=sumx/len(u)
Var_mean=0
for i in u:
    Var_mean+=((i-mean)**2)/len(t)
Var_fit=0
for i in range(len(t)):
    Var_fit+=((u[i]-np.poly1d([l[1][0],l[0][0]])(t[i]))**2)/len(t)
Rsq=1-Var_fit/Var_mean
    
print("Linear : coefficients   ",round(l[0][0],3)," ",round(l[1][0],3),"\n      : R-sq = ",round(Rsq,3),"\n")
plt.plot(x,np.poly1d([l[1][0],l[0][0]])(x),c="yellow",label="Linear")

l=Fit(2,pairs)

Var_fit=0
for i in range(len(t)):
    Var_fit+=((u[i]-np.poly1d([l[2][0],l[1][0],l[0][0]])(t[i]))**2)/len(t)
Rsq=1-Var_fit/Var_mean

print("Quadratic : coefficients   ",round(l[0][0],3)," ",round(l[1][0],3)," ",round(l[2][0],3),"\n      : R-sq = ",round(Rsq,3),"\n")
plt.plot(x,np.poly1d([l[2][0],l[1][0],l[0][0]])(x),c="green",label="Quadratic")

l=Fit(3,pairs)

Var_fit=0
for i in range(len(t)):
    Var_fit+=((u[i]-np.poly1d([l[3][0],l[2][0],l[1][0],l[0][0]])(t[i]))**2)/len(t)
Rsq=1-Var_fit/Var_mean

print("Cubic : coefficients   ",round(l[0][0],3)," ",round(l[1][0],3)," ",round(l[2][0],3)," ",round(l[3][0],3),"\n      : R-sq = ",round(Rsq,3),"\n")
plt.plot(x,np.poly1d([l[3][0],l[2][0],l[1][0],l[0][0]])(x),c="purple",label="Cubic")

l=Fit(4,pairs)

Var_fit=0
for i in range(len(t)):
    Var_fit+=((u[i]-np.poly1d([l[4][0],l[3][0],l[2][0],l[1][0],l[0][0]])(t[i]))**2)/len(t)
Rsq=1-Var_fit/Var_mean

print("Quartic : coefficients   ",round(l[0][0],3)," ",round(l[1][0],3)," ",round(l[2][0],3)," ",round(l[3][0],3)," ",round(l[4][0],3),"\n      : R-sq = ",round(Rsq,3),"\n")
plt.plot(x,np.poly1d([l[4][0],l[3][0],l[2][0],l[1][0],l[0][0]])(x),c="blue",label="Quartic")

plt.legend();

