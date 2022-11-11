#!/usr/bin/env python
# coding: utf-8

# In[31]:


import numpy as np
import matplotlib.pyplot as plt
from sympy import*

def trapezoidal(func,a,b,h):
    a=float(a)
    b=float(b)
    integration=h*func(a)/2
    a+=h
    
    while a<b:
        integration+=h*func(a)
        a+=h
    integration+=h*func(b)/2
    
    return integration

# when a list of integrations is passed to this function (h_o_i),of same order of accuracy, it returns the highest order accuracy
# integration using the Richardson Interpolation.

def h_o_i(l):
#     we got O(h2) integrations
    n=len(l)
    num=n
    
    for i in range(2,2*n,2): # i the order of accuracy at a given time
        new_l=[]
        for j in range(num-1):
            x=((2**i)*l[j+1] - l[j])/(2**i-1)
            new_l.append(x)
        l=new_l
        num-=1
    return l[0]
    

def Romberg(func,a,b,max_err):
    a=float(a)
    b=float(b)
    x=symbols('x')
    func=lambdify(x,func)
    intervals=0
    h=b-a
    
    this_list=[trapezoidal(func,a,b,h)] # first it contains the values with O(h2) error
    prev_integration=h_o_i(this_list)
    err=10.2

    while abs(err)>abs(max_err):
        h=h/2
        this_list.append(trapezoidal(func,a,b,h))
        
        integration=h_o_i(this_list)
        err=integration-prev_integration
        prev_integration=integration
        intervals+=1
        
    return (prev_integration,intervals,err)

def Gauss(func,a,b,max_err):
    a=float(a)
    b=float(b)
    x=symbols('x')
    f=lambdify(x,func)
    intervals=0
    
    o=(b-a)/2
    p=(b+a)/2
    
    # 1-point quadrature
    prev=2*f(0*o+p)
    intervals+=1
    # 2-point quadrature
    integration=f(-0.57735*o+p)+f(0.57735*o+p)
    intervals+=1
    err=integration-prev
    prev=integration
    if abs(err)<abs(max_err):
        return (prev*o,intervals,err)
    # 3-point quadrature
    integration=0.88889*f(0*o+p) + 0.55556*f(-0.77460*o+p) + 0.55556*f(0.77460*o+p)
    intervals+=1
    err=integration-prev
    prev=integration
    if abs(err)<abs(max_err):
        return (prev*o,intervals,err)
    # 4-point quadrature
    integration=0.65215*f(0.33998*o+p) + 0.65215*f(-0.33998*o+p) + 0.34785*f(0.86114*o+p) + 0.34785*f(-0.86114*o+p)
    intervals+=1
    err=integration-prev
    prev=integration
    if abs(err)<abs(max_err):
        return (prev*o,intervals,err)
    # 5-point quadrature
    integration=0.56889*f(0*o+p) + 0.47863*f(0.53847*o+p) + 0.47863*f(-0.53847*o+p) + 0.23693*f(0.90618*o+p) + 0.23693*f(-0.90618*o+p)
    intervals+=1
    err=integration-prev
    prev=integration
    if abs(err)<abs(max_err):
        return (prev*o,intervals,err)
    return (prev*o,intervals,err)

with open("q1.txt", "r") as f:
        func = f.readline()
        a, b = map(float, f.readline().split(' '))
        max_err = float(f.readline())

method=input("Which method you want to use?\n1. Romberg Integration\n2. Gauss-Legendre Quadrature\n")
I=0 # Integration value
i=0 # Number of Intervals
apre=0 # Approximate relative error (%)

if method=="1":
    (I,i,apre)=Romberg(func,a,b,max_err/100)
elif method=="2":
    (I,i,apre)=Gauss(func,a,b,max_err/100)
else:
    print("Incorrect input")

x=symbols('x')
func=lambdify(x,func)
x_vals=np.linspace(a,b+0.00001,5)
y_vals=[func(i) for i in x_vals]
plt.scatter(x_vals,y_vals,color='black')
plt.grid()
plt.xlabel("x")
plt.ylabel("y")

print(f"\nI = {I}\nNumber of Intervals = {i}\nApproximate relative error (%) = {round(apre*100,6)}")

