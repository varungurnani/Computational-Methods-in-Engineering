#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np
import matplotlib.pyplot as plt
from sympy import*

def forward_Euler(func,t_0,y_0,t_f,h):
    t=symbols('t')
    y=symbols('y')
    func=lambdify([t,y],func)
    
    t_values=[]
    y_values=[]
    t_values.append(t_0)
    y_values.append(y_0)
    
    while t_0<t_f:
        print(f"{round(t_0,1)} {round(y_0,8)}\n")
        y_0 = y_0 + h*func(t_0,y_0)
        t_0 = t_0+h
        t_values.append(t_0)
        y_values.append(y_0)
    plt.scatter(t_values,y_values,color='black')
        
#     return y_0

def RK_2nd(func,t_0,y_0,t_f,h):
    t=symbols('t')
    y=symbols('y')
    func=lambdify([t,y],func)
    
    t_values=[]
    y_values=[]
    t_values.append(t_0)
    y_values.append(y_0)
    
    while t_0<t_f:
        print(f"{round(t_0,1)} {round(y_0,8)}\n")
        y_0 = y_0 + h*func(t_0+h/2,y_0+(h/2)*func(t_0,y_0))
        t_0 = t_0+h
        t_values.append(t_0)
        y_values.append(y_0)
    plt.scatter(t_values,y_values,color='black')
#     return y_0
    
def RK_4th(func,t_0,y_0,t_f,h):
    t=symbols('t')
    y=symbols('y')
    func=lambdify([t,y],func)
    
    t_values=[]
    y_values=[]
    t_values.append(t_0)
    y_values.append(y_0)
    
    while t_0<t_f-h/2:
        k1=func(t_0,y_0)
        k2=func(t_0+h/2,y_0+(h/2)*k1)
        k3=func(t_0+h/2,y_0+(h/2)*k2)
        k4=func(t_0+h,y_0+h*k3)
    
        y_0=y_0+(h/6)*(k1 + 2*k2 + 2*k3 + k4)
        t_0=t_0+h
        t_values.append(t_0)
        y_values.append(y_0)
        print(f"{round(t_0,1)} {round(y_0,8)}\n")
    plt.scatter(t_values,y_values,color='black')
    
#     return y_0


with open("q2.txt", "r") as f:
        func = f.readline()
        t_0, y_0 = map(float, f.readline().split(' '))
        t_f = float(f.readline())
        h = float(f.readline())

method=int(input("Which method you want to use?\n1. Forward Euler method\n2. 2nd order RK method (Midpoint method)\n3. 4th order RK method\n"))
print("\nt, y\n")
print(f"{round(t_0,1)} {round(y_0,8)}\n")

if method==1:
    forward_Euler(func,t_0,y_0,t_f,h)
elif method==2:
    RK_2nd(func,t_0,y_0,t_f,h)
elif method==3:
    RK_4th(func,t_0,y_0,t_f,h)
else:
    print("Incorrect input")

plt.grid()
plt.xlabel('t')
plt.ylabel('y')

