#!/usr/bin/env python
# coding: utf-8

# In[46]:


import numpy
import cmath
from sympy import *
from numpy import *
import math
import matplotlib.pyplot as plt
numpy.set_printoptions(precision=4,threshold=10000,linewidth=150,suppress=True)

def plotof(eq):
    x=Symbol('x')
    f=lambdify(x,eq,"numpy")
    x_vals = linspace(float(input('Enter the lower limit of graph of x axis: ')),float(input('Enter the upper limit of graph of y axis: ')), 100)
    y_vals = f(x_vals)
    plt.plot(x_vals, y_vals,'o-m')
    plt.title('Plot of F(x) vs x')
    plt.legend(['y='+ eq])
    plt.grid()
    plt.show()
    return 0

def muller(eq):
    x=Symbol('x')
    f=lambdify(x,eq,"numpy")
    x=numpy.array([float(input('Please enter value of x0: \n')),float(input('Please enter value of x1: \n')),float(input('Please enter the value of x2: \n'))])
    fx=f(x)
    max_iter = int(input('Enter the Maximum number of Iterations:'))
    max_error = float(input('Enter the Maximum Error:'))
    errors=[]
    for iteration in range(0,max_iter):
        h0=x[1]-x[0]
        h1=x[2]-x[1]
        delta0=(fx[1]-fx[0])/(h0)
        delta1=(fx[2]-fx[1])/(h1)
        a=(delta1-delta0)/(h1+h0)
        b=a*h1+delta1
        x[0]=x[1]
        fx[0]=fx[1]
        x[1]=x[2]
        fx[1]=fx[2]
        x[2]=x[2]+(-2*fx[2])/(b+sign(b)*math.sqrt(b*b-4*a*fx[2]))
        fx[2]=f(x[2])
        err=abs(((x[2]-x[1])/x[2])*100) #Relative percentage error
        errors.append(err)
        if fx[2]==0:
            break
        elif (err<max_error ):
            break
    print(f"The root of provided equation is {x[2]} ")
    plotof(eq)
    return 0

rootsbairstow=[] #List to store the roots of the given polynomial.
def bairstow(equation,r,s,max_iter,max_error):
    x=Symbol('x')
    f=lambdify(x,equation,"numpy")
    errorsR=[]
    errorsS=[]
    degreefx=0
    for i in range(0,max_iter):
        b=[]
        c=[]
        fx=Poly(equation,x)
        degreefx=degree(fx)
        alpha,beta=fx.div(Poly(x**2-r*x-s,x))# Function to calucate the values of B0 to Bn through division.
        b=alpha.all_coeffs()
        b1b0=beta.all_coeffs()
        b1b0[1]=b1b0[1]+b1b0[0]*r
        b=b+b1b0
        for j in range(0,degree(fx)):
            if len(c)==0:
                c.append(b[j]) #calcualting value of Cn.
            elif len(c)==1:
                c.append(b[j]+r*c[0]) #calcualting value of Cn-1.
            else:
                c.append(b[j]+r*c[j-1]+s*c[j-2]) #calcualting values of Cn-2 to C1.
        A=numpy.array([[c[len(c)-2],c[len(c)-3]],[c[len(c)-1],c[len(c)-2]]], dtype='float')
        B=numpy.array([-b[len(b)-2],-b[len(b)-1]], dtype='float')
        deltaR,deltaS=numpy.linalg.solve(A,B)
        r=r+deltaR
        s=s+deltaS
        EaR=abs(deltaR/r)*100
        EaS=abs(deltaS/s)*100
        errorsR.append(EaR)
        errorsS.append(EaS)
        if EaS<max_error and EaR<max_error:
            break
    rootsbairstow.append((r+cmath.sqrt(r*r+4*s))/2)
    rootsbairstow.append((r-cmath.sqrt(r*r+4*s))/2)
    if degreefx>4:
        bairstow(alpha,r,s,max_iter,max_error) # A recursive funtion.
    elif degreefx==4:
        rootsbairstow.append(complex((-b[1]+cmath.sqrt(b[1]*b[1]-4*b[2]*b[0]))/2*b[0]))
        rootsbairstow.append(complex((-b[1]-cmath.sqrt(b[1]*b[1]-4*b[0]*b[2]))/2*b[0]))
    elif degreefx==3:
        al,be=fx.div(Poly(x**2-r*x-s,x))
        bl=al.all_coeffs()
        rootsbairstow.append(complex(-1*bl[1]/bl[0],0))
    return 0

function = input('Please enter your  function in following format (eg. 7*x**4 - 75*x**3 + 620*x**2 - 29*x - 10): \n')
method=input('Please select the method (by typing the number from below) you want to use:\nMuller : 1\nBairstow : 2\n')

if method.upper()=='1':
    muller(function)
elif method.upper()=='2':
    first_estimate=float(input('Please enter value of first_estimate: '))
    second_estimate=float(input('Please enter value of second_estimate: '))
    maxi_iter = int(input('Enter the Maximum number of Iterations: '))
    maxi_error= float(input('Enter the Maximum Error: '))
    bairstow(function,first_estimate,second_estimate,maxi_iter,maxi_error)
    print('Roots of given',function,'are: ')
    for item in rootsbairstow:
        print(item)
    plotof(function)
else:
    print('Error: The option is not correct')

