# Programming Assignment 1: Non-Linear Equations

from sympy import *
import numpy as np
import matplotlib.pyplot as plt
def derivative(f,x):
    dx = 1e-6
    df = f(x+dx)-f(x-dx)
    return df/(2*dx)
def bisection(f, xl, xu, max_iter, max_err):
    iter = 0
    err = max_err*2
    xr = xl
    x=[]
    y=[]
    while(iter<max_iter and err>max_err):
        xr_prev = xr
        xr = (xl+xu)/2
        iter = iter+1
        if(xr!=0):
            err=abs(xr-xr_prev)/xr * 100
        x.append(iter)
        y.append(err)
        test = f(xl) * f(xr)
        if(test<0):
            xu = xr
        elif (test>0):
            xl = xr
        elif (test==0):
            err = 0
        if(iter>max_iter or err<max_err):
            break
    x1 = np.linspace(-10, 10, 100)
    y1 = f(x1)
    plt.plot(x1, y1,'o-m')
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("Plot of f(x) vs x")
    plt.show()
    plt.plot(x,y,'o-m')
    plt.xlabel("Iteration Number")
    plt.ylabel("Relative approximate error")
    plt.title("Relative approximate error vs Iteration Number")
    plt.show()
    return xr  
def false_position(f, xl, xu, max_iter, max_err):
    iter = 0
    err = max_err*2
    xr = xl
    x=[]
    y=[]
    while(iter<max_iter and err>max_err):
        xr_prev = xr
        xr = xu - f(xu)*(xl-xu)/(f(xl)-f(xu))
        iter = iter+1
        if(xr!=0):
            err=abs(xr-xr_prev)/xr * 100
        x.append(iter)
        y.append(err)
        test = f(xl) * f(xr)
        if(test<0):
            xu = xr
        elif (test>0):
            xl = xr
        elif (test==0):
            err = 0
        if(iter>max_iter or err<max_err):
            break
    x1 = np.linspace(-5, 5, 100)
    y1 = f(x1)
    plt.plot(x1, y1,'o-m')
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("Plot of f(x) vs x")
    plt.show()
    plt.plot(x,y,'o-m')
    plt.xlabel("Iteration Number")
    plt.ylabel("Relative approximate error")
    plt.title("Relative approximate error vs Iteration Number")
    plt.show()
    return xr 
def fixed_point(f, phi, xl, max_iter, max_err):
    a=xl
    rae=[] #relative approximate error values
    for i in range(max_iter):
        error=phi(a)-a
        rae.append(error)
        if(abs(error)<=max_err):
            a=phi(a)
            break
        a=phi(a)

    x1 = np.linspace(-5, 5, 100)
    y1 = f(x1)
    plt.plot(x1, y1,'o-m')
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("Plot of f(x) vs x")
    plt.show()
    plt.plot(np.arange(1,len(rae)+1,1),rae,'o-m')
    plt.xlabel("Iteration Number")
    plt.ylabel("Relative approximate error")
    plt.title("Relative approximate error vs Iteration Number")
    plt.show()
    return a
def newton_raphson(f, x0, max_iter, max_err, derivative_of_f):
    iter = 0
    err = max_err*2
    x2=[]
    y2=[]
    while(iter<max_iter and err>max_err):
        x0_prev = x0
        x0 = x0 - f(x0)/derivative(f, x0)
        iter = iter+1
        if(x0!=0):
            err=abs(x0-x0_prev)/x0 * 100
        x2.append(iter)
        y2.append(err)
        if(iter>max_iter or err<max_err):
            break
    x1 = np.linspace(-5, 5, 100)
    y1 = f(x1)
    plt.plot(x1, y1,'o-m')
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("Plot of f(x) vs x")
    plt.show()
    plt.plot(x2,y2,'o-m')
    plt.xlabel("Iteration Number")
    plt.ylabel("Relative approximate error")
    plt.title("Relative approximate error vs Iteration Number")
    plt.show()
    return x0
def secant(f, x_minus1, x0, max_iter, max_err):
    iter = 0
    err = max_err*2
    x2=[]
    y2=[]
    while(iter<max_iter and err>max_err):
        x_prev = x0
        x0 = x0 - f(x0)*(x_minus1-x0)/(f(x_minus1)-f(x0))
        x_minus1 = x_prev
        iter = iter+1
        if(x0!=0):
            err=abs(x0-x_minus1)/x0 * 100
        x2.append(iter)
        y2.append(err)
        if(iter>max_iter or err<max_err):
            break
    x1 = np.linspace(-5, 5, 100)
    y1 = f(x1)
    plt.plot(x1, y1,'o-m')
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title("Plot of f(x) vs x")
    plt.show()
    plt.plot(x2,y2,'o-m')
    plt.xlabel("Iteration Number")
    plt.ylabel("Relative approximate error")
    plt.title("Relative approximate error vs Iteration Number")
    plt.show()
    return x0 
func = input("Enter the function f(x) = ")
x = Symbol('x')
f = lambdify(x, func)
print("1. Bisection Method")
print("2. False-position Method")
print("3. Fixed-Point Method")
print("4. Newton-Raphson Method")
print("5. Secant Method")
serial = int(input("Enter the serial number of any of the above methods to use them: "))
if(serial==1):
    xl = float(input("Enter lower limit: "))
    xu = float(input("Enter upper limit: "))
    max_iter = int(input("Enter maximum number of iterations: "))
    max_err = float(input("Enter maximum permissible error: "))
    bisect = bisection(f, xl, xu, max_iter, max_err)
    print("Approximate solution of the equation is", bisect)
if(serial==2):
    xl = float(input("Enter lower limit: "))
    xu = float(input("Enter upper limit: "))
    max_iter = int(input("Enter maximum number of iterations: "))
    max_err = float(input("Enter maximum permissible error: "))
    false = false_position(f, xl, xu, max_iter, max_err)
    print("Approximate solution of the equation is", false)
if(serial==3):
    func = input("Enter the function phi(x) = ")
    phi = lambdify(x, func)
    xl = float(input("Enter initial estimate: "))
    max_iter = int(input("Enter maximum number of iterations: "))
    max_err = float(input("Enter maximum permissible error: "))
    fixed = fixed_point(f, phi, xl, max_iter, max_err)
    print("Approximate solution of the equation is", fixed)
if(serial==4):
    func = input("Enter the function f'(x) = ")
    derivative_of_f = lambdify(x, func)
    x0 = float(input("Enter initial guess value: "))
    max_iter = int(input("Enter maximum number of iterations: "))
    max_err = float(input("Enter maximum permissible error: "))
    nr = newton_raphson(f, x0, max_iter, max_err, derivative_of_f)
    print("Approximate solution of the equation is", nr)
if(serial==5):
    x_minus1 = float(input("Enter first value of initial estimate: "))
    x0 = float(input("Enter second value of initial estimate: "))
    max_iter = int(input("Enter maximum number of iterations: "))
    max_err = float(input("Enter maximum permissible error: "))
    sec = secant(f, x_minus1, x0, max_iter, max_err)
    print("Approx. solution of the equation f(x) is", sec)

