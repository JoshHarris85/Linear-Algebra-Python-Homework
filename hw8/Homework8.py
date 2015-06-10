""" 
by Josh Harris for CS 3513
Homework8.py 
Homework 8
Due 4/16/14
Created on Sat Apr 12 00:21:23 2014
@author: JMH
"""
import numpy as np
import math
from romberg import * 
import random
from trapezoid import *

def f1(x): return (1 / (1 + x**2))
def f2(x): return (x * (math.sqrt(x**2 + 9.0)))
#Had help understanding the wiki formula from my Brother in law to know what I should code.
#Source: https://en.wikipedia.org/wiki/Monte_Carlo_integration#Overview

def MonteCarlo(f, a, b, n):
    if a < b:
        samples = []
        tot = 0
        yMax = 0
        for i in range(n):
            samples.append((random.random() * (b - a) + a))
        for y in samples: 
            yMax = f(y)
            tot += yMax
        return (b - a) * tot / float(n)
    
result11 = MonteCarlo(f1,0,1,1000)
result12 = MonteCarlo(f1,0,1,10000)
result13 = MonteCarlo(f1,0,1,1000000)
result21 = MonteCarlo(f2,0,4,1000)
result22 = MonteCarlo(f2,0,4,10000)
result23 = MonteCarlo(f2,0,4,1000000)
print("The solution for the 1st equation using Monte Carlo at 10^2 is", result11)
print("The solution for the 1st equation using Monte Carlo at 10^4 is", result12)
print("The solution for the 1st equation using Monte Carlo at 10^6 is", result13)
print("The solution for the 2nd equation using Monte Carlo at 10^2 is", result21)
print("The solution for the 2nd equation using Monte Carlo at 10^4 is", result22)
print("The solution for the 2nd equation using Monte Carlo at 10^6 is", result23, "\n")

def check1():
    I1, n1 = romberg(f1, 0, 1)
    I2, n2 = romberg(f2, 0, 4)
    print("Comparing Solutions\n")
    print("The integral for the 1st equation using romberg's is:", I1, "and the number of evaluations are:", n1 )
    print("The integral for the 2nd equation using romberg's is:", I2, "and the number of evaluations are:", n2, "\n" )
    print("The percentage error in the 1st solution for 10^2 is", math.fabs((((I1 - result11) / I1 )*100)))
    print("The percentage error in the 1st solution for 10^4 is", math.fabs((((I1 - result12) / I1 )*100)))
    print("The percentage error in the 1st solution for 10^6 is", math.fabs((((I1 - result13) / I1 )*100)))
    print("The percentage error in the 2nd solution for 10^2 is", math.fabs((((I2 - result21) / I2 )*100)))
    print("The percentage error in the 2nd solution for 10^4 is", math.fabs((((I2 - result22) / I2 )*100)))
    print("The percentage error in the 2nd solution for 10^6 is", math.fabs((((I2 - result23) / I2 )*100)))
    print("\n")
check1()
    

 # For equation 1 
print("Output has been set to a precision to make the table neat. Problem 2 for Equation One:", "\n")   
def problem21(): 
    Iold = 0.0
    a = 0.0 
    b = 1.0
    integral, npan = romberg(f1, 0, 1)
    print("n", "\t", "h", "\t", "T_h(f)", "\t", "e_h", "\t", "e_2h/e_h")
    for k in range(1,9):
        Inew = trapezoid(f1, a, b, Iold, k)
        if (k > 1) and (abs(Inew - Iold)) < 1.0e-10: break
        eh2 = ((integral - Iold) / (integral - Inew))
        Iold = Inew
        n1 = (2**(k-1))
        print(2**(k-1), "\t", '{0:.4f}'.format(((b-a)/n1)), "\t", '{0:.4f}'.format(Iold), "\t",'{0:.4f}'.format(math.fabs(integral - Iold)), "\t", '{0:.4f}'.format((eh2)), "\n") 
problem21()
        
print("Output has been set to a precision to make the table neat. Problem 2 for Equation Two:", "\n")          
# For equation 2
def problem22(): 
    Iold = 0.0
    a = 0.0 
    b = 1.0
    integral, npan = romberg(f2, 0, 4)
    print("n", "\t", "h", "\t", "T_h(f)", "\t", "e_h", "\t", "e_2h/e_h")
    for k in range(1,9):
        Inew = trapezoid(f1, a, b, Iold, k)
        if (k > 1) and (abs(Inew - Iold)) < 1.0e-10: break
        eh2 = ((integral - Iold) / (integral - Inew))
        Iold = Inew
        n1 = (2**(k-1))
        print(2**(k-1), "\t", '{0:.4f}'.format(((b-a)/n1)), "\t", '{0:.4f}'.format(Iold), "\t",'{0:.4f}'.format(math.fabs(integral - Iold)), "\t", '{0:.4f}'.format((eh2)), "\n") 


problem22()        




def rombergFix(f,a,b,tol=1.0e-6):
    
    def richardson(r,k):
        for j in range(k-1,0,-1):
            const = 4.0**(k-j)
            r[j] = (const*r[j+1] - r[j])/(const - 1.0)
        return r
     
    r = np.zeros(21)
    r[1] = trapezoid(f,a,b,0.0,1)
    r_old = r[1]
    for k in range(2,23):
        r[k] = trapezoid(f,a,b,r[k-1],k)
        r = richardson(r,k)
        if abs(r[1]-r_old) < tol*max(abs(r[1]),1.0):
            for i in range(k):
                return r,2**(k-1)
        r_old = r[1]
    print("Romberg quadrature did not converge")
    
xArr, Ite = rombergFix(f1, 0, 1, 10**-12)       
print("The Original Array from Romberg for equation 1\n", xArr, "\n")


err = np.array([[ (xArr[8] - xArr[7]), 0, 0],
                [(xArr[8] - xArr[6]), (xArr[8] - xArr[5]), 0], 
                [(xArr[8] - xArr[4]), (xArr[8] - xArr[3]), (xArr[8] - xArr[2])]])
                
err2h = np.array([[ (err[0,1]/err[1,0]), 0, 0],
                [(err[1,0]/err[2,0]), (err[1,1]/err[2,1]), 0], 
                [(err[2,0]/xArr[8]), (err[2,1]/xArr[8]), (err[2,2]/xArr[8])]])
print("The error ratios of equation 1\n", err2h)       


xArr2, Ite2 = rombergFix(f2, 0, 4, 10**-12)  
print("\nThe Original Array from Romberg for equation 2\n", xArr2, "\n")      
err2 = np.array([[ (xArr2[8] - xArr2[7]), 0, 0],
                [(xArr2[8] - xArr2[6]), (xArr2[8] - xArr2[5]), 0], 
                [(xArr2[8] - xArr2[4]), (xArr2[8] - xArr2[3]), (xArr2[8] - xArr2[2])]])
                
err2h2 = np.array([[ (err2[0,1]/err2[1,0]), 0, 0],
                [(err2[1,0]/err2[2,0]), (err2[1,1]/err2[2,1]), 0], 
                [(err2[2,0]/xArr[8]), (err2[2,1]/xArr2[8]), (err2[2,2]/xArr2[8])]])
print("The error ratios of equation 2\n",err2) 
        
