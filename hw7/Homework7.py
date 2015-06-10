""" 
by Josh Harris for CS 3513
Homework7.py 
Homework 7
Due 3/12/14
Created on Sun Mar 9 14:12:22 2014
@author: JMH
"""


import numpy as np
from rootsearch import *
from ridder import * 
from newtonRaphson import *
from newtonRaphson2 import *

#F(x) = 0 in a given interval with ridder. Test by finding the roots of x sin x + 3 cos x - x = 0 (-6,6)
#Problem 10 Page 166
def f(x): return x*math.sin(x) + 3*math.cos(x) - x
def problem10():
    a,b,dx = (-6, 6, 0.01)
    print("The roots for problem 10 are:")
    while True:
        x1,x2 = rootsearch(f,a,b,dx)
        if x1 != None:
            a = x2
            root = ridder(f,x1,x2)
            if root != None: print(root)
        else:
            print("Done\n")
            break
    #input("Press return to exit")
problem10()

#Problem 11 on page 166
def df11(x): return -2*math.sin(x) + x*math.cos(x) - 1
def problem11():
    a,b,dx = (-6, 6, 0.01)
    print("The roots for problem 11 are:")
    while True:
        x1,x2 = rootsearch(f,a,b,dx)
        if x1 != None:
            a = x2
            root = newtonRaphson(f,df11,x1,x2)
            if root != None: print(root)
        else:
            print("Done\n")
            break
    #input("Press return to exit")
problem11()

#Extra Credit
#Solving Taylor's inside the extra credit update formula equation.txt included with this code
def f1(x): return x*math.sin(x) + 3*math.cos(x) - x
def df(x): return -2*math.sin(x) + x*math.cos(x) - 1
def df2(x): return -x*math.sin(x) - math.cos(x)
def extraCredit(x,tol=1.0e-9):
    for i in range(30):
        dx = -((f1(x) * 2*df(x)) / (2*df(x)**2-f1(x)*df2(x)))
        x = x + dx
        if abs(dx) < tol: return x, i
    print("Too many iterations\n")
root, numIter = extraCredit(-5)
print("Extra Credit: Using the data from problem 11:\nThe root is", root, "\nNumber of iterations are", numIter)
root2, numIter2 = extraCredit(-3)
print("The root is", root2, "\nNumber of iterations are", numIter2)
root3, numIter3 = extraCredit(2)
print("The root is", root3, "\nNumber of iterations are", numIter3, "\n")

#Through algebra set y = equation then substitute opposite equation into the first
def f23(x): return 2704*x**4 -13520*x**3 + 22984*x**2 -14820*x +2625  
def df23(x): return 52*(208*x**3 - 780*x**2 + 884*x - 285)
def problem23():
    a,b,dx = (-2, 2, 0.001)
    print("Problem 23: The coordinates where the two points of the circle intersect are:")
    while True:
        x1,x2 = rootsearch(f23,a,b,dx)
        if x1 != None:
            a = x2
            root = newtonRaphson(f23,df23,x1,x2)
            if root != None: print(root)
        else:
            print("Done\n")
            break
problem23()


# I did not have enough time to complete problem 26 

