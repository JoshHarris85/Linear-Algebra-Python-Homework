""" 
Hnum.py 
by Josh Harris for CS 3513
Homework 1
Due 1/29/14
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 22:57:19 2014

@author: JMH
"""
from decimal import *


# Part A 
def standardForwardSummation():
    i = 0 
    for j in range(1, 10001):   
        i = (1.0 / j) + i
    return i
    
#Print for function above
print (standardForwardSummation())

# Part B 
getcontext().prec = 5
def fivePrecisionSummation():
    i = 0
    for j in range(1, 10001):   
        i = ( 1 / Decimal(j) + i)
    return (i)
    
#Print for function above
print(fivePrecisionSummation())


# Part C 

def fivePrecisionBackwards():
    i = 0
    for j in range(10000, 0, -1):
       i = (1 / Decimal(j) + i)
    return (i)

#Print for function above
print(fivePrecisionBackwards())


# Part D 
v = []
for i in range(1, 10001):
    v.append(1 / Decimal(i))


def compensatedSum(x):
    n = len(x)
    c = 0
    s = x[0]
    for i in range(1, n-1):
        y = c + x[i]
        t = s + y 
        c = (s - t) + y
        s = t
    return s 


print (compensatedSum(v))
#Print for function above
