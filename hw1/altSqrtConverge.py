"""
altSqrtConverge.py
by Josh Harris for CS 3513
Homework 1
Due 1/29/14

Created on Tue Jan 28 21:38:58 2014
Title: sqrtConverge.py
Author:Doug Heisterkamp
Date: 1/29/14
Availability: https://webmail.cs.okstate.edu/svn/cs3513/2014Spring/Examples/background/sqrtConverge.py
All code used from Doug Heisterkamp
"""

import altpysqrt      # our module that implements the iterative square root algorithm
import decimal     # for higher precision floating point
import numpy as np # for arrays

# set precision of decimals to a high number 
decimal.getcontext().prec = 300
a = decimal.getcontext().exp(2) # e^2
b = decimal.getcontext().exp(1) # e^1 = sqrt(e^2)
tol = decimal.Decimal(1e-299)

print("Calculating sqrt(e^2).")
slist = []
s,flag,c = altpysqrt.sqrtIter(a,a,eps=tol,maxIter=100, iterSeq=slist)
if not flag :
   print("Failed to converge.")
else:
   print("Error in final value is {}".format(b-s))
   # convert list to array for manipulations.  
   # Note: sl is an array of Decimal objects
   sl = np.array(slist)
   # use the true value to calc convergence
   e_n = abs(sl - b)    
   zflag = False
   print("Using true value to calc error.")
   print("{:4} {:>15}  {:>15}  {:>15}  {:>15}".format("iter","|e_n|", 
         "|e_n+1|/|e_n|^1","|e_n+1|/|e_n|^2","|e_n+1|/|e_n|^3"))
   for i in range(0,len(e_n)-1):
      if e_n[i+1] == 0:
         print("{:4d} {:>15.2e}".format(i,e_n[i]))
         print("e_n[{}] is zero, skipping rest of list.".format(i+1))
         zflag = True
         break
      e2 = e_n[i] * e_n[i]
      e3 = e2 * e_n[i]
      print("{:4d} {:>15.3e}  {:>15.3e}  {:>15.3e}  {:>15.3e}".format(i,e_n[i], 
             e_n[i+1]/e_n[i],e_n[i+1]/e2,e_n[i+1]/e3))
   if not zflag:
      print("{:4d} {:>15.3e}".format(len(e_n)-1,e_n[-1]))


# a snippet of code to wait for a key press before exiting in MS Windows
import platform
if "Win" in platform.uname()[0] :
   print("Press a key to exit")
   input()
