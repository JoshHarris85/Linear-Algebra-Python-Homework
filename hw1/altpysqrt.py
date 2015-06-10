"""
altpysqrt.py
by Josh Harris for CS 3513
Homework 1
Due 1/29/14

Created on Tue Jan 28 21:38:58 2014
Title: pysqrt.py
Author:Doug Heisterkamp
Date: 1/29/14
Availability: https://webmail.cs.okstate.edu/svn/cs3513/2014Spring/Examples/background/pysqrt.py
All code used from Doug Heisterkamp except line 31
@author: JMH
"""

def sqrtIter(a,xinit,eps=1e-12,maxIter=15,iterSeq=None):
   """Example function that calculates the square root of a value
      a --- function calculates square root of a
      xinit --initial guess
      eps --- error tolerance using absolute error
      maxIter -- maximum number of iterations
   """
   count = 0
   xnew = xinit
   if iterSeq != None:
      iterSeq.append(xnew)
   successFlag = False
   while count < maxIter :
      count += 1
      xold = xnew
      xnew = (((xold * xold * xold) + (3 * xold * a)) / (3 * (xold * xold) + a))
      if iterSeq != None:
         iterSeq.append(xnew)
      if abs(xnew-xold) <  eps:
         successFlag = True
         break
   return xnew,successFlag,count

if __name__ == "__main__" :
   print("Test output")
   s,flag,c = sqrtIter(9,5,1e-10,100)
   print("sqrtIter(9,5,1e-10,100) = {} with flag {} and iteration count {}".format(s,flag,c))
   s,flag,c = sqrtIter(2,2,1e-10,100)
   print("sqrtIter(2.0,2.0) = {} with flag {} and iteration count {}".format(s,flag,c))
   slist = []
   s,flag,c = sqrtIter(15,5,iterSeq=slist)
   print("sqrtIter(15,5,iterSeq=slist) = {} with flag {} and iteration count {}".format(s,flag,c))
   print("    slist = {}".format(slist))

   # a snippet of code to wait for a key press before exiting in MS Windows
   import platform
   if "Win" in platform.uname()[0] :
      print("Press a key to exit")
      input()
