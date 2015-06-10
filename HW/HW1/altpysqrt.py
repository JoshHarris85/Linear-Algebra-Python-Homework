"""
Homework 1, CS 3513, Spring 2014

altpysrt uses x_{n+1} = (x_n^3 + 3 x_n a)/( 3 x_n^2 + 1)

"""

# CS 3513
# Spring 2014

# Written by : Doug Heisterkamp
# last modified: 01/26/2014

# modified version of pysqrt

from compat27 import  * # imports to make python version 2.7 behave like 3.2

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
      #xnew = (xold + a/xold)/2
      xnew = (xold**3 + 3*xold *a)/ (3*xold**2 + a)
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
