"""Testing LU implementations on problem 15 from HW1"""

# Doug Heisterkamp
# last modified: 2/1/2014

import numpy as np
import numpy.linalg as la

def hilbert(n, dtype=float) :
   return np.fromfunction( lambda i,j :1/( i+j+1), (n,n),dtype=dtype)

# repeat problem 15 using LUdecomp, LUpivot, and LUpivotOP

def problem15():
   "problem set 2.1, problem 15, page 57"

   import LUdecomp
   def specializedLU(n):
      A = hilbert(n)
      b = np.sum(A,axis=1) 
      LU = LUdecomp.LUdecomp(A)
      return LUdecomp.LUsolve(LU,b)
 
   x_err = 0
   n=1
   print("Problem 15 using LUdecomp")
   print("{:>5} : {:<12s}".format("n","error, inf norm"))
   while x_err < 1e-6:
      n += 1
      x_true = np.ones(n,dtype=float)
      x_approx = specializedLU(n)
      x_err = la.norm(x_true - x_approx, np.inf)
      print("{:5d} : {:<6.4e}".format(n,x_err))
   print("The largest n with solution within six significant digits is {}".format(n-1))

def problem15a():
   "problem set 2.1, problem 15, page 57"
   import LUpivot

   def specializedLU(n):
      A = hilbert(n)
      b = np.sum(A,axis=1) 
      LU,seq = LUpivot.LUdecomp(A,1e-12)
      return LUpivot.LUsolve(LU,b,seq)
 
   x_err = 0
   n=1
   print("Problem 15 using LUpivot")
   print("{:>5} : {:<12s}".format("n","error, inf norm"))
   while x_err < 1e-6:
      n += 1
      x_true = np.ones(n,dtype=float)
      x_approx = specializedLU(n)
      x_err = la.norm(x_true - x_approx, np.inf)
      print("{:5d} : {:<6.4e}".format(n,x_err))
   print("The largest n with solution within six significant digits is {}".format(n-1))


def problem15b():
   "problem set 2.1, problem 15, page 57"
   import LUpivotOP

   def specializedLU(n):
      A = hilbert(n)
      b = np.sum(A,axis=1) 
      LU,seq = LUpivotOP.LUdecomp(A,1e-12)
      return LUpivotOP.LUsolve(LU,b,seq)
 
   x_err = 0
   n=1
   print("Problem 15 using LUpivotOP")
   print("{:>5} : {:<12s}".format("n","error, inf norm"))
   while x_err < 1e-6:
      n += 1
      x_true = np.ones(n,dtype=float)
      x_approx = specializedLU(n)
      x_err = la.norm(x_true - x_approx, np.inf)
      print("{:5d} : {:<6.4e}".format(n,x_err))
   print("The largest n with solution within six significant digits is {}".format(n-1))

def problem15c():
   "problem set 2.1, problem 15, page 57"
   import LUpivotComplete

   def specializedLU(n):
      A = hilbert(n)
      b = np.sum(A,axis=1) 
      LU,seqR,seqC = LUpivotComplete.LUdecomp(A,1e-12)
      return LUpivotComplete.LUsolve(LU,b,seqR,seqC)
 
   x_err = 0
   n=1
   print("Problem 15 using LUpivotComplete")
   print("{:>5} : {:<12s}".format("n","error, inf norm"))
   while x_err < 1e-6:
      n += 1
      x_true = np.ones(n,dtype=float)
      x_approx = specializedLU(n)
      x_err = la.norm(x_true - x_approx, np.inf)
      print("{:5d} : {:<6.4e}".format(n,x_err))
   print("The largest n with solution within six significant digits is {}".format(n-1))


if __name__ == "__main__" :
   problem15()
   problem15a()
   problem15b()
   problem15c()


