"""
CS 3513, Homework 2
Problem Set 2.1
Example Solution
"""

# Doug Heisterkamp
# Last modified: 2/5/2014

import numpy as np
import numpy.linalg as la

def problem8():
   "problem set 2.1, problem 8, page 56"
   import LUdecomp
   A = np.array([[-3,6,-4],[9,-8,24],[-12,24,-26]],dtype=float)
   A_orig = A.copy()
   LU = LUdecomp.LUdecomp(A)
   b = np.array([-3,65,-42],dtype=float)
   b_orig = b.copy()
   x = LUdecomp.LUsolve(LU,b)
   # extract L and U for verification
   U = np.triu(LU)  # 
   L = np.tril(LU)
   L[ np.diag_indices_from(L) ] = 1.0 
   print("""
Problem 8:
A = 
{}
LU decomposition A = LU, LU (in one matrix) = 
{}
Solving Ax=b, with b = {}
Solution x = {}
Verifying solution: 
     residual ||Ax-b||_2 = {}
     ||A - dot(L,U)||_inf = {}
""".format(A_orig,LU,b_orig,x, 
   la.norm(np.dot(A_orig,x)-b_orig,2), 
   la.norm(A_orig - np.dot(L,U),np.inf))
   )

def problem11():
   "problem set 2.1, problem 11, page 56"
   import choleski
   A = np.array([[1,1,1],[1,2,2],[1,2,3]],dtype=float)
   A_orig = A.copy()
   L = choleski.choleski(A)
   b = np.array([1,3/2,3],dtype=float)
   b_orig = b.copy()
   x = choleski.choleskiSol(L,b)
   print("""
Problem 11:
A = 
{}
Choleski decomposition A = LL^T, L = 
{}
Solving Ax=b, with b = {}
Solution x = {}
Verifying solution: 
     residual ||Ax-b||_2 = {}
     ||A - LL^T||_inf = {}
""".format(A_orig,L,b_orig,x, 
   la.norm(np.dot(A_orig,x)-b_orig,2), 
   la.norm(A_orig - np.dot(L,L.transpose()),np.inf))
   )


def problem15():
   "problem set 2.1, problem 15, page 57"

   import LUdecomp
   def hilbert(n, dtype=float) :
      return np.fromfunction( lambda i,j :1/( i+j+1), (n,n),dtype=dtype)

   def specializedLU(n):
      A = hilbert(n)
      b = np.sum(A,axis=1) 
      LU = LUdecomp.LUdecomp(A)
      return LUdecomp.LUsolve(LU,b)
 
   x_err = 0
   n=1
   print("Problem 15")
   print("{:>5} : {:<12s}".format("n","error, inf norm"))
   while x_err < 1e-6:
      n += 1
      x_true = np.ones(n,dtype=float)
      x_approx = specializedLU(n)
      x_err = la.norm(x_true - x_approx, np.inf)
      print("{:5d} : {:<6.4e}".format(n,x_err))
   print("The largest n with solution within six significant digits is {}".format(n-1))


def problem21():
   "problem set 2.1, problem 21, page 57"
   A=np.array([[1,-1,-1],[0,1,-2],[0,0,1]],dtype=float)
   invA = la.inv(A)
   norm_A_F = la.norm(A,'fro')
   norm_invA_F = la.norm(invA,'fro')
   norm_A_inf = la.norm(A,np.inf)
   norm_invA_inf = la.norm(invA,np.inf)
   print("""
Problem 21:
A =
{}
inverse(A) = 
{}
   Frobenius norm of A, ||A||_F = {}
   Frobenius norm of A^-1, ||A^-1||_F = {}
   Condition number using Frobenius norm, k(A)_F = {}
   inf norm of A, ||A||_inf = {}
   ing norm of A^-1, ||A^-1||_inf = {}
   Condition number using inf norm , k(A)_inf = {}
""".format(A,invA,norm_A_F,norm_invA_F, norm_A_F * norm_invA_F,
              norm_A_inf, norm_invA_inf, norm_A_inf*norm_invA_inf))


def problem25():
   "problem set 2.1, problem 25, pages 58-59"
   import LUdecomp
   m = np.array([10,4,5,6],dtype=float)
   u = np.array([0.25,0.3,0.2])
   g = 9.82  # m/s^2
   theta = np.pi/4  # radians = 45 degrees
   A = np.array( [[1,0,0,m[0]],
                  [-1,1,0,m[1]],
                  [0,-1,+1,m[2]],
                  [0,0,-1,m[3]]] ,dtype=float)
   A_orig = A.copy()
   st,ct = np.sin(theta),np.cos(theta)
   b = np.array([ m[0]*g*(st-u[0]*ct),
                  m[1]*g*(st-u[1]*ct),
                  m[2]*g*(st-u[2]*ct),
                  -m[3]*g ])
   b_orig = b.copy()
   LU = LUdecomp.LUdecomp(A)
   x = LUdecomp.LUsolve(LU,b)
   print("""
Problem 25
   T_1 = {:9.6f}
   T_2 = {:9.6f}
   T_3 = {:9.6f}
   a   = {:9.6f} 
with residual ||Ax-b||_2 = {:6.4e}""".format(
      x[0],x[1],x[2],x[3],la.norm(np.dot(A_orig,x)-b_orig),2)) 
 

if __name__ == "__main__" :
   print(__doc__)  # module doc from beginning of file
   problem8()
   problem11()
   problem15()
   problem21()
   problem25()


