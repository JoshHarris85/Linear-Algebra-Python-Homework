"""
Using a couple of problems from HW2 to test LUpivotComplete.
"""

# Doug Heisterkamp
# Last modified: 2/3/2014

import numpy as np
import numpy.linalg as la

def problem8():
   "problem set 2.1, problem 8, page 56"
   import LUpivotComplete
   A = np.array([[-3,6,4],[9,-8,24],[-12,24,-26]],dtype=float)
   A_orig = A.copy()
   LU,sr,sc = LUpivotComplete.LUdecomp(A)
   b = np.array([-3,65,-42],dtype=float)
   b_orig = b.copy()
   x = LUpivotComplete.LUsolve(LU,b,sr,sc)
   # extract L and U for verification
   U = np.triu(LU)  # 
   L = np.tril(LU)
   L[ np.diag_indices_from(L) ] = 1.0 
   print("""
Problem 8:
A = 
{}
Row permutation P = {}
Col permutation Q = {}
PAQ = 
{}
LU decomposition PAQ= LU, LU (in one matrix) = 
{}
Solving Ax=b, with b = {}
Solution x = {}
Verifying solution: 
     residual ||Ax-b||_2 = {}
     ||PAQ - dot(L,U)||_inf = {}
""".format(A_orig,sr,sc, A_orig[sr][:,sc],LU,b_orig,x, 
   la.norm(np.dot(A_orig,x)-b_orig,2), 
   la.norm(A_orig[sr][:,sc] - np.dot(L,U),np.inf))
   )

def problem25():
   "problem set 2.1, problem 25, pages 58-59"
   import LUpivotComplete
   m = np.array([10,4,5,6],dtype=float)
   u = np.array([0.25,0.3,0.2])
   g = 9.82  # m/s^2
   theta = 45  # degrees
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
   LU,sr,sc = LUpivotComplete.LUdecomp(A)
   x = LUpivotComplete.LUsolve(LU,b,sr,sc)
   # extract L and U for verification
   U = np.triu(LU)  # 
   L = np.tril(LU)
   L[ np.diag_indices_from(L) ] = 1.0 
   print("""
Problem 25
A = 
{}
Row permutation P = {}
Col permutation Q = {}
PAQ = 
{}
Solving Ax=b, with b = {}
Solution x = {}
Verifying solution: 
     residual ||Ax-b||_2 = {}
     ||PAQ - dot(L,U)||_inf = {}
Unpacking x into variable names from problem:
   T_1 = {:9.6f}
   T_2 = {:9.6f}
   T_3 = {:9.6f}
   a   = {:9.6f} 
""".format(
      A_orig,sr,sc, A_orig[sr][:,sc], b_orig, x,
      la.norm(np.dot(A_orig,x)-b_orig,2),
      la.norm(A_orig[sr][:,sc] - np.dot(L,U),np.inf),
      x[0],x[1],x[2],x[3])) 
 

if __name__ == "__main__" :
   print(__doc__)  # module doc from beginning of file
   problem8()
   problem25()


