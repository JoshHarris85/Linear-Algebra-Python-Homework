"""
CS 3513, Homework 3
One step of iterative refinement
Example Solution
"""

# Doug Heisterkamp
# Last modified: 2/9/2014

import numpy as np
import numpy.linalg as la
import LUdecomp
import LUpivot
import LUpivotComplete

def condm(a):
   return la.norm(a,np.inf)*la.norm(la.inv(a),np.inf)

def calcError(x_tuple,x_true):
   return [ la.norm(i-x_true,2) for i in x_tuple]

def refineLUdecomp(A,b):
   A_work = A.copy()
   b_work = b.copy()
   LU = LUdecomp.LUdecomp(A_work)
   x_orig = LUdecomp.LUsolve(LU,b_work)
   r = b - np.dot(A,x_orig)
   z = LUdecomp.LUsolve(LU,r)
   x_ref = x_orig + z
   return x_orig,x_ref

def refineLUpivot(A,b):
   A_work = A.copy()
   b_work = b.copy()
   LU,seq = LUpivot.LUdecomp(A_work,1e-14)
   x_orig = LUpivot.LUsolve(LU,b_work,seq)
   r = b - np.dot(A,x_orig)
   z = LUpivot.LUsolve(LU,r,seq)
   x_ref = x_orig + z
   return x_orig,x_ref


def refineLUpivotComplete(A,b):
   A_work = A.copy()
   b_work = b.copy()
   LU,seqR,seqC = LUpivotComplete.LUdecomp(A_work,1e-14)
   x_orig = LUpivotComplete.LUsolve(LU,b_work,seqR,seqC)
   r = b - np.dot(A,x_orig)
   z = LUpivotComplete.LUsolve(LU,r,seqR,seqC)
   x_ref = x_orig + z
   return x_orig,x_ref

def hilbert(n, dtype=float) :
   return np.fromfunction( lambda i,j :1/( i+j+1), (n,n),dtype=dtype)

def M(n, dtype=float) :
   a = np.zeros((n,n),dtype=dtype)
   a[np.tril_indices_from(a,-1)] = -1.0
   a[np.diag_indices_from(a)] = 1.0
   a[:,-1] = 1.0
   return a


def iterRefExer():
   results = np.zeros((7,17),dtype=float)
   for n in range(3,10):
      r = n-3
      results[r,0] = n
      h = hilbert(n)
      x_true_h = np.ones(n)
      bh = np.sum(h,axis=1) 
      m = M(10*n)
      x_true_m = np.ones(m.shape[1])
      bm = np.sum(m,axis=1) 
      results[r,[1,2]] = calcError(refineLUdecomp(h,bh),x_true_h)
      results[r,[3,4]] = calcError(refineLUpivot(h,bh),x_true_h)
      results[r,[5,6]] = calcError(refineLUpivotComplete(h,bh),x_true_h)
      results[r,[7,8]] = calcError(refineLUdecomp(m,bm),x_true_m)
      results[r,[9,10]] = calcError(refineLUpivot(m,bm),x_true_m)
      results[r,[11,12]] = calcError(refineLUpivotComplete(m,bm),x_true_m)
      results[r,13] = la.norm((la.solve(h,bh) - x_true_h),2)
      results[r,14] = la.norm((la.solve(m,bm) - x_true_m),2)
      results[r,[15,16]] = [ condm(h), condm(m)]
   #print("results = \n",results)
   return results
       
def pe(e1,e2):
   if e2 != 0:
      print("{:>10.2e}".format(e1/e2),end="")
   else:
      print("{:>10s}".format('--------'),end="")
          
def pe2(e1,e2):
   print("{:>9.2e}/{:<9.2e}".format(e1,e2),end="")

def printRes(r,label,c1,c2):
   print("{:<12s}: ".format(label),end="")
   for i in range(r.shape[0]):
      pe(r[i,c1],r[i,c2])
   print()


def printResE(r,label,c1,c2):
   print("{:<12s}: ".format(label),end="")
   for i in range(r.shape[0]):
      pe2(r[i,c1],r[i,c2])
   print()


def printCol(r,label,c,fmt):
   print("{:<12s}: ".format(label),end="")
   for i in range(r.shape[0]):
      print(fmt.format(r[i,c]),end="")
   print()



if __name__ == "__main__" :
   print(__doc__)  # module doc from beginning of file
   r = iterRefExer()

   
   print(" Condition numbers : ||A||_2 ||inv(A)||_2")
   printCol(r,"n",0,"{:>10.0f}")
   printCol(r,"cond h",15,"{:>10.2e}")
   printCol(r,"cond m",16,"{:>10.2e}")
   print()
   print("Errors: ||x - hat(x)||_2 / ||x - x_refined||_2 ")
   printCol(r,"n",0,"{:>18.0f}")
   printResE(r,"decomp h",1,2)
   printResE(r,"decomp m",7,8)
   printResE(r,"pivot h",3,4)
   printResE(r,"pivot m",9,10)
   printResE(r,"comp h",5,6)
   printResE(r,"comp m",11,12)
   printCol(r,"solve h",13,"{:>18.2e}")
   printCol(r,"solve m",14,"{:>18.2e}")


   print()
   print("Error Ratios : ||x - hat(x)||_2 / ||x - x_refined||_2 ")
   printCol(r,"n",0,"{:>10.0f}")
   printRes(r,"decomp h",1,2)
   printRes(r,"decomp m",7,8)
   printRes(r,"pivot h",3,4)
   printRes(r,"pivot m",9,10)
   printRes(r,"comp h",5,6)
   printRes(r,"comp m",11,12)

   print("\nNote: M matrices used 10*n for dimensions")





