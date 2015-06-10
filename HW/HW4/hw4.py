"""
CS 3513, Homework 4
Problem Set 2.3
Example Solution
"""

# Doug Heisterkamp
# Last modified: 2/17/2014

import numpy as np
import numpy.linalg as la
import gaussSeidel
import conjGrad


def genAexer_2_17(n):
   a = np.zeros((n,n),float)
   t  = np.diag_indices_from(a)
   a[t] = 4.0
   # shift diagonal indices for off diagonal indices 
   upper = (t[0][:-1],(t[1]+1)[:-1])
   lower = (t[0][1:],(t[1]-1)[1:])
   a[upper] = -1
   a[lower] = -1
   a[0,-1] = 1
   a[-1,0] = 1
   return a


def genAp19():
   a = np.zeros((9,9),float)
   t  = np.diag_indices_from(a)
   a[t] = -4.0
   # shift diagonal indices for off diagonal indices 
   upper = (t[0][:-1],(t[1]+1)[:-1])
   a[upper] = [1,1,0,1,1,0,1,1]
   lower = (t[0][1:],(t[1]-1)[1:])
   a[lower] = [1,1,0,1,1,0,1,1]
   upper = (t[0][:-3],(t[1]+3)[:-3])
   a[upper] = 1
   lower = (t[0][3:],(t[1]-3)[3:])
   a[lower] = 1
   return a


class gsClass:
   "functor class for gaussSiedel"
   def __init__(self,matrixA,vectorB):
      self.A = matrixA.copy()
      self.b = vectorB.copy()
   def __call__(self,x,w):
      # using residual form 
      for i in range(self.A.shape[0]):
         x[i] = x[i] +(w/self.A[i,i])*(self.b[i] 
                        - np.dot(self.A[i,:],x))
      return x

class cgClass:
   "functor class for conjGrad"
   def __init__(self,matrixA):
      self.A = matrixA.copy()
   def __call__(self,x):
      return np.dot(self.A,x)



def problem17():
   "problem set 2.3, problem 17, page 100"
   def p17eqn(x,w):
      n = len(x)
      y = 1.0-w
      x[0] = w*(x[1] - x[n-1])/4.0 +y*x[0]
      for i in range(1,n-1):
         x[i] = w*(x[i-1] + x[i+1])/4.0 + y*x[i]
      x[n-1] = w*(100.0 - x[0] +x[n-2])/4.0 + y*x[n-1]
      return x
   n  = 20
   b = np.zeros(n,float)
   b[-1] = 100.0
   x = b.copy()
   x,i,w = gaussSeidel.gaussSeidel(p17eqn,x)
   a = genAexer_2_17(n)
   # print("p17, a=\n",a)
   exC = gsClass(a,b)
   altX = b.copy()
   altX,altI,altW = gaussSeidel.gaussSeidel(exC,altX)
   print("""Problem 17, SOR results 
     x={}, 
     iterations = {}, 
     w={}
     ||x - altX||_2 = {}""".format(x,i,w,la.norm(x-altX)))



def problem18():
   "problem set 2.3, problem 18, page 101"
   def p18Ax(x):
      n = len(x)
      res = np.zeros(n)
      res[0] = 4*x[0] -x[1] +x[-1]
      for i in range(1,n-1):
         res[i] = -x[i-1] + 4*x[i] -x[i+1]
      res[n-1] = +x[0] - x[n-2] + 4 * x[n-1]
      return res
   n  =20 
   b = np.zeros(n,float)
   b[-1] = 100.0
   x = b.copy()
   x,i = conjGrad.conjGrad(p18Ax,x,b)
   a = genAexer_2_17(n)
   # print("p18, a=\n",a)
   exC = cgClass(a)
   altX = b.copy()
   altX,altI = conjGrad.conjGrad(exC,altX,b)
   print("""Problem 18, Conjugate Gradient results 
     x={}, 
     iterations = {}, 
     ||x - altX||_2 = {}""".format(x,i,la.norm(x-altX)))

def problem19():
   "problem set 2.3, problem 19, page 101"
   def p19Ax(x):
      assert(len(x)==9)
      res = np.zeros(9)
      res[0] =-4*x[0] +x[1]       +x[3]
      res[1] =-4*x[1] +x[2] +x[0] +x[4]
      res[2] =-4*x[2]       +x[1] +x[5]
      res[3] =-4*x[3] +x[4]       +x[6] + x[0]
      res[4] =-4*x[4] +x[5] +x[3] +x[7] + x[1]
      res[5] =-4*x[5]       +x[4] +x[8] + x[2]
      res[6] =-4*x[6] +x[7]             + x[3]
      res[7] =-4*x[7] +x[8] +x[6]       + x[4]
      res[8] =-4*x[8]       +x[7]       + x[5]
      return res
   b = np.array([0,0,100,0,0,100,200,200,300],float)
   T = b.copy()
   T,i = conjGrad.conjGrad(p19Ax,T,b)
   a = genAp19()
   # print("p19, a=\n",a)
   exC = cgClass(a)
   altT = b.copy()
   altT,altI = conjGrad.conjGrad(exC,altT,b)
   print("""Problem 19, Conjugate Gradient results 
     T={}, 
     iterations = {}, 
     ||T - altT||_2 = {}""".format(T,i,la.norm(T-altT)))


if __name__ == "__main__" :
   print(__doc__)  # module doc from beginning of file
   problem17()
   print()
   problem18()
   print()
   problem19()





