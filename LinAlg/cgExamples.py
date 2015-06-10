"""Example Conjugate Gradient calculations"""

# Doug Heisterkamp
# Last modified: 2/10/2014
import numpy as np
import conjGrad

# first example is for Ax=b with
A=np.array([[2,-1, 1],[-1,3,-2],[1,-2,4]],dtype=float)
b=np.array([10,-8,16],dtype=float)

def ex1(x):
   assert(len(x)==3)
   # A is only 3x3, write out explicitly
   Ax = np.zeros(3)
   Ax[0] =  2*x[0] -x[1] + x[2]
   Ax[1] =  -x[0] +3*x[1] -2* x[2]
   Ax[2] =  x[0] -2*x[1] + 4*x[2]
   return Ax

# initial x = b
x = b.copy()
x,i = conjGrad.conjGrad(ex1,x,b)

print("""Example 1, conjGrad results 
     x={}, 
     iterations = {} """.format(x,i))

# second example, use a functor class for dense A
# a class with a "call" method so it acts like a function
class AxClass:
   def __init__(self,matrixA):
      self.A = matrixA.copy()
   def __call__(self,x):
      return np.dot(self.A,x)

# now use the class definition to create 
# an instance with A 
ex2 = AxClass(A)
x = b.copy()
x,i = conjGrad.conjGrad(ex2,x,b)

print("""Example 2, conjGrad results 
     x={}, 
     iterations = {}""".format(x,i))


