"""Example SOR calculations"""

# Doug Heisterkamp
# Last modified: 2/10/2014
import numpy as np
import gaussSeidel

# first example is for Ax=b with
A=np.array([[2,-1, 1],[-1,3,-2],[1,-2,4]],dtype=float)
b=np.array([10,-8,16],dtype=float)

def ex1(x,w):
   # A is only 3x3, write out explicitly
   x[0] = (1-w)*x[0] + (w/2)*(10 + x[1] - x[2])
   x[1] = (1-w)*x[1] + (w/3)*(-8 + x[0] + 2*x[2])
   x[2] = (1-w)*x[2] + (w/4)*(16 - x[0] + 2*x[1])
   return x

# initial x = b
x = b.copy()
x,i,w = gaussSeidel.gaussSeidel(ex1,x)

print("""Example 1, SOR results 
     x={}, 
     iterations = {}, 
     w={}""".format(x,i,w))

# second example, use a functor class for dense A,b
# a class with a "call" method so it acts like a function
class exClass:
   def __init__(self,matrixA,vectorB):
      self.A = matrixA.copy()
      self.b = vectorB.copy()
   def __call__(self,x,w):
      # using residual form 
      for i in range(self.A.shape[0]):
         x[i] = x[i] +(w/self.A[i,i])*(self.b[i] 
                        - np.dot(self.A[i,:],x))
      return x

# now use the class definition to create 
# an instance with A and b
exC = exClass(A,b)
x = b.copy()
x,i,w = gaussSeidel.gaussSeidel(exC,x)

print("""Example 2, SOR results 
     x={}, 
     iterations = {}, 
     w={}""".format(x,i,w))


