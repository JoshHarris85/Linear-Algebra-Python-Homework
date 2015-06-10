""" 
by Josh Harris for CS 3513
Homework4.py 
Homework 4
Due 2/19/14

@author: JMH
"""

#Functions Ax and iterEqs taken from book and revised


import numpy as np
from gaussSeidel import *
from conjGrad import *

#Problem 17
def iterEqs(x, omega):
    n = len(x)
    x[0] = omega*(x[1] - x[n-1])/4.0 + (1.0 - omega)*x[0]
    for i in range(1,n-1):
        x[i] = omega*(x[i-1] + x[i+1])/4.0 + (1.0 - omega)*x[i]
    x[n-1] = omega*(100 - x[0] + x[n-2])/4.0 + (1.0 - omega)*x[n-1]
    return x


def problem17():   
    n = eval(input("Number of equations ==> "))
    x = np.zeros(n)
    x,numIter,omega = gaussSeidel(iterEqs,x)
    print("\nNumber of iterations =", numIter)
    print("\nRelaxation factor =", omega)
    print("\nThe solution for 17 is:\n", x)
    input("\nPress return to exit")
problem17()

#Problem 18
def Ax(v):
    n = len(v)
    Ax = np.zeros(n)
    Ax[0] = 4.0*v[0] - v[1]+ v[n-1]
    Ax[1:n-1] = -v[0:n-2] + 4.0*v[1:n-1] -v[2:n]
    Ax[n-1] = -v[n-2] + 4.0*v[n-1] + v[0]
    return Ax
    

def problem18():
    n = eval(input("Number of equations ==> "))
    b = np.zeros(n)
    b[n-1] = 100.00
    x = np.zeros(n)
    x, numIter = conjGrad(Ax,x,b)
    print("\nNumber of iterations =", numIter)
    print("\nThe solution for 18 is:\n", x)
    input("\nPress return to exit")
problem18()

#Problem 19
def Ax2(v):
    n = len(v)
    Ax2 = np.zeros(n)
    Ax2[0] = -4.0 * v[0] + v[1] + v[3]
    Ax2[1] = v[0] - 4.0 * v[1] + v[2] + v[4]
    Ax2[2] = v[1] - 4.0 * v[2] + v[5]
    Ax2[3] = v[0] - 4.0 * v[3] + v[4] + v[6]
    Ax2[4] = v[1] + v[3] - 4.0 * v[4] + v[5] + v[7]
    Ax2[5] = v[2] + v[4] - 4.0 * v[5] + v[8]
    Ax2[6] = v[3] - 4.0 * v[6] + v[7]
    Ax2[7] = v[4] + v[6] - 4.0 * v[7] + v[8]
    Ax2[8] = v[5] + v[7] - 4.0 * v[8]
    return Ax2
    
def problem19():   
    b = np.array([0.0, 0.0, 100.0, 0.0, 0.0, 100.0, 200.0, 200.0, 300.0],dtype=float)
    n = 9
    x = np.zeros(n)
    x, numIter = conjGrad(Ax2, x, b)
    print("\nNumber of iterations =",numIter)
    print("\nThe solution for 19 is:\n",x) 
    input("\nPress return to exit")
problem19()
#Creation of matrices for verification
n=20

A1 = np.array([[4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
               [-1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1],
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4]],dtype=float)               
B1 = np.zeros(n)
B1[n-1] = 100.0


A2 = np.array([[4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
               [-1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -1],
               [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 4]],dtype=float)             
B2 = np.zeros(n)
B2[n-1] = 100.0


A3 = np.array([[-4, 1, 0, 1, 0, 0, 0, 0, 0],
               [1, -4, 1, 0, 1, 0, 0, 0, 0],
               [0, 1, -4, 0, 0, 1, 0, 0, 0],
               [1, 0, 0, -4, 1, 0, 1, 0, 0],
               [0, 1, 0, 1, -4, 1, 0, 1, 0],
               [0, 0, 1, 0, 1, -4, 0, 0, 1],
               [0, 0, 0, 1, 0, 0, -4, 1, 0],
               [0, 0, 0, 0, 1, 0, 1, -4, 1],
               [0, 0, 0, 0, 0, 1, 0, 1, -4]],dtype=float)
B3 = np.array([0.0, 0.0, 100.0, 0.0, 0.0, 100.0, 200.0, 200.0, 300.0],dtype=float)



# Class code taken from iteration article
class AxClass:
    def __init__(self,matrixA):
        self.A = matrixA.copy()
    def __call__(self,x):
        return np.dot(self.A,x)

# Class code taken from iteration article
class exClass:
    def __init__(self,matrixA,vectorB):
        self.A = matrixA.copy()
        self.b = vectorB.copy()
    def __call__(self,x,w):
        for i in range(self.A.shape[0]):
            x[i] = x[i] +(w/self.A[i,i])*(self.b[i]- np.dot(self.A[i,:],x))
        return x

def problem17Verify(A1, B1):
    exC = exClass(A1, B1)
    x = B1.copy()
    x, i, w = gaussSeidel(exC,x)
    print("\nVerification for problem 17:")
    print("\nNumber of iterations =", i)
    print("\nRelaxation factor =", w)
    print("\nThe solution is:\n", x)
problem17Verify(A1, B1)
print("\nThe error between 8.65000954e-14 and -1.88823725e-13 is so small it is essentially zero.")
print ("8.65000954e-14 - (-1.88823725e-13) =",8.65000954e-14 - (-1.88823725e-13), "= 0.00000000000027532382")

def problem18Verify(A2, B2):
    ex2 = AxClass(A2)
    x = B2.copy()
    x, i = conjGrad(ex2,x,B2)
    print("\nVerification for problem 18:")
    print("\nNumber of iterations =", i)
    print("\nThe solution is:\n", x)
problem18Verify(A2,B2)


def problem19Verify(A3, B3):
    ex2 = AxClass(A3)
    x = B3.copy()
    x, i = conjGrad(ex2,x,B3)
    print("\nVerification for problem 19:")
    print("\nNumber of iterations =", i)
    print("\nThe solution is:\n", x)
problem19Verify(A3,B3)
