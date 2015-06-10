""" 
by Josh Harris for CS 3513
Homework3.py 
Homework 2
Due 2/12/14

@author: JMH
"""


import numpy as np
import numpy.linalg as la
import math
import LUpivot

#Problem 16              
def problem16():
    s = math.degrees(math.sin(53))
    c = math.degrees(math.cos(53))
    A = np.array([[c,1,0,0,0],[0,s,0,0,1],[0,0,2*s,0,0],[0,-c,c,1,0],[0,s,s,0,0]],dtype=float)
    B = np.array([[0],[0],[1],[0],[0]],dtype=float)

    x, seq = LUpivot.LUdecomp(A)
    print("The answer to problem 16 is:") 

    print(LUpivot.LUsolve(x,B,seq))
    print ("\n")
problem16()

#problem 17
def problem17():
    R = 5 

    while(R <= 20):
        A = np.array([[50+R,-R,-30],[-R,65+R,-15],[-30,-15,45]],dtype=float)
        B = np.array([[0],[0],[120]],dtype=float)
        x, seq = LUpivot.LUdecomp(A)
        print("The answer to problem 17 with Omega at",R)
        print(LUpivot.LUsolve(x,B,seq))
        R = R * 2
        
problem17()

#Problem 20
def problem20(): 
    A = np.array([[-8,4,0,0,0],[8,-10,2,0,0],[0,6,-11,5,0],[0,0,3,-7,4],[0,0,0,2,-4]],dtype=float)
    A_orig = A.copy()
    B = np.array([[-80],[0],[0],[0],[-30]],dtype=float)
    b_orig = B.copy()
    LU, seq = LUpivot.LUdecomp(A)
    x = LUpivot.LUsolve(A,B,seq)
    U = np.triu(LU)  # 
    L = np.tril(LU)
    L[ np.diag_indices_from(L) ] = 1.0 
    #Code taken from Dr. Heisterkamp's hw2.py for verifying a solution.
    print("""
Problem 20:
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

    
problem20()


#Final problem
#various code below taken from hw2.py
""" I unfortunately could not figure out how to complete this problem. 
I went ahead and tried my best to implement the algorithm for iterative refinement, and 
I tried to create a Wilkinson matrix but failed. 
"""
def hilbert(n, dtype=float) :
      return np.fromfunction( lambda i,j :1/( i+j+1), (n,n),dtype=dtype)

#Failed attempt to make a Wilkinson Matrix    
def createM(n):
    Mn = np.array([[]])
    for i in range(n):
        for j in range(n):
            if(i == j):
                np.append(Mn,1)
            if(i > j):
                np.append(Mn,-1)
            if(j == n-1):
                np.append(Mn,1)
            else:
                np.append(Mn,0)
    return Mn
Mn = createM(5)
print (Mn)            
       
      
def iterativeRefinement():
    errRef = 0
    n = 3
    while n <= 9:
        Hm = hilbert(n)
        b = np.sum(Hm,axis=1) 
        LU, seq = LUpivot.LUdecomp(Hm)
        x = LUpivot.LUsolve(Hm,b,seq)
        #find the residual
        residual = np.dot(Hm,x)-b,2
        #Solve the system using the residual as B values
        z = LUpivot.LUsolve(Hm,residual,seq)
        #Update the correction
        #Note show err = || x - xhat || euclidean 2 
        errOrig = abs(b) - abs(LU)
        errRef = LU  + z
        errRatio = errOrig / errRef
        print(errRatio)
        n = n + 1
        
        



