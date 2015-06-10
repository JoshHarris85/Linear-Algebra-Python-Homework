""" 
Homework2.py 
by Josh Harris for CS 3513
Homework 2
Due 2/5/14

Created on Mon Jan 27 22:57:19 2014

@author: JMH
"""
#Problem 8
import numpy as np
import error
import math
import numpy.linalg as LA
import scipy
#From discussion board OSU d2l
from scipy.linalg import hilbert

a1 = np.matrix([[-3,  6,  -4],
              [9,  -8,  24],
              [-12,  24,  -26]])
              
b1 = np.matrix([[-3], 
              [65],
              [-42]])

#Code Taken From Numerical Methods in Engineering with Python 3
def doolittleDecomp(a):
    n = len(a)
    for k in range(0, n-1):
        for i in range(k+1, n):
            if a[i,k] != 0.0:
                lam = a[i, k]/a[k, k]
                a[i, k+1:n] = a[i,k+1:n] - lam*a[k, k+1:n]
                a[i, k] = lam
    return a
                
                
#Code Taken From Numerical Methods in Engineering with Python 3
def doolittlesolve(a,b):
    n = len(a)
    for k in range(1,n):
        b[k] = b[k] - np.dot(a[k, 0:k],b[0:k])
    b[n-1] = b[n-1]/a[n-1, n-1]    
    for k in range(n-2, -1, -1):
       b[k] = (b[k] - np.dot(a[k, k+1:n],b[k+1:n]))/a[k, k]
    return b
        
doolittleDecomp(a1)
print ("The answer for problem 8 is:")
print (doolittlesolve(a1,b1))
print("\n")







#NOTE: The type is set to float and there are no errors in the numbers of the array.
# Yet, I still get Matrix is not Positive Definite error. I saw someone else asked this on the discussion,
# and I did exactly as Dr. Heisterkamp informed the person asking without luck. 
#So I decided to complete the problem through cho_factor and cho_solve
#Problem 11
a2 = np.matrix([[1, 1, 1],  
              [1, 2, 2],  
              [1, 2, 3]], dtype=np.float)
              
b2 = np.matrix([[1], 
              [3/2],
              [3]], dtype=np.float)             

a22 = scipy.linalg.cho_factor(a2)
print("The answer to problem 11 is: ")
print(scipy.linalg.cho_solve(a22, b2))
print("\n")

"""
#Code Taken From Numerical Methods in Engineering with Python 3             
def choleski(a):
    n = len(a)
    for k in range(n):
        try:
            a[k,k] = math.sqrt(a[k,k] - np.dot(a[k,0:k],a[k,0:k]))
        except ValueError:
            error.err('Matrix is not positive definite')
        for i in range(k+1,n):
            a[i,k] = (a[i,k] - np.dot(a[i,0:k],a[k,0:k]))/a[k,k]
    for k in range(1,n): a[0:k,k] = 0.0
    return a
#Code Taken From Numerical Methods in Engineering with Python 3   
def choleskiSol(L,b):
    n = len(b)
  # Solution of [L]{y} = {b}  
    for k in range(n):
        b[k] = (b[k] - np.dot(L[k,0:k],b[0:k]))/L[k,k]
  # Solution of [L_transpose]{x} = {y}      
    for k in range(n-1,-1,-1):
        b[k] = (b[k] - np.dot(L[k+1:n,k],b[k+1:n]))/L[k,k]
    return b

choleski(a2)
print ("The answer for problem 11 is:")
print (choleskiSol(a2,b2))
"""





#Problem 15
#Code is not working correctly. I did not have enough time to finish it.
a3 = hilbert(2) 
b3 = []
x = []
i = 0
actNorm = 1*10**-6

while True:
    b3.append(a3[i].sum())
    x.append([1])
    xNorm = np.linalg.norm(x, np.inf)
    xApprox = np.linalg.solve(a3, b3)
    xHat = np.linalg.norm(xApprox, np.inf)
    approxErr = xNorm - xHat
    if abs(approxErr) < actNorm:
        break
    i = i + 1

print("The answer for 15 is: ")
print(approxErr)
print("\n")






#Problem 21
a4 = np.matrix([[1,  -1,  -1],
              [0,  1,  -2],
              [0,  0,  1]])


print("The answers for 21 are: ")
print("The condition of the euclidean norm is")
print(np.linalg.cond(a4))
print("The condition of the infinity norm is")
print(np.linalg.cond(a4, np.inf))
print("\n")


    
    
    
    

#Problem 25
m = [10, 4, 5, 6]

u = [0.25, 0.3, 0.2]

thet = 45

x = np.array([[1, 0, 0, m[0]], [-1, 1, 0, m[1]], [0, -1, 1, m[2]], [0, 0, -1, m[3]]])

y = np.array([m[0] * 9.82 * (math.sin(thet) - u[0] * math.cos(thet)), m[1] * 9.82 * (math.sin(thet) - u[1] * math.cos(thet)), m[2] * 9.82 * (math.sin(thet) - u[2] * math.cos(thet)), -m[3] * 9.82])

answer = np.linalg.solve(x, y)

print("The answers for 25 are: ")
print(answer) 





























    