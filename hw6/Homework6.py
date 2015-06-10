""" 
by Josh Harris for CS 3513
Homework6.py 
Homework 6
Due 3/5/14
@author: JMH
"""

import numpy as np
from polyFit import * 
import choleski as cho


#Problem 0 
def problem0():
    x = np.array([-0.04, 0.93, 1.95, 2.90, 3.83, 5.00, 5.98, 7.05, 8.21, 9.08, 10.09]) 
    y = np.array([-8.66, -6.44, -4.36, -3.27, -0.88, 0.87, 3.31, 4.63, 6.19, 7.40, 8.85])
    m = 2
    Answer = polyFit(x,y,m)
    print("The answer to problem 0 using polyfit is: ", Answer)
problem0()

#Problem 1 using Cholesky's decomp
def problem1():
    x = np.array([-0.04, 0.93, 1.95, 2.90, 3.83, 5.00, 5.98, 7.05, 8.21, 9.08, 10.09]) 
    y = np.array([-8.66, -6.44, -4.36, -3.27, -0.88, 0.87, 3.31, 4.63, 6.19, 7.40, 8.85])
    m = 2
    #Taken from polyFit to form matrix A and B 
    a = np.zeros((m+1,m+1))
    b = np.zeros(m+1)
    s = np.zeros(2*m+1)
    for i in range(len(x)):
        temp = y[i]
        for j in range(m+1):
            b[j] = b[j] + temp
            temp = temp*x[i]
        temp = 1.0
        for j in range(2*m+1):
            s[j] = s[j] + temp
            temp = temp*x[i]
    for i in range(m+1):
        for j in range(m+1):
            a[i,j] = s[i+j]
    print("\nFormed A:")
    print(a)
    print("\nFormed B:")
    print(b)
    
    decomp = cho.choleski(a)
    Answer = cho.choleskiSol(decomp, b)
    print("\nThe answer to Problem 1 using Choleski is: ", Answer)
problem1()

def backSub(a, b): 
    n = len(b)
    for k in range(n-1,-1,-1):
        b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
    return b

#Problem 2 using QR decomp
def problem2():
    a = np.array([[1.10000000e+01, 5.49800000e+01, 3.89869400e+02],
                  [5.49800000e+01, 3.89869400e+02, 3.10728457e+03],
                  [3.89869400e+02, 3.10728457e+03, 2.63808673e+04]])
    b = np.array([7.64, 237.0955, 2236.003879])
    q,r = np.linalg.qr(a, mode='complete')
    transposeQ = np.transpose(q)
    y = np.dot(transposeQ, b)
    Answer = backSub(r, y)
    print("\nThe answer to problem 2 using QR backsub is:", Answer)
problem2()



#Problem 3 using SVD decomp 
def problem3():
    a = np.array([[1.10000000e+01, 5.49800000e+01, 3.89869400e+02],
                  [5.49800000e+01, 3.89869400e+02, 3.10728457e+03],
                  [3.89869400e+02, 3.10728457e+03, 2.63808673e+04]])
    b = np.array([7.64, 237.0955, 2236.003879])
    U, s, V = np.linalg.svd(a)
    uTranspose = np.transpose(U)
    vTranspose = np.transpose(V)
    #code source from http://sukhbinder.wordpress.com/2013/03/26/solving-axb-by-svd/
    c = np.dot(uTranspose,b)
    w = np.linalg.solve(np.diag(s),c)
    x = np.dot(vTranspose,w)
    print("\nThe answer to problem 3 using SVD is:" ,x)
problem3()





#Problem 4 using numpy.linalg.lstsq
def problem4():
    a = np.array([[1.10000000e+01, 5.49800000e+01, 3.89869400e+02],
                  [5.49800000e+01, 3.89869400e+02, 3.10728457e+03],
                  [3.89869400e+02, 3.10728457e+03, 2.63808673e+04]])
    b = np.array([7.64, 237.0955, 2236.003879])
    Answer = np.linalg.lstsq(a, b)[0]
    print("\nThe answer to problem 4 using numpy's lstsq is: ",Answer)
problem4()
    

    

    