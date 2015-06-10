""" 
by Josh Harris for CS 3513
Homework5.py 
Homework 5
Due 2/26/14
Created on Sun Feb 23 18:55:38 2014
@author: JMH
"""

import numpy as np
import scipy.interpolate as ip
import neville
import cubicSpline 
import rational
import matplotlib.pyplot as plt

def lagrangeA(xArr,yArr):
    y = 0
    L0 = ((y - yArr[1]) * (y - yArr[2])) / ((yArr[0] - yArr[1]) * (yArr[0] - yArr[2]))
    L1 = ((y - yArr[0]) * (y - yArr[2])) / ((yArr[1] - yArr[0]) * (yArr[1] - yArr[2]))
    L2 = ((y - yArr[0]) * (y - yArr[1])) / ((yArr[2] - yArr[0]) * (yArr[2] - yArr[1]))
    x = (xArr[0] * L0) + (xArr[1] * L1) + (xArr[2] * L2)
    return x
    
def lagrangeB(xArr,yArr):
    y = 0
    L0 = ((y - yArr[1]) * (y - yArr[2]) * (y - yArr[3])) / ((yArr[0] - yArr[1]) * (yArr[0] - yArr[2]) * (yArr[0] - yArr[3]))
    L1 = ((y - yArr[0]) * (y - yArr[2]) * (y - yArr[3])) / ((yArr[1] - yArr[0]) * (yArr[1] - yArr[2]) * (yArr[1] - yArr[3]))
    L2 = ((y - yArr[0]) * (y - yArr[1]) * (y - yArr[3])) / ((yArr[2] - yArr[0]) * (yArr[2] - yArr[1]) * (yArr[2] - yArr[3]))
    L3 = ((y - yArr[0]) * (y - yArr[1]) * (y - yArr[2])) / ((yArr[3] - yArr[0]) * (yArr[3] - yArr[1]) * (yArr[3] - yArr[2]))
    x = (xArr[0] * L0) + (xArr[1] * L1) + (xArr[2] * L2) + (xArr[3] * L3)
    return x 

def problem2A():
    xArr = np.array([2, 2.5, 3])
    yArr = np.array([0.8509, -0.4112, -1.5727])
    print("The answer for 2 Part A is:", lagrangeA(xArr,yArr))
    
def problem2B():
    xArr = np.array([1.5, 2, 2.5, 3])
    yArr = np.array([1.9047, 0.8509, -0.4112, -1.5727])
    print("The answer for 2 Part B is:", lagrangeB(xArr,yArr))    
problem2A()
problem2B()


def problem3():
    x = 0.7692
    xArr = np.array([0, 0.5, 1, 1.5])
    yArr = np.array([1.8421, 2.4694, 2.4921, 1.9047])
    Answer = neville.neville(xArr,yArr,x)
    print("\nThe answer for 3 is:", Answer)
problem3()


def problem9():
    x = np.array([0, 3, 6])
    y = np.array([1.225, 0.905, 0.652])
    p2 = ip.lagrange(x,y)
    print("\nThe answer to problem 9 is: ", p2)  
problem9()
    
    
def problem12():
    yArr = np.array([1.0, 0.8, 0.6, 0.4, 0.2])
    xArr = np.array([-1.049, -0.266, 0.377, 0.855, 1.150])
    k = cubicSpline.curvatures(xArr, yArr)
    print ("\nThe answer to problem 12 is:", cubicSpline.evalSpline(xArr, yArr, k, 0.0))
problem12()


def problem15():
    xArr = np.array([-250, -200, -100, 0, 100, 300])
    yArr = np.array([0.0163, 0.318, 0.699, 0.870, 0.941, 1.04])
    #Code taken from pg 118 in the book.
    x = np.arange(-250, 500)
    n = len(x)
    y = np.zeros((n,2))
    for i in range(n):
        y[i,0] = rational.rational(xArr,yArr,x[i])
        y[i,1] = neville.neville(xArr,yArr,x[i])
    plt.plot(xArr, yArr, 'o', x, y[:,0], '-', x, y[:,1],'--')
    plt.xlabel('x')
    plt.legend(('Data', 'Rational', 'Neville'), loc = 0)
    plt.show()
    print("\nOn problem 15, the graph shows the rational interpolent is smoother around the 100 to 300 x range, meaning the rational is a superior interpolant")
problem15()


def problem17():
    re = np.array([0.2, 2, 20, 200, 2000, 20000])
    cD = np.array([103, 13.9, 2.72, 0.800, 0.401, 0.433])
    x = np.log(re)
    y = np.log(cD)
    k = cubicSpline.curvatures(x,y)
    x1 = cubicSpline.evalSpline(x,y,k,np.log(5))
    x2 = cubicSpline.evalSpline(x,y,k,np.log(50))
    x3 = cubicSpline.evalSpline(x,y,k,np.log(500))
    x4 = cubicSpline.evalSpline(x,y,k,np.log(5000))
    x1Resto = np.exp(x1)
    x2Resto = np.exp(x2)
    x3Resto = np.exp(x3)
    x4Resto = np.exp(x4)
    print("\nThe answers to problem 17 are:")
    print("Re at 5 is", x1Resto)
    print("Re at 50 is", x2Resto)
    print("Re at 500 is", x3Resto)
    print("Re at 5000 is", x4Resto)
problem17()