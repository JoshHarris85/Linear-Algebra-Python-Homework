"""
CS 3513, Homework 3
Problem Set 2.2
Example Solution
"""

# Doug Heisterkamp
# Last modified: 2/8/2014

import numpy as np
import numpy.linalg as la
import LUpivot

def problem16():
   "problem set 2.2, problem 16, page 81"
   def symTruss(theta):
      "input : theta in degree, returns forces P_1 to P_5 as an array"
      r = theta*np.pi/180 
      c = np.cos(r)
      s = np.sin(r)
      A = np.array([[c,1,0,0,0],
                    [0,s,0,0,1],
                    [0,0,2*s,0,0],
                    [0,-c,c,1,0],
                    [0,s,s,0,0]],dtype=float)
      b = np.array([0,0,1,0,0],dtype=float)
      LU,seq = LUpivot.LUdecomp(A)
      x = LUpivot.LUsolve(LU,b,seq)
      return x
   
   x = symTruss(53)
   print("Problem 16: using theta = 53.")
   print("  P1={:5.4f}, P2={:5.4f}, P3={:5.4f}, P4={:5.4f}, P5={:5.4f}".format( 
      x[0], x[1], x[2], x[3], x[4])) 


def problem17():
   "problem set 2.2, problem 17, page 82"

   print("Problem 17, assuming typo in last equation.")
   for R in [5,10,20] :

      A = np.array([[50+R,-R,-30],
                    [-R,65+R,-15],
                    [-30,-15,45],  # typo in book? Using i1 for first i2
                   ],dtype=float)
      b = np.array([0,0,120],dtype=float)
      LU,seq = LUpivot.LUdecomp(A)
      i = LUpivot.LUsolve(LU,b,seq)
      print("   For R={:2d} ohms, Loop currents (amps): i1={:5.4f}, i2={:5.4f}, i3={:5.4f}".format(
             R,i[0],i[1],i[2]))


def problem20():
   "problem set 2.2, problem 20, page 83"
   import LUdecomp3

   print("Problem 20: solving with LUdecomp3")
   # matrix is tridiagonal, so using LUdecomp3 to
   # do something different.

   td_c = np.array([8,6,3,2],dtype=float)
   td_d = np.array([-8,-10,-11,-7,-4],dtype=float)
   td_e = np.array([4,2,5,4],dtype=float)
   b = np.array([-80,0,0,0,-30],dtype=float)
   td_c,td_d,td_e = LUdecomp3.LUdecomp3(td_c,td_d,td_e)
   c = LUdecomp3.LUsolve3(td_c,td_d,td_e,b)
   print("   Concentrations: c1={:5.3f}, c2={:5.3f}, c3={:5.3f}, c4={:5.3f}, c5={:5.3f}".format(
        c[0],c[1],c[2],c[3],c[4]))



if __name__ == "__main__" :
   print(__doc__)  # module doc from beginning of file
   problem16()
   print()
   problem17()
   print()
   problem20()





