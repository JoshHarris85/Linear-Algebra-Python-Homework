"""
CS 3513, Homework 5
Problem Set 3.1
Example Solution
"""

# Doug Heisterkamp
# Last modified: 2/23/2014

import numpy as np
import numpy.linalg as la


def makeLagrangePoly(i,dataX):
   "Returns the ith Lagrange polynomial of array dataX as a function"
   def Li(x,d=dataX):
      # brute force implementation 
      r = 1.0;
      for j in range(len(d)):
         if i == j :
            continue
         r *= (x-d[j])/(d[i]-d[j])
      return r
   return Li

def problem2():
   "problem set 3.1, problem 2, page 126"
   print("Problem 2")
   y = np.array([1.8421,2.4694,2.4921,1.9047,0.8509,-0.4112,-1.5727])
   x = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0])
   # part a, use three points, inverse interpolation, so using y
   L0 = makeLagrangePoly(0,y[4:])
   L1 = makeLagrangePoly(1,y[4:])
   L2 = makeLagrangePoly(2,y[4:])
   xroot = x[4]*L0(0) + x[5]*L1(0) + x[6] *L2(0) 
   print("    Using three points for Lagrange : y=0 at x = {:6.4f}".format(xroot))
   # part b, use four points, inverse interpolation, so using y
   L0 = makeLagrangePoly(0,y[3:])
   L1 = makeLagrangePoly(1,y[3:])
   L2 = makeLagrangePoly(2,y[3:])
   L3 = makeLagrangePoly(3,y[3:])
   xroot = x[3]*L0(0) + x[4]*L1(0) + x[5] *L2(0) + x[6]*L3(0) 
   print("    Using four points for Lagrange  : y=0 at x = {:6.4f}".format(xroot))


def problem3():
   "problem set 3.1, problem 3, page 126"
   import neville
   print("Problem 3")
   y = np.array([1.8421,2.4694,2.4921,1.9047,0.8509,-0.4112,-1.5727])
   x = np.array([0.0,0.5,1.0,1.5,2.0,2.5,3.0])
   p = 0.7692
   print("   Maximum y at x={:6.4f} is {:6.4f} using Neville's algorithm and four points".format(
      p,neville.neville(x[0:4],y[0:4],p)))
   

def problem9():
   "problem set 3.1, problem 9, page 127"
   def calcLC(d):
      "helper function to calc coeff of quadratic Langrange polynomials"
      assert(len(d)==3)  # only for quadratic
      a = np.zeros((3,3),float)
      a[0,0] = 1.0/((d[0]-d[1])*(d[0]-d[2]))  # L0 denominator, als0 h^2 coeff
      a[1,0] = 1.0/((d[1]-d[0])*(d[1]-d[2]))  # L1 denominator
      a[2,0] = 1.0/((d[2]-d[0])*(d[2]-d[1]))  # L2 denominator
      a[0,1] = (-d[1] -d[2]) * a[0,0]    # L0, h coeff
      a[1,1] = (-d[0] -d[2]) * a[1,0]    # L1, h coeff
      a[2,1] = (-d[0] -d[1]) * a[2,0]    # L2, h coeff
      a[0,2] = d[1]*d[2] * a[0,0]    # L0, h^0 coeff
      a[1,2] = d[0]*d[2] * a[1,0]    # L1, h^0 coeff
      a[2,2] = d[0]*d[1] * a[2,0]    # L2, h^0 coeff
      return a

   print("Problem 9")
   h = np.array([0,3,6],float)
   p = np.array([1.225,0.905,0.652])
   a = calcLC(h)
   c = p[0]*a[0,:] + p[1]*a[1,:] + p[2]*a[2,:]
   print("   Using Lagrange interpolation.")
   print("      L0(h) = {:6.4f}*h^2 + {:6.4f}*h + {:6.4}".format(a[0,0],a[0,1],a[0,2]))
   print("      L1(h) = {:6.4f}*h^2 + {:6.4f}*h + {:6.4}".format(a[1,0],a[1,1],a[1,2]))
   print("      L2(h) = {:6.4f}*h^2 + {:6.4f}*h + {:6.4}".format(a[2,0],a[2,1],a[2,2]))
   print("      p(h) = {:6.3f}*L0(h) + {:6.3f}*L1(h) + {:6.3f}*L2(h)".format(p[0],p[1],p[2]))
   print("   p(h) = {:8.6f}*h^2 + {:7.4f}*h + {:6.4}\n".format(c[0],c[1],c[2]))


def problem12():
   "problem set 3.1, problem 12, page 127"
   print("Problem 12")
   import cubicSpline
   # for cubic spline, need data in ascending order
   a = np.array( [ [0.2, 1.150],  # [x,y] pairs
                   [0.4, 0.855],
                   [0.6, 0.377],
                   [0.8,-0.266],
                   [1.0,-1.049]])
   b   = a[ :, [1,0]].tolist()  # swap x,y location
   b.sort()  # sort on y
   c = np.array(b)
   k = cubicSpline.curvatures(c[:,0],c[:,1])
   xroot = cubicSpline.evalSpline(c[:,0],c[:,1],k,0.0)
   print("    Using cubic spline for inverse interpolation.")
   print("    y=0 at x = {:6.4f}".format(xroot))

def problem15():
   "problem set 3.1, problem 15, page 128"
   print("Problem 15")
   import newtonPoly
   import rational
   import numpy as np
   data = np.array( [ [-250,0.0163], 
                      [-200,0.318],
                      [-100,0.699],
                      [   0,0.870],
                      [ 100,0.941],
                      [ 300,1.04]])
   # using Newton's method for polynomial interpolation
   t = np.linspace(-250,500,200)
   cc = newtonPoly.coeffts(data[:,0],data[:,1])
   np = [ newtonPoly.evalPoly(cc,data[:,0],i) for i in t]
   rp = [ rational.rational(data[:,0],data[:,1],i) for i in t]
   import matplotlib.pyplot as plt
   plt.plot(t,np,"b-",label="Polynomial interpolation")
   plt.plot(t,rp,"r--",label="Rational polynomial interpolation")
   plt.plot(data[:,0],data[:,1],"go",label="data points")
   plt.xlabel("Temperature ($^{\circ}$C)")
   plt.ylabel(r"Specific heat $c_{p}$   $\left(\frac{kJ}{kg\cdot K}\right)$")
   plt.title("Specific heat of Aluminum")
   plt.legend(loc='upper left', fancybox=True, shadow=True)
   plt.savefig("hw5p15.pdf")
   print("   Figure saved as hw5p15.pdf")



def problem17():
   "problem set 3.1, problem 17, page 128"
   print("Problem 17")
   data = np.array( [ [0.2  , 103 ], #[Re,c_d] pairs
                      [2    , 13.9],
                      [20   , 2.72],
                      [200  , 0.80],
                      [2000 , 0.401],
                      [20000, 0.433]],float)
   # switch data to log-log scale
   logdata = np.log(data)
   import cubicSpline
   k = cubicSpline.curvatures(logdata[:,0],logdata[:,1])
   print("    Using cubic spline over log-log data.")
   for r in [5,50,500,5000]:
      # note: map r to log(r), eval spline, map results back using exp
      print("    At Re={:5}, c_D = {:6.4f}".format(r,
             np.exp(cubicSpline.evalSpline(logdata[:,0],logdata[:,1],k,np.log(r)))))
   # for more info, plot the interpolation
   import matplotlib.pyplot as plt
   plt.clf()
   def logSpline(r,lx=logdata[:,0],ly=logdata[:,1],curv=k):
      return np.exp(cubicSpline.evalSpline(lx,ly,curv,np.log(r)))
   vlSpline = np.vectorize(logSpline)
   lt = np.linspace(logdata[0,0],logdata[-1,0],2000)
   t = np.exp(lt)
   lcs = vlSpline(t)
   plt.loglog(t,lcs,"b-",label="cubic spline,log-log")
   k2 = cubicSpline.curvatures(data[:,0],data[:,1])
   def spline(r,x=data[:,0],y=data[:,1],curv=k2):
      return cubicSpline.evalSpline(x,y,curv,r)
   vSpline = np.vectorize(spline)
   cs = vSpline(t)
   plt.loglog(t,cs,"r--",label="cubic spline")
   plt.loglog(data[:,0],data[:,1],"go",label="data points")
   r = [5,50,500,5000]
   plt.loglog(r,vlSpline(r),"r^",label="Estimated points")
   plt.xlabel("Reynold's number")
   plt.ylabel(r"drag coefficient $c_{D}$ of sphere")
   plt.title("Drag coefficient of a sphere")
   plt.legend(loc='upper right', fancybox=True, shadow=True)
   plt.savefig("hw5p17a.pdf")

   # redo plot in original data space, not log-log
   plt.clf()
   plt.plot(t,lcs,"b-",label="cubic spline,log-log")
   plt.plot(t,cs,"r--",label="cubic spline")
   plt.plot(data[:,0],data[:,1],"go",label="data points")
   plt.plot(r,vlSpline(r),"r^",label="Estimated points")
   plt.xlabel("Reynold's number")
   plt.ylabel(r"drag coefficient $c_{D}$ of sphere")
   plt.title("Drag coefficient of a sphere")
   plt.legend(loc='upper right', fancybox=True, shadow=True)
   plt.savefig("hw5p17b.pdf")



if __name__ == "__main__" :
   print(__doc__)  # module doc from beginning of file
   problem2()
   print()
   problem3()
   print()
   problem9()
   print()
   problem12()
   print()
   problem15()
   print()
   problem17()
   print()

