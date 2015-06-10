"""
CS 3513, Homework 7
Problem Set 4.1
Example Solution
"""

# Doug Heisterkamp
# Last modified: 3/10/2014

import numpy as np
import numpy.linalg as la
import rootsearch
import ridder
import newtonRaphson
import newtonRaphson2
import matplotlib.pyplot as plt

def example_f(x):
   return x*np.sin(x) + 3*np.cos(x) - x

def example_df(x):
   return x*np.cos(x) -2*np.sin(x) - 1


def problem10():
   "problem set 4.1, problem 10, page 166"
   print("Problem 10: root search with Ridder's Method")
   lower,upper,dx = -6.0,6.0,0.1
   roots = []
   while True:
      print("      searching range {:+5.4f} ... {:+5.4f}".format(lower,upper))
      x1,x2 = rootsearch.rootsearch(example_f,lower,upper,dx)
      if x1 != None :
         lower = x2
         r = ridder.ridder(example_f,x1,x2)
         if r != None:
            roots.append(r)
      else:
         break
   if roots :
      print("   roots = {}".format(roots))
   else:
      print("    found no roots")


def problem11():
   "problem set 4.1, problem 11, page 166"
   print("Problem 11: root search with Newton's Method")
   lower,upper,dx = -6.0,6.0,0.1
   roots = []
   while True:
      print("      searching range {:+5.4f} ... {:+5.4f}".format(lower,upper))
      x1,x2 = rootsearch.rootsearch(example_f,lower,upper,dx)
      if x1 != None :
         lower = x2
         r = newtonRaphson.newtonRaphson(example_f,example_df,x1,x2)
         if r != None:
            roots.append(r)
      else:
         break
   if roots :
      print("   roots = {}".format(roots))
   else:
      print("    found no roots")


def problem23():
   "problem set 4.1, problem 23, page 170"
   print("Problem 23: intersecting circles")
   def f(x):
      "circles (x-2)^2 + y^2 -4 = 0, and x^2 + y-3)^2 -4 =0"
      f = np.zeros((2,),dtype=x.dtype)
      f[0] = (x[0] -2.0)**2 + x[1]**2 - 4
      f[1] = x[0]**2 + (x[1] -3)**2 - 4 
      return f
   # two points of intersection, just start in area of intersection
   x1 = np.array([2.0,3.0])
   res1 = newtonRaphson2.newtonRaphson2(f,x1)
   print("   intersection point = {} ".format(res1))  
   x2 = np.array([0.0,0.0])
   res2 = newtonRaphson2.newtonRaphson2(f,x2)
   print("   intersection point = {} ".format(res2))  
   plt.clf()
   # use parametric form for plotting circle
   plt.plot([x1[0],x2[0]],[x1[1],x2[1]],"g^",label="initial points")
   plt.plot([res1[0],res2[0]],[res1[1],res2[1]],"ro",label="intersection points")
   r = np.linspace(0,2*np.pi,200)
   xc1 = [ 2*np.cos(a)+2 for a in r]
   yc1 = [ 2*np.sin(a)+0 for a in r]
   plt.plot(xc1,yc1,"b-",label="$(x-2)^2+y^2=4$")
   xc2 = [ 2*np.cos(a)+0 for a in r]
   yc2 = [ 2*np.sin(a)+3 for a in r]
   plt.plot(xc2,yc2,"m-",label="$x^2+(y-3)^2=4$")
   plt.xlim([-3,6])
   plt.ylim([-3,6])
   plt.axes().set_aspect('equal')
   plt.title("Intersection of circles")
   plt.legend(loc="upper right", fancybox=True, shadow=True, fontsize=12)
   plt.savefig("problem23.pdf")

   

 

def problem26():
   "problem set 4.1, problem 26, page 170"
   print("Problem 26: circle through three points")
   def f(t):
      f = np.zeros((3,),dtype=t.dtype)
      f[0] = (8.21-t[0])**2 + t[1]**2 - t[2]**2
      f[1] = (0.34-t[0])**2 + (6.62-t[1])**2 - t[2]**2
      f[2] = (5.96-t[0])**2 + (-1.12-t[1])**2 - t[2]**2
      return f
   t = np.array([4.0,3.0,5.0]) # inital guess, center is rough midpoint of first two points
   res = newtonRaphson2.newtonRaphson2(f,t)
   print("  Circle : (x - {:7.6f})^2 + (y - {:7.6f})^2 = {:7.6f}^2".format(res[0],res[1],res[2]))
   plt.clf()
   # use parametric form for plotting circle
   try:
      seq = []
      res = newtonRaphson2.newtonRaphson2(f,t,s=seq,) # using modified version to return intermediate results
   except TypeError:
      print("local newtonRaphson2 is not modified.  Not calcuating intermediate results")
      # continue on.  seq = [] so it will be ignored
   points = np.array([[8.21,0.0],[0.34,6.62],[5.96,-1.12]])
   plt.plot(points[:,0],points[:,1],"ro",label="Data points")
   r = np.linspace(0,2*np.pi,200)
   xc = [ res[2]*np.cos(a)+res[0] for a in r]
   yc = [ res[2]*np.sin(a)+res[1] for a in r]
   plt.plot(xc,yc,"b-",label="final circle")
   init_xc = [ t[2]*np.cos(a)+t[0] for a in r]
   init_yc = [ t[2]*np.sin(a)+t[1] for a in r]
   plt.plot(init_xc,init_yc,"r--",label="initial circle")
   for j in seq:
      init_xc = [ j[2]*np.cos(a)+j[0] for a in r]
      init_yc = [ j[2]*np.sin(a)+j[1] for a in r]
      plt.plot(init_xc,init_yc,"g--")
   plt.title("Circle through three points")
   plt.xlim([-3,12])
   plt.ylim([-3,10])
   plt.axes().set_aspect('equal')
   plt.legend(loc="upper right", fancybox=True, shadow=True,fontsize=12)
   plt.savefig("problem26.pdf")

   

 

if __name__ == "__main__" :
   print(__doc__)  # module doc from beginning of file
   problem10()
   print()
   problem11()
   print()
   problem23()
   print()
   problem26()
   print()

