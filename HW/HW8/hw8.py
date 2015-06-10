"""
CS 3513, Homework 8
Example Solution
"""

# Doug Heisterkamp
# Last modified: 4/16/2014

import numpy 

import trapezoid

def example1(x):
   return 1.0/(1+x*x)

def example2(x):
   return x*numpy.sqrt(x*x + 9.0)

def example0(x):
   return numpy.sin(x)



def problem1(runs):
   "Evaluate with Monte Carlo method"
   print("Problem 1: solving integrals using Monte Carlo methods")

   def mcfsum(vf,lower,upper,n,b=100):
      "helper function to evaluate vectoized function vf n times over range lower..upper in groups of size b"
      fsum = 0.0
      n_div_b,n_remainder = divmod(n,b)
      for i in range(n_div_b):
         # evaluate groups of b samples at a time
         r = numpy.random.random_sample((b,))*(upper-lower)+lower
         fsum += vf(r).sum()
      if n_remainder > 0 :
         r = numpy.random.random_sample((n_remainder,))*(upper-lower)+lower
         fsum += vf(r).sum()
      return fsum
 
   # assumes random number generator has all ready been seeded
   b = 100
   samples = [100,10000,1000000]
   for f,lower,upper,mes,value in runs:
      vf = numpy.vectorize(f)
      print("   {}".format(mes))
      for n in samples:
         area = (upper-lower)*mcfsum(vf,lower,upper,n)/n
         print("       = {:16.14f} using {:7d} samples.  Relative error = {:5.3e}".format(area,n,abs((area-value)/value)))


def problem2(runs):
   "Numerical verification of rate of convergence for trapezoidal integration"

   print("Problem 2: Numerical verification of rate of convergence for trapezoidal integration")
   for f,lower,upper,mes,value in runs:
      print("   {}".format(mes))
      print("     {:>5} {:>8} {:>15} {:>12} {:>12}".format('n','h','T_h(f)', 'e_h', 'e_2h/e_h'))
      old_t = 0
      e2h = 0
      eh = 0
      for k in range(1,9):
         n = 2**(k-1)
         h = (upper-lower)/n
         t = trapezoid.trapezoid(f,lower,upper,old_t,k)
         old_t = t
         eh = abs(value-t)
         if k>1 :
            print("     {:5d} {:8.6f} {:15.13f} {:12.5e} {:12.5e}".format(n,h,t,eh,e2h/eh))
         else:
            print("     {:5d} {:8.6f} {:15.13f} {:12.5e}".format(n,h,t,eh))
         e2h = eh



def altRomberg(f,a,b,tol=1.0e-12):
   r = numpy.zeros((21,21))
   r[0,0] = trapezoid.trapezoid(f,a,b,0,1)
   for k in range(1,21):
      r[k,0] = trapezoid.trapezoid(f,a,b,r[k-1,0],k+1)
      for j in range(1,k+1):
         r[k,j] = (4**j * r[k,j-1] - r[k-1,j-1])/(4**j -1)
      if abs(r[k,k-1]-r[k-1,k-1]) < tol*max(abs(r[k,k-1]),1.0):
         return r[k,k],True,r[:k+1,:k+1].copy()
   return r[-1,-1],False,r


def problem3(runs):
   print("Problem 3: Romberg")
   printopts = numpy.get_printoptions() # get current settings
   numpy.set_printoptions(linewidth=250)
   for f,lower,upper,mes,value in runs:
      print("   {}".format(mes))
      v,flag,r = altRomberg(f,lower,upper)
      numpy.set_printoptions(precision=13)
      print('   Value = {}, Converge flag = {} , Table R = \n{}'.format(v,flag,r))
      er = numpy.tril(numpy.abs(r - v ))
      ratios = numpy.zeros((er.shape[0]-1,er.shape[0]-1))
      for col in range(ratios.shape[1]):
         for row in range(col,ratios.shape[0]):
            if er[row+1,col] > 0 :
               ratios[row,col] = er[row,col]/er[row+1,col]
            else :
               ratios[row,col] = numpy.inf
      numpy.set_printoptions(precision=4)
      print('   Error ratios = \n{}'.format(ratios))
   numpy.set_printoptions(printopts) # restore settings


if __name__ == "__main__" :
   print(__doc__)  # module doc from beginning of file

   print("Using scipy.integrate.romberg to calculate integral values")
   print("Note: Integral 0 is not part of the assignment.")
   import scipy
   import scipy.integrate
   labels = ['Integral 0 : int_0^pi sin(x) dx :',
             'Integral 1 : int_0^1  1/(1+x^2)  dx : ',
             'Integral 2 : int_0^4  x * sqrt(x^2 + 9) dx : ' ]
   area0 = scipy.integrate.romberg(example0, 0, numpy.pi,  tol=1e-12, rtol=1.e-13, show=False, divmax=50, vec_func=False)
   print("{} = {:15.12f} using scipy".format(labels[0],area0))

   area1 = scipy.integrate.romberg(example1, 0, 1,  tol=1e-12, rtol=1.e-13, show=False, divmax=50, vec_func=False)

   print("{} = {:15.12f} using scipy".format(labels[1],area1))
   area2 = scipy.integrate.romberg(example2, 0, 4,  tol=1e-12, rtol=1.e-13, show=False, divmax=21, vec_func=False)
   print("{} = {:15.12f} using scipy".format(labels[2],area2))

   runs = [ (example0,0,numpy.pi,labels[0] ,area0),(example1,0,1,labels[1],area1),(example2,0,4,labels[2],area2) ]
   problem1(runs)
   print()
   problem2(runs)
   print()
   problem3(runs)
   print()
