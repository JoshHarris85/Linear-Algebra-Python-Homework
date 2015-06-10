print(
"""Homework 1,
CS 3513, Spring 2014

Question 2: Harmonic Number
""")

# CS 3513 Spring 2014
# Written by : Doug Heisterkamp
# last modified: 01/26/2014

from compat27 import  * # imports to make python version 2.7 behave like 3.2
import decimal     # for higher precision floating point

def CompensatedSum(x):
   'returns the compensated sum of a list, x, of numbers'
   n = len(x)
   c = 0
   s = x[0]
   for i in range(1,n):
      y = c + x[i]
      t = s + y 
      c = (s-t) + y
      s = t
   return s


decimal.getcontext().prec = 5
D = decimal.Decimal  # a shortcut name
# set precision of decimals to a low number for deomstation

psum = 0.0
for i in range(1,10001):
   psum = psum + 1.0/i

print("H_10000 = {:11.9f} using python floats and forward sum".format(psum))

dsum = D(0)
d1 = D(1)
for i in range(1,10001):
   dsum = dsum + d1/i 

bsum = D(0)
for i in range(10000,0,-1):
   bsum = bsum + d1/i 


print("H_10000 = {} using decimal floats with 5 digits precision and forward sum".format(dsum))
print("H_10000 = {} using decimal floats with 5 digits precision and reverse sum".format(bsum))

x = [ d1/i for i in range(1,10001) ] 

csum = CompensatedSum(x)
print("H_10000 = {} using decimal floats with 5 digits precision and compensated sum".format(csum))

# a snippet of code to wait for a key press before exiting in MS Windows
import platform
if "Win" in platform.uname()[0] :
   print("Press a key to exit")
   input()
