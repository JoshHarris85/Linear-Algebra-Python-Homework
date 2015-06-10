## module LUpivotOP
''' Alternative implementation using outer product

    a,seq = LUdecomp(a,tol=1.0e-9).
    LU decomposition of matrix [a] using scaled row pivoting.
    The returned matrix [a] = contains [U] in the upper
    triangle and the nondiagonal terms of [L] in the lower triangle.
    Note that [L][U] is a row-wise permutation of the original [a];
    the permutations are recorded in the vector {seq}.
    
    x = LUsolve(a,b,seq).
    Solves [L][U]{x} = {b}, where the matrix [a] = and the
    permutation vector {seq} are returned from LUdecomp.
'''

# modified version of LUpivot.py from textbook's source code
# An outer product is used to update the matrix
# Doug Heisterkamp
# last modifed: 2/1/14

import numpy as np
import swap
import error

def LUdecomp(a,tol=1.0e-9):
    n = len(a)
    seq = np.array(range(n))
    
  # Set up scale factors
    s = np.zeros((n))
    for i in range(n):
        s[i] = max(abs(a[i,:]))        
    
    for k in range(0,n-1):
        
      # Row interchange, if needed
        p = np.argmax(np.abs(a[k:n,k])/s[k:n]) + k
        if abs(a[p,k]) <  tol: error.err('Matrix is singular')
        if p != k:
            swap.swapRows(s,k,p)
            swap.swapRows(a,k,p)
            swap.swapRows(seq,k,p)
      # Elimination using outer product       
        a[k+1:,k] /=a[k,k]
        a[k+1:,k+1:] -=  np.outer(a[k+1:,k],a[k,k+1:])
    return a,seq

def LUsolve(a,b,seq):
    n = len(a)
    
  # Rearrange constant vector; store it in [x]
    x = b.copy()
    for i in range(n):
        x[i] = b[seq[i]]
        
  # Solution
    for k in range(1,n):
        x[k] = x[k] - np.dot(a[k,0:k],x[0:k])
    x[n-1] = x[n-1]/a[n-1,n-1]    
    for k in range(n-2,-1,-1):
       x[k] = (x[k] - np.dot(a[k,k+1:n],x[k+1:n]))/a[k,k]
    return x

