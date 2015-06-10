## module LUpivotComplete
''' LU with complete pivoting.

    Alternative implementation using outer product

    a,seqR,seqC = LUdecomp(a,tol=1.0e-9).
    LU decomposition of matrix [a] using scaled row pivoting.
    The returned matrix [a] = contains [U] in the upper
    triangle and the nondiagonal terms of [L] in the lower triangle.
    Note that [L][U] is a row-wise, col-wise permutation of the original [a];
    the permutations are recorded in the vectors seqR, seqC.
    
    x = LUsolve(a,b,seqR, seqC).
    Solves [L][U]{x} = {b}, where the matrix [a] = and the
    permutation vectors seqR, seqC are returned from LUdecomp.
'''

# modified version of LUpivotOP.py which is 
# modified version of LUpivot.py from textbook's source code
# An outer product is used to update the matrix
# Complete pivoting is used instead of partial pivoting
#
# Doug Heisterkamp
# last modifed: 2/11/14

import numpy as np
import swap
import error

def LUdecomp(a,tol=1.0e-9):
    n = len(a)
    seqR = np.array(range(n))
    seqC = np.array(range(n))
    
  # don't need scale factors in complete pivoting
    
    for k in range(0,n-1):
        
      # find the largest element in remaining matrix and
        t = np.argmax(np.abs(a[k:n,k:n]))
      # t is index into 1D version of a[k:n,k:n], so unravel back to 2D
        r,c = np.unravel_index(t, a[k:n,k:n].shape)  # location in submatrix
        r += k  # convert location to full matrix a
        c += k
      # swap a[r,c] to pivot location a[k,k]
      # Row interchange, if needed
        if abs(a[r,c]) <  tol: error.err('Matrix is singular')
        if r != k:
            swap.swapRows(a,k,r)
            swap.swapRows(seqR,k,r)
      # Col interchange, if needed
        if c != k:
            swap.swapCols(a,k,c)
            swap.swapRows(seqC,k,c)

      # Elimination using outer product       
        a[k+1:,k] /=a[k,k]
        a[k+1:,k+1:] -=  np.outer(a[k+1:,k],a[k,k+1:])
    return a,seqR,seqC

def LUsolve(a,b,seqR,seqC):
    n = len(a)
    
  # Rearrange constant vector; store it in [x]
    x = b[seqR]
        
  # Solution
    for k in range(1,n):  
        x[k] = x[k] - np.dot(a[k,0:k],x[0:k])
    x[n-1] = x[n-1]/a[n-1,n-1]    
    for k in range(n-2,-1,-1):
       x[k] = (x[k] - np.dot(a[k,k+1:n],x[k+1:n]))/a[k,k]

  # Rearrange x to undo column swaps 
    invseqC = np.zeros(n,dtype=seqC.dtype)
    for i in range(n):
       invseqC[seqC[i]] = i
    return x[invseqC]

