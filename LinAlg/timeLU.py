import numpy as np
import numpy.linalg as la
import LUdecomp
import LUpivot
import LUpivotOP
import LUpivotComplete

def timeLUdecomp(a):
   for m in a:
      LU = LUdecomp.LUdecomp(m)

def timeLUpivot(a):
   for m in a:
      LU,seq = LUpivot.LUdecomp(m)

def timeLUpivotOP(a):
   for m in a:
      LU,seq = LUpivotOP.LUdecomp(m)

def timeLUpivotComplete(a):
   for m in a:
      LU,seqR,seqC = LUpivotComplete.LUdecomp(m)


if __name__ == "__main__" :
   import timeit
   np.random.seed(3513801)
   k = 10
   n = 200 
   s = 10
   print("Creating matrices : {} random matrices of {}x{}".format(s,n,n))
   a = np.random.random((s,n,n))*100
   print("starting timings")
   t1 = min(timeit.repeat('timeLUdecomp(a)',setup = 'from __main__ import timeLUdecomp,a',number=1,repeat=k))
   print("timeLUdecomp best time {} of {} runs".format(t1,k))
   t2 = min(timeit.repeat('timeLUpivot(a)',setup = 'from __main__ import timeLUpivot,a',number=1,repeat=10))
   print("timeLUpivot best time {} of {} runs".format(t2,k))
   t3 = min(timeit.repeat('timeLUpivotOP(a)',setup = 'from __main__ import timeLUpivotOP,a',number=1,repeat=10))
   print("timeLUpivotOP best time {} of {} runs".format(t3,k))
   t4 = min(timeit.repeat('timeLUpivotComplete(a)',setup = 'from __main__ import timeLUpivotComplete,a',number=1,repeat=10))
   print("timeLUpivotComplete best time {} of {} runs".format(t4,k))


