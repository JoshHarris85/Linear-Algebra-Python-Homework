Josh Harris
Homework 1 
CS 3513
Due 1/29/14
Note: All code was in working order when submitted. If you have any problems feel free to contact me. All output is listed below. 


Part 1: 
Output:
Calculating sqrt(e^2).
Error in final value is 0E-299
Using true value to calc error.
iter           |e_n|  |e_n+1|/|e_n|^1  |e_n+1|/|e_n|^2  |e_n+1|/|e_n|^3
   0        4.671e+0         1.274e-1         2.729e-2         5.842e-3
   1        5.953e-1         8.786e-3         1.476e-2         2.480e-2
   2        5.230e-3         9.228e-7         1.764e-4         3.374e-2
   3        4.826e-9        7.882e-19        1.633e-10         3.383e-2
   4       3.804e-27        4.896e-55        1.287e-28         3.383e-2
   5       1.862e-81       1.174e-163        6.301e-83         3.383e-2
   6       2.19e-244
e_n[7] is zero, skipping rest of list.
Press a key to exit


The estimate of the order of convergence is Alpha α = 4
The asymptotic error constant is Lambda λ = 0.03383



Part 2:
As we can see from the ouput I got below:

(A) 9.787606036044348
(B) 9.7506
(C) 9.7873
(D) 9.7875

Since the precision in A is the highest it is obviously going to be the closest to a correct answer. 
The precision in B is lowered dramatically for each iteration of 1/j. It is a worse answer since there are less digits on each j.
The Precision in C is a bit better than B since we went from 10000 to 1 
Precison in D is in between because you are taking the current 1/j with the last 1/j and adding them each iteration to create a new s variable.