/*
 * FILE: find_zero.h
 * FUNCTION: find zero of a function
 * HISTORY:
 * July 2003 -- port from FORTRAN to C
 */

double FindZero (double (*func)(double), double bound, int nsig, 
          double a, 
          double b, int maxit);
