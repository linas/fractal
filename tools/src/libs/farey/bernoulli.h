/* 
 * bernoulli.h
 *
 * Bernoulli number function prototypes 
 */
#ifndef BERNOULLI_H__
#define BERNOULLI_H__

// return the nth' bernoulli number
// this is a fairly fast algorithm ...
// must have n>=0
long double bernoulli (int n);

// return the value of the n'th bernoulli polynomial
// at the value ''x''
long double bernoulli_poly (int n, long double x);

#endif /* BERNOULLI_H__ */
