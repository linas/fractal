/*
 * euler.c
 *
 * Implement Euler summation for some complex-valued arithmetic function.
 * viz, given f(n) for 0 < n return
 *
 * sum_{n=0}^\infty 2^{-(n+1)} \sum_{k=0}^n {n choose k} f(k+1)
 *
 * See wikipedia for more about Euler summation.
 *
 * March 2019  -- Linas Vepstas
 */
#ifndef EULER_H__
#define EULER_H__

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

typedef complex arithmetic(int);

complex euler_sum(arithmetic);


#ifdef __cplusplus
};
#endif

#endif /* EULER_H__ */

complex euler_sum(arithmetic f)
{
	complex sum = 0.0;
	double tn = 2.0;
	int n = 0;
	while (1)
	{
		n++;
		tn *= 2.0;
	}
	return 0.0;
}
