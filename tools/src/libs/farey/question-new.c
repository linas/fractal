/*
 * question-new.c
 *
 * Straight-C routines for computing the Minkowski question mark function.
 *
 * Linas Vepstas May 2007
 */

#include <math.h>

#include "flt.h"
#include "question.h"

/**
 * dyadic_to_stern_brocot -- convert real number to stern-brocot
 *
 * This function maps the unit interval, expressed as a dyadic
 * tree, to the full stern brocot tree.  The input is 0<x<1.
 * Perform binary expansion of x. Use L,R for binary digits 0,1. 
 * Concatenate the L and R to form a matrix. Padd on the right with
 * zeros. The resulting fractinal linear transform will map the 
 * entire complex plane to b/d. This function returns b/d.
 *
 * Basically, this function maps the dyadic tree to the 
 * Stern-Brocot tree.
 */
long double dyadic_to_stern_brocot (long double x)
{
	if (0.0L>x) x -= floorl(x) - 1.0L;
	if (1.0L<x) x -= floorl(x);

	mobius_t ell = mobius_set (1.0L, 0.0L, 1.0L, 1.0L);
	mobius_t are = mobius_set (1.0L, 1.0L, 0.0L, 1.0L);

	mobius_t acc = mobius_set (1.0L, 0.0L, 0.0L, 1.0L);
	int i;

	/* plain double has a 45-bit mantissa */
	/* long double has 64 bits in mantissa (in gcc, right now) */
	for (i=0; i<64; i++)
	{
		if (0.5L <= x)
		{
			acc = mobius_mul (acc, are);
			x -= 0.5L;
		}
		else
		{
			acc = mobius_mul (acc, ell);
		}
		x *= 2.0L;
	}
	return acc.b.re / acc.d.re;
}

/**
 * question_inverse - return the inverse of the question mark function
 *
 * This implements a rapid algorithm to compute the inverse of 
 * the question mark function.
 */
long double question_inverse (long double x)
{
	return dyadic_to_stern_brocot (0.5L*x);
}

