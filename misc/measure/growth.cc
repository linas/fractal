
/*
 * Growth to measure at 1/3
 * This is the spot where the minkowsi measure grows the fastest.
 *
 * Linas Vepstas September 2015
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "stern.h"

int main(int argc, char *argv[])
{

	unsigned long p, q, pm, qm;

	for (int level=0; level < 64; level += 2)
	{
		// numerator = 1 + 4 + 16 + 64 + 256 + ...
		unsigned long numerator = 0;
		unsigned long four = 1;
		for (int j=0; j< level/2; j++)
		{
			numerator += four;
			four <<= 2;
		}

		unsigned long denominator = 1UL<<level;
		double norm = 1.0 / (double)(denominator);

		// Converges to 1/3 from below
		double third = 1.0 / 3.0 - ((double) numerator) * norm;

		stern_brocot_tree(numerator, level, p, q);
		stern_brocot_tree(numerator+1, level, pm, qm);
#if 0
		// a and b both converge to (2 - golden mean)
		double a = ((double) p) / (double) q;
		double b = ((double) pm) / (double) qm;
#endif

		double delta = norm * q * qm;
		double lim = log(delta) / (double) level;
		lim /= log(2.0);

		printf("%d	%9.6g	%lu	%lu	%18.16g	%18.16g\n", level, third, p, q, delta, lim);
		fflush(stdout);
	}
}
