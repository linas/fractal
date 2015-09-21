
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

	unsigned __int128 p, q, pm, qm;

	for (int level=0; level < 128; level += 2)
	{
		// numerator = 1 + 4 + 16 + 64 + 256 + ...
		unsigned __int128 numerator = 0;
		unsigned __int128 four = 1;
		for (int j=0; j< level/2; j++)
		{
			numerator += four;
			four <<= 2;
		}

		unsigned __int128 denominator = 1UL<<level;
		long double norm = 1.0l / (long double)(denominator);

		// Converges to 1/3 from below
		// double third = 1.0 / 3.0 - ((double) numerator) * norm;

		// Interesting: p on a given row is q from previous row
		// i.e. p at level == q at level-2
		stern_brocot_tree(numerator, level, p, q);
		stern_brocot_tree(numerator+1, level, pm, qm);
#if 0
		// a and b both converge to (2 - golden mean)
		double a = ((double) p) / (double) q;
		double b = ((double) pm) / (double) qm;
#endif

		long double delta = q;
		delta *= qm;
		delta *= norm;
		long double lim = logl(delta) / (long double) level;
		lim /= M_LN2l;

		// printf("%d	%lu	%lu	%18.16g	%18.16g\n", level, q, qm, delta, lim);
		printf("%d\t", level);
		pr128(q);
		printf("\t");
		pr128(qm);
		printf("%26.24Lg	%18.16Lg\n", delta, lim);
		fflush(stdout);
	}
}
