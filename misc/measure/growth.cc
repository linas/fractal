
/*
 * Growth to measure at 1/3
 * This is the spot where the minkowsi measure grows the fastest.
 *
 * Linas Vepstas September 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include "stern.h"

int main(int argc, char *argv[])
{

	unsigned long p, q, pm, qm, pmid, qmid;

	for (int level=0; level < 20; level += 2)
	{
		// numerator = 1 + 4 + 16 + 64 + 256 + ...
		unsigned long numerator = 0;
		unsigned long four = 1;
		for (int j=0; j< level/2; j++)
		{
			numerator += four;
			four <<= 2;
		}

		unsigned long denominator = 1<<level;
		double norm = 1.0 / (double)(denominator);

		// Converges to 1/3 from below
		double third = ((double) numerator) * norm;

		stern_brocot_tree(numerator, level, p, q);
		stern_brocot_tree(numerator+1, level, pm, qm);
		double a = ((double) p) / (double) q;
		double b = ((double) pm) / (double) qm;
		stern_brocot_tree(2*numerator+1, level+1, pmid, qmid);
		double y = ((double) pmid) / (double) qmid;

		double delta = norm * q * qm;

		printf("%d	%g	%lu	%lu	%g	%g\n", level, third, p, q, delta, y);
	}
}
