
/*
 * A recursive expression for the Minkowski measure
 * This is the same old question mark, just a different API.
 * Specifically, since the thing is nicely, recursively defined, its more
 * tractable.
 *
 * Linas Vepstas September 2015
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "stern.h"


int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s level\n", argv[0]);
		exit(1);
	}

	int level = atoi(argv[1]);

	unsigned __int128 p, q, pm, qm, pmid, qmid;
	double sum = 0.0;

	// Golden ratio
	double phi = 0.5 * (1.0 + sqrt(5.0));

	// Growth rate
	double grow = pow(0.5 * phi *phi, level) + 0.2 * phi * phi * phi;

	double norm = 1.0 / (double)(1<<level);
	for (int i=0; i< (1<<level); i++)
	{
		stern_brocot_tree(i, level, p, q);
		stern_brocot_tree(i+1, level, pm, qm);
		double x = i * norm;
		double a = ((double) p) / (double) q;
		double b = ((double) pm) / (double) qm;
		stern_brocot_tree(2*i+1, level+1, pmid, qmid);
		double y = ((double) pmid) / (double) qmid;
		// unsigned long det = pm * q - p * qm;
		double delta = norm * q * qm;
		delta /= grow;

		unsigned long pp = p;
		unsigned long qq = q;

		sum += (b-a)* delta;
		printf("%d	%g	%lu	%lu	%g	%g	%g\n", i, x, pp, qq, y, delta, sum);
	}
}
