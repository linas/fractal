
/*
 * Box-counting dimension.
 * An integral over the upper bound.
 *
 * A self-bound for recursive expression for the Minkowski measure
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

int maxlev = 0;

double upper_bound(int lowlev, int uplev, int i)
{
	unsigned __int128 p, q, pm, qm;

	// Golden ratio
	double phi = 0.5 * (1.0 + sqrt(5.0));

	double bound = 0.0;

	for (int ill = lowlev; ill <= uplev; ill++)
	{
		int level = ill;

		// Growth rate
		double grow = pow(0.5 * phi *phi, level) * 0.2 * phi * phi * phi;
		double norm = 1.0 / (double)(1<<level);

		int j = i >> (uplev-ill);
		int k = j;

		stern_brocot_tree(k, level, p, q);
		stern_brocot_tree(k+1, level, pm, qm);
		double delta = norm * q * qm / grow;

		if (bound < delta)
		{
			bound = delta;
			if (maxlev < level) maxlev = level;
		}
	}
	return bound;
}

double content(int level)
{
	double sum = 0.0;

	unsigned __int128 p, q, pm, qm;
	int uplev = level+2;
	for (int i=0; i< (1<<uplev); i++)
	{
		double bound = upper_bound(level, uplev, i);

		stern_brocot_tree(i, uplev, p, q);
		stern_brocot_tree(i+1, uplev, pm, qm);

		double a = ((double) p) / (double) q;
		double b = ((double) pm) / (double) qm;
		sum += bound * (b-a);
	}
	return sum;
}

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s level\n", argv[0]);
		exit(1);
	}

	int maxlev = atoi(argv[1]);

	for (int lev=0; lev< maxlev; lev++)
	{
		double cont = content(lev);
		printf("%d\t%g\n", lev, cont);
		fflush(stdout);
	}
}
