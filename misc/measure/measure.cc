
/*
 * Stern-Brocot tree, new algo
 * This is the same old question mark, just a different API.
 *
 * XXX This has been copied to stern-brocot.c  This version
 * is dead.
 *
 * Linas Vepstas September 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include "stern.h"


int main(int argc, char *argv[])
{
	int level = atoi(argv[1]);

	unsigned long p, q, pm, qm, pmid, qmid;
	double sum = 0.0;

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

		sum += (b-a)* delta;
		printf("%d	%g	%lu	%lu	%g	%g	%g\n", i, x, p, q, y, delta, sum);
	}
}
