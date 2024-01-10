/*
 * visits.c
 *
 * Visitation function, from the diary.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "selfie.c"
#include "selfie-util.c"

/*
 * Arrange midpoint iterations into a tree, so that each midpoint
 * splits a prior midpoint interval into two.  This generates a kind
 * of prefered Borel set (prefered base of a topology?), where the
 * boundaries are always at midpoints.
 *
 * visit: array of `depth` entries, filled in with location of the
 *        i'th iteration in the midpoint tree. The zeroth iteration
 *        m_0=beta/2 is always 0 and the first m_1=beta(beta-1)/2
 *        is always 1.
 *
 * midp:  array of `depth` entries, filled in with midpoint values.
 *
 * depth: How far to iterate.
 * epsilon: if greater than zero, then stop iteration after obn**steps
 *          is smaller than epsilon. If negative, then iterate forever,
 *          until depth is reached or overflow occurs,
 */
int visitation_tree(double beta,
                    unsigned long *visit, // unsigned long visit[DEPTH];
                    double *midp,         // double midp[DEPTH];
                    int depth,
                    double epsilon)
{
	unsigned long rigb[depth];

	// Setup
	for (int i=0; i<depth; i++) visit[i] = NEG_ONE;
	visit[0] = 0;
	rigb[0] = 1;

	double midpnt = 0.5*beta;
	midp[0] = midpnt;
	double obn = 1.0;
	int tot = 0;
	for (int i=1; i<depth; i++)
	{
		// Step to next midpoint
		if (0.5 < midpnt) midpnt -= 0.5;
		midpnt *= beta;
		obn /= beta;

		if (0.0 < epsilon && obn < epsilon) break;

		midp[i] = midpnt;

		// Search
		double rmost = 4200;
		int rig = -1;
		for (int j=0; j<i; j++)
		{
			if (rmost < midp[j]) continue;
			if (midpnt < midp[j] && midp[j] < rmost)
			{
				rmost = midp[j];
				rig = j;
			}
		}

		if (MAXIDX < rigb[rig]) break;

		// Record right-hand boundary, and update it.
		visit[i] = rigb[rig];
		rigb[i] = 2*rigb[rig];
		rigb[rig] = 2*rigb[rig]+1;
		// printf("i=%d mp=%g least=%g rig=%d -- visit= %ld\n",
		//	i, midpnt, rmost, rig, visit[i]);
		tot = i;
	}
	tot++;

	return tot;
}

// Same as above, but beta is gold midpoint given by index.
int visit_tree(unsigned long idx,
               unsigned long *visits, // unsigned long visit[DEPTH];
               double *midpts)        // double midp[DEPTH];
{
	double gold = golden_beta(idx);
	int order = bitlen(idx) + 1;

	return visitation_tree(gold, visits, midpts, order, -1.0);
}
