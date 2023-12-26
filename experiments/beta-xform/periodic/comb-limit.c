/*
 * comb-limit.c
 * Explore convergence of comb integral
 *
 * December 2023
 */

#include <math.h>
#include <stdio.h>

#include "necklace.h"

int main()
{
	for (int i=2; i<62; i++)
	{
		long mn = necklace(i);
		long twon = 1UL << i;
		double symp = ((double) mn) / ((double) twon);
		double asymp = 1.0/symp-i;
		double limit = pow(((double) mn), -symp);
		double tn = twon;
		double ll = pow(tn, 1.0/tn) - 1.0;

		printf("%d Mn=%ld	sym=%g  lim=%g ll=%g\n", i, mn, asymp, limit, ll);
	}

}
