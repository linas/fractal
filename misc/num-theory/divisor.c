
/*
 * divisor.c
 *
 * Compute the divisor arithmetic function
 *
 * Linas Vepstas December 2005
 */

#include "divisor.h"

int divisor (int n)
{
	int acc = 0;
	int d;

	for (d=1; d<=n; d++)
	{
		if (n%d) continue;
		acc ++;
	}

	return acc;
}
