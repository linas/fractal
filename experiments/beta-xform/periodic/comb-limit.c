/*
 * comb-limit.c
 * Explore convergence of comb integral
 *
 * December 2023
 */

#include <stdio.h>

#include "necklace.h"

int main()
{
	for (int i=2; i<62; i++)
	{
		long mn = necklace(i);
		printf("%d Mn=%ld\n", i, mn);
	}
}
