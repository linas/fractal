/**
 * stirling.c
 *
 * Stirling numbers of the first and second kind
 *
 * Linas Vepstas <linasvepstas@gmail.com> Dec 2014
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "stirling.h"

long int stirling_first (int n, int k)
{
	return 0;
}

long int stirling_second (int n, int k)
{
	if (1 == k) return 1;
	if (1 == n) return 1;
	return 0;
}


#define TEST 1
#if TEST

main()
{
}

#endif
