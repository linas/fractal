/**
 * ffbasis.c
 *
 * Fallig factorial basis for the divisor function
 *
 * Linas Vepstas Decmber 2014
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "stirling.h"

long long int ortho(int p, int k)
{
	int n;
	long long int sum = 0LL;
	for (n=k; n<=p; n++)
	{
		sum += sitrling_first(p, n) * stirling_second(n, k);
	}
	return sum;
} 

int
main (int argc, char * argv[])
{
	int p, k;

	for (k=1; k<10 ; k++)
	{
		for (p=k; p<10; p++)
		{
			printf("duude k=%d p=%d sum=%lld\n", k, p, ortho(p, k));
		}
	} 

}

