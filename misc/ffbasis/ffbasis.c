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
		sum += stirling_first(p, n) * stirling_second(n, k);
	}
	return sum;
} 

int
main (int argc, char * argv[])
{
	int p, k;

#define MAX 30
	for (k=1; k<MAX ; k++)
	{
		for (p=k; p<MAX; p++)
		{
			long long int sum = ortho(p, k);
			if (p != k && 0 != sum) printf("duude k=%d p=%d sum=%lld\n", k, p, sum);
			if (p == k && 1 != sum) printf("duude k=%d p=%d sum=%lld\n", k, p, sum);
		}
		printf("--- k=%d\n", k);
	} 
	printf("a-K\n");
	return 0;
}

