/*
 * Compute the divisor operator and the mobius operator
 *
 * Linas Vepstas November 2014
 */

#include <stdio.h>
#include <stdlib.h>

/* return 1 if d divides n, else return zero */
/* This is the divisor operator */
static inline int divides(int d, int n)
{
	if (d <= 0) return 0;
	if (d == 1) return 1;
	if (n%d == 0) return 1;
	return 0;
}

/* Inverse of the divisor function
 * Very slow recursive algo
 */
int moebius_oper_slow(int n, int k)
{
	int m;
	if (n == k) return 1;
	if (n > k) return 0;
	int sum = 0;
	for (m=n; m<k; m++)
	{
		sum += moebius_oper_slow(n, m) * divides(m, k);
	}

	return -sum;
}


/* Faster version of above. */
int moebius_oper(int n, int k)
{
#define SZ 1900
	static char cache[SZ][SZ];
	static bool done[SZ][SZ];
	static bool is_init = false;
	if (not is_init)
	{
		is_init = true;
		for (int i=0; i<SZ; i++)
			for (int j=0; j<SZ; j++)
				done[i][j] = false;
	}

	if (n == k) return 1;
	if (n > k) return 0;

	if (done[n][k]) return cache[n][k];

	int sum = 0;
	for (int m=n; m<k; m++)
	{
		sum += moebius_oper(n, m) * divides(m, k);
	}

	cache[n][k] = -sum;
	done[n][k] = true;
	return -sum;
}


main(int argc, char *argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s k\n", argv[0]);
		exit(1);
	}
	int n = atoi(argv[1]);


	for (n=1; n<SZ; n++)
	{
		int k;
		for (k=1; n*k<SZ; k++)
		{
			for (int j=1; j<n; j++)
			{
				if (0 != moebius_oper(n, n*(k-1)+j)) printf("Error: bad zero n=%d k=%d j=%d\n", n, k, j);
			}
			if (moebius_oper(1,k) != moebius_oper(n, n*k)) printf("Error: bad scale n=%d k=%d\n", n, k);
		}
	}
#if 0
	for (int k=1; k<20; k++)
	{
		printf("mu(%d) = %d vs %d\n", k, moebius_oper(2,k), moebius_oper(1, k/2));
	}
#endif

	printf("done\n");
}
