/*
 * Compute the divisor operator and the mobius operator
 *
 * Linas Vepstas November 2014
 */

#include <stdio.h>

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
int moebius_oper(int n, int k)
{
	int m;
	if (n == k) return 1;
	if (n > k) return 0;
	int sum = 0;
	for (m=n; m<k; m++)
	{
		sum += moebius_oper(n, m) * divides(m, k);
	}

	return -sum;
}


main()
{
	int k;
	for (k=1; k<20; k++)
	{
		printf("mu(%d) = %d\n", k, moebius_oper(1,k));
	}
}
