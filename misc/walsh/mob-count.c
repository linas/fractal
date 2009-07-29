/*
 * mob-count.c
 *
 * Mobius counting function
 *
 * Linas Vepstas July 2009
 */

#include <stdio.h>

int primes[] = {1,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,51,53,61,67,71,73,79,83,89,91,97};

int mob_n(int n)
{
	int acc = 1;
	int shift = 1;
	while (n != 0)
	{
		if (n & 0x1) acc *= primes[shift];
		n >>= 1;
		shift ++;
	}
	return acc;
}

main (int argc, char * argv[])
{
	int i;
	for (i=1; i<50123; i++)
	{
		int m = mob_n (i);
		printf("%d	%d\n", i, m);
	}
}
