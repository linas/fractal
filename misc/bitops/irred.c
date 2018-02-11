/*
 * irred.c
 *
 * Find and verify irreducible golden polynomials.
 * February 2018
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

bool divides(int biststr, int denom)
{
	
}

/* Return true only if bistr is irreducible */
bool is_irreducible(int biststr)
{
	for (int n=1; n<bitstr; n +=2)
	{
		if (divides(bitstr, n)) return false;
	}
	return true;
}

/* Return length of bitstr, length in bits */
int len(int bitstr)
{
	int len=0;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

int main(int argc, char* argv[])
{
	for (int n=1; n<20; n +=2)
	{
		printf("%d	%d	%d\n", n, len(n), is_irreducible(n));
	}
}
