
/*
 * cat-algebra.c
 *
 * Categorical algebra mobius to num-theo mobius
 *
 * Linas Vepstas November 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"

int zee(unsigned int i, unsigned int j)
{
	if (j%i == 0) return 1;
	return 0;
}

int en(unsigned int i, unsigned int j)
{
	if (i==j) return 0;
	return -zee(i,j);
}

int en_pow_k(unsigned int k, unsigned int i, unsigned int j)
{
	unsigned int r;
	int dot = 0;

	if (k==0) return (i==j)?1:0;
	if (k==1) return en(i,j);
	
	/* recursive matrix product */
	dot = 0;
	for(r=1; r<j; r++)
	{
		/* N(i,j) == 0 when j <= i*/
		/* But this is a very conservative estimate, since N^k
		 * flattens out even faster */
		dot += en_pow_k(k-1, i, r) * en(r, j);
	}

	return dot;
}

int mu(unsigned int i, unsigned int j)
{
	unsigned int k;

	int sum = 0;
	for(k=0; k<j; k++)
	{
		sum += en_pow_k(k,i,j);
	}

	return sum;
}


int main (int argc, char * argv[])
{
	int i,j;
	if (argc < 2)
	{
		fprintf (stdout, "Usage: %s <k>\n", argv[0]);
		exit(1);
	}

	i = atoi(argv[1]);
	for (j=i; j ; j += i)
	{
		int m = mu(i,j);
		int mm = moebius_mu(j/i);

		printf("duude %d \t%d %d\n", j, m, mm);
	}

	return 0;
}
