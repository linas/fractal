
/*
 * padic.C
 *
 * p-adic number games
 *
 * Linas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int
padic (int prime, int num, int denom)
{
	int cnt = 0;
	if (0 == num) return 0;
	if (0 == denom) return 0;

	while (0 == num%prime)
	{
		cnt ++;
		num /= prime;
	}
	while (0 == denom%prime)
	{
		cnt --;
		denom /= prime;
	}
	return cnt;
}


main (int argc, char * argv[]) 
{
	int i;
	int prime=3;
	int nmax = 720;

	if (3>argc)
	{
		printf ("Usage: %s <prime>  <nmax>\n", argv[0]);
		exit (1);
	}
	prime = atoi (argv[1]);
	nmax = atoi(argv[2]);
	
	for (i=0; i<nmax; i++)
	{
		double x = ((double) i) / ((double) nmax);

		int pp = rand();
		int qq = rand();
		pp %= qq;
		int v = padic (prime, pp, qq);

		double val = v;

		if (0.0 < val)
		{
			val = 0.5 - 1.0/val + sqrt (0.25+1.0/(val*val));
		}
		else if (0.0 > val)
		{
			val = 0.5 - 1.0/val - sqrt (0.25+1.0/(val*val));
		}
		if (0 == v) val = 0.5;
		
		x = ((double) pp) / ((double) qq);

		printf ("%8.6g	%d	%8.6g\n", x, v, val);
	}
}
