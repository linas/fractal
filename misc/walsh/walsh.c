/*
 * walsh.c
 *
 * Walsh functions provide a basis for square-integrable functions.
 * Walsh functions can be enumerated using integers that have non-zero
 * Mobius mu values. Walsh functions are even or odd, depending on 
 * whether mu is positive or negative.
 *
 * The code here creates the Walsh functions, indexing them with integers.
 * This allows various series to be explored.
 *
 * The 'prime factors' of the walsh functions are the Rademacher functions.
 *
 * Linas Vepstas July 2009
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Generate the Rademacher functions r_n(x).
 * The counting used here is such that r_0(x) = 1 
 * r_1(x) = -1 for x<0.5 and 1 for x>0.5
 * r_2(x) = r_1(2*x)  or r_{n+1}(x) = r_n(2*x)
 */
double rademacher (double x, int n)
{
	if (0 == n) return 1.0;
	if (1 == n)
	{
		if (x < 0.5) return -1.0;
		return 1.0;
	}
	
	int shift = 1 << (n-1);
	x *= shift;
	x -= floor(x);
	return rademacher(x, 1);
}

/**
 * Generate the walsh basis functions w_n(x)
 * The index is proceeds by binary counting. 
 * w_0(x) = r_0(x) = 1
 * w_1(x) = r_1(x)
 * w_2(x) = r_2(x)
 * w_3(x) = r_1(x) r_2(x)
 * w_4(x) = r_3(x)
 * w_5(x) = r_3(x) r_1(x)
 */
double walsh (double x, int n)
{
	if (0 == n) return 1.0;
	double prod = 1.0;
	int shift = 1;
	
	while(n != 0)
	{
		if (n & 0x1) prod *= rademacher(x, shift);
		n >>= 1;
		shift += 1;
	}
	return prod;
}

/*
 * Assorted series summation functions
 * of the form sum_n a_n r_n(x) where 
 * a_n are geometric series, etc. 
 */

double geo_series(double x, double g)
{
	int n;

	double acc = 0.0;
	double gn = 1.0;
	for (n=0; n<33; n++)
	{
		acc += gn * walsh(x, n);
		gn *= g;
	}

	return acc;
}


main (int argc, char * argv[])
{
	if (2 > argc)
	{
		printf ("Useage: %s n\n", argv[0]);
		exit(1);
	}

	// int n = atoi(argv[1]);
	double g = atof(argv[1]);
	
	printf("#\n# Rademacher/Walsh functions\n#\n");
	int nsteps = 300;
	int i;
	for (i=0; i < nsteps; i++)
	{
		double x = i / ((double) nsteps);
		// double y = rademacher (x, 2);
		// double y = walsh (x, n);
		double y = geo_series (x, g);
		
		printf("%d	%f	%f\n", i, x, y);
	}
}
