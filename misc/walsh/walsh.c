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

/**
 * Generate the Rademacher functions r_n(x).
 * The counting used here is such that r_0(x) = 1 
 * r_1(x) = -1 for x<0.5 and 1 for x>0.5
 * r_2(x) = r_1(2*x)  or r_{n+1}(x) = r_n(2*x)
 *
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

main (int argc, char * argv[])
{
	
	printf("#\n# Rademacher/Waalsh functions\n#\n");
	int nsteps = 300;
	int i;
	for (i=0; i<= nsteps; i++)
	{
		double x = i / ((double) nsteps);
		double y = rademacher (x, 2);
		
		printf(%d	%f	%f\n", i, x, y);
	}
}
