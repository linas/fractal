
/* 
 * build and solve comb-like eignevalue equation
 *
 * Linas Vepstas September 2004
 */

#include <math.h>

/*
 * Build a matrix corresponding to the differential equation
 * y^4r'' + y^3r' (2-4ay(y-b)) + y^3r[4a^2y(y-b)^2 -2ay -4a(y-b)] + rsin(2piy)
 * where
 * the matrix expresses this in terms of the taylor expansion, 
 * with r(y) = sum a_n y^n
 */

double
get_matrix_elt (int p, int n, double a, double b)
{
	double elt = 0.0;
	int k;

	// multiply by sin (2piy)
	k = p-n;
	if ((k>=0) && (0 == k%2))
	{
		int i;
		double term = 1.0;
		if (1 == k%2) term = -1.0;
		for (i=0; i<k; i++)
		{
			term *= 2.0*M_PI / ((double) (i+1));
		}
		elt += term;
	}
	
	if (p == n+2)
	{
		elt += n*(n-1);   // y^4 r''
		elt += 2*n;       // y^3 r' 2		
	}
	else if (p == n+3)
	{
		elt += 4.0*a*b* ((double) n);  // y^3r' 4ayb
		elt += 4.0*a*b;                // y^3r (-4ab)
	}
	else if (p == n+4)
	{
		elt -= 4.0*a* ((double) n);  // y^3r' (-4ay^2)
		elt += 4.0*a*a*b*b;          // y^3 r  4a^2y b^2
		elt -= 2.0*a;                // y^3 r  (-2ay)
		elt -= 4.0*a;                // y^3 r  (-4ay)
	}
	else if (p == n+5)
	{
		elt -= 8.0*a*a*b;       // y^3 r 4a^2y (-2by)
	}
	else if (p == n+6)
	{
		elt += 4.0*a*a;       // y^3 r 4a^2y y^2
	}

	return elt;
}

main ()
{
	int i, j;
	int dim = 13;

	double a = 0.16;
	double b = 6.0;

	for (i=0; i<dim; i++)
	{
		for (j=0; j<dim; j++)
		{
			double v = get_matrix_elt (i,j, a, b);
			printf ("%g\t", v);
		}
		printf ("\n");
	}
}
