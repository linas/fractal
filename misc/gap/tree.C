
/* tree.C
 *
 * Explore other tree structures
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"

// per driebe chapter 3
double 
e_re (int n, int l, double x)
{
	int i;
	x *=  2*l+1;
	for (i=0; i<n; i++)
	{
		x *= 2.0;
	}
	x *= 2.0 * M_PI;
	double enl = cos (x); 
	return enl;
}

double 
e_im (int n, int l, double x)
{
	int i;
	x *=  2*l+1;
	for (i=0; i<n; i++)
	{
		x *= 2.0;
	}
	x *= 2.0 * M_PI;
	double enl = sin (x); 
	return enl;
}

main (int argc, char *argv[])
{
	double z;
	int i;
	ContinuedFraction f;
	f.SetEvenize();

	int m = 1;
	int n = 1;
	
	int nmax = 523;
	for (i=0; i<nmax; i++)
	{

		int p = i;
		int q = nmax;
		double x = ((double) p)/ ((double) q);
		
#ifdef QUADRATIC_PARABOLA
		// quadratic parabiola y=4(x-1/2)^2
		p += nmax;
		q *= 2;
		int pp =(2*p-q)*(2*p-q);
		int qq = q*q;
		f.SetRatio (pp,qq);
#endif

		double z = x-0.5;
		z *= 2.0;
		z = 2*z*z*z-z;
		z += 1.0;
		z *= 0.5;

		// z = (2x-1)^3 -x +1
		// z(0) = 0  and z(1) = 1

		int pp = 2*p-q;
		pp = pp*pp*pp + (q-p)*q*q;
		int qq = q*q*q;
		f.SetRatio (pp,qq);
		
		double y = f.ToFarey();
		
		printf("%5d	%8.6g	%8.6g	%8.6g\n", i,x,y, z);

	}
}

