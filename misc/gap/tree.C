/* 
 * tree.C
 *
 * Explore the eigenfunctions of bournoulli operator
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

// per driebe chapter 3
double 
enl_re (int n, int l, double x)
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
enl_im (int n, int l, double x)
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

double 
e_zl_re (double z_re, double z_im, int l, double x)
{
	int i;
	double znre = 1.0;
	double znim = 0.0;
	double acc = 0.0;
	i = 0;
	while (1)
	{
		acc += znre*enl_re(i,l,x) - znim * enl_im (i,l,x);
		double tmp = znre * z_re - znim * z_im;
		znim = z_re *znim + z_im * znre;
		znre = tmp;
		i++;
		double zabs = znre*znre+znim*znim;
		if (1.0e-16 > zabs) break;
	}
	return acc;
}

main (int argc, char *argv[])
{
	double z;
	int i;
	ContinuedFraction f;
	f.SetEvenize();

	if (2>argc)
	{
		printf ("Usage: %s <z-value>\n", argv[0]);
		exit (1);
	}

	z = 0.0;
	z = atof (argv[1]);
	int ell = 0;
	
	int nmax = 531;
	for (i=0; i<nmax; i++)
	{

		int p = i;
		int q = nmax;
		double x = ((double) p)/ ((double) q);
		
		double e = e_zl_re (z,0.0,ell, x);

		e *= sqrt(1.0-z);

		// f.SetReal (x);
		// double y = x-f.ToFarey();
		// double y = e_zl_re (x,0.0, ell, 0.33);
		// double e = y;
		
		printf("%5d	%8.6g	%8.6g\n", i,x,e);

	}
}

