/*
 * brat-gap-hair.C
 *
 * FUNCTION:
 * Graph the bernoulli-map continuous-spectrum fractal eigenfunctions
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

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

double 
e_zl_im (double z_re, double z_im, int l, double x)
{
	int i;
	double znre = 1.0;
	double znim = 0.0;
	double acc = 0.0;
	i = 0;
	while (1)
	{
		acc += znim*enl_re(i,l,x) + znre * enl_im (i,l,x);
		double tmp = znre * z_re - znim * z_im;
		znim = z_re *znim + z_im * znre;
		znre = tmp;
		i++;
		double zabs = znre*znre+znim*znim;
		if (1.0e-16 > zabs) break;
	}
	return acc;
}


void 
MakeHisto (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm)
{
   int		i,j, k;

   int globlen = sizex*sizey;
   for (i=0; i<globlen; i++) {
      glob [i] = 0.0;
   }

	double zre, zim;

	zre = -0.5;
   zim = 0.0;
	int ell = 0;

	double rescale = 0.25;
	double imscale = 0.25;

	double recent = re_center;
	double imcent = im_center;
	
	for (k=0; k<itermax; k++)
	{
		int p = k;
		int q = itermax;
		double x = ((double) p)/ ((double) q);
		
		double er = e_zl_re (zre,zim,ell, x);
		double ei = e_zl_im (zre,zim,ell, x);

		// printf("%5d	%8.6g	%8.6g	%8.6g\n", i,x,er, ei);

		er += recent;
		ei += imcent;

		i = (int) (rescale * er * (double) sizex);
		if (0>i) continue;
		if (i>=sizex) continue;

		j = (int) (imscale * ei * (double) sizey);
		if (0>j) continue;
		if (j>=sizey) continue;

		glob [j*sizex +i] ++;
	}

   /* renormalize */
	double r = ((double) sizex) / ((double) itermax);
   for (i=0; i<sizex*sizey; i++) 
	{
		glob [i] *= r;
   }
}

/* --------------------------- END OF LIFE ------------------------- */
