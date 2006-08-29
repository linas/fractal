/* 
 * affine.c
 *
 * Draw de Rham curves by iteration of affine matrices
 *
 * Linas Vepstas may 2005, august 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

typedef double affine[2][3];

affine d0;
affine d1;
affine result;

static inline copy (affine r, affine a)
{
	r[0][0] = a[0][0];
	r[0][1] = a[0][1];
	r[0][2] = a[0][2];
	
	r[1][0] = a[1][0];
	r[1][1] = a[1][1];
	r[1][2] = a[1][2];
}

static inline mult (affine r, affine a, affine b)
{
	r[0][0] = a[0][0]*b[0][0] + a[0][1] * b[1][0];
	r[0][1] = a[0][0]*b[0][1] + a[0][1] * b[1][1];
	r[0][2] = a[0][0]*b[0][2] + a[0][1] * b[1][2] + a[0][2];
	
	r[1][0] = a[1][0]*b[0][0] + a[1][1] * b[1][0];
	r[1][1] = a[1][0]*b[0][1] + a[1][1] * b[1][1];
	r[1][2] = a[1][0]*b[0][2] + a[1][1] * b[1][2] + a[1][2];
}

void fixpt (double val)
{
	affine tmp;
	
	int i = 0;
	val *= (double) (1<<30);
	unsigned int nt = (int) val;
	if (nt & 0x1) 
	{
		copy(result,d1);
	}
	else
	{
		copy(result,d0);
	}
	nt >>= 1;

	for (i=1; i<30; i++)
	{
		if (nt & 0x1) 
		{
			mult (tmp, d1, result);
			copy(result,tmp);
		}
		else
		{
			mult (tmp, d0, result);
			copy(result,tmp);
		}
		nt >>= 1;
	}
} 

static double affine_series (double re_q, double im_q, int itermax, double param)
{
	int i;
	int p,q;

	double ax, ay, d,e,f,g;
	q  = 23;
	ax = 0.5;
	ay = 1.0;
	d = 0.0;
	e = 0.6;
	f = 0.18;
	g = 0.6;

	e = re_q;
	f = im_q;
	g = param;

	d0[0][0] = ax;
	d0[1][0] = ay;
	d0[0][1] = d;
	d0[1][1] = e;
	d0[0][2] = 0.0;
	d0[1][2] = 0.0;
	
	d1[0][0] = 1.0-ax;
	d1[1][0] = -ay;
	d1[0][1] = f;
	d1[1][1] = g;
	d1[0][2] = ax;
	d1[1][2] = ay;
	
	q = itermax;
	double dist = 0.0;
	for (p=0; p<q; p++) 
	{
		double val = (double) p / (double) q;
		fixpt (val);

		double x = result[0][2];
		double y = result[1][2];
		dist += sqrt (x*x+y*y);
	}
	dist /= itermax;
	return dist;
}


DECL_MAKE_HEIGHT(affine_series);

/* --------------------------- END OF LIFE ------------------------- */
