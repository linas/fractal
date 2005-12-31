/*
 * takagi.C
 *
 * FUNCTION:
 * display takagi on the complex unit disk
 * for a "generic" point.
 *
 * HISTORY:
 * Linas December 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

static int is_init = 0;

#define MAX_TERMS 2300
#define MANTISSA_BITS 48

static char bits [MAX_TERMS+MANTISSA_BITS];

// return (2^n x) - floor (2^n x)
//
static double get_x (int n)
{
	int i;
	double x=0.0;
	double tn = 0.5;
	for (i=0; i<MANTISSA_BITS; i++)
	{
		x += bits[i+n] * tn;
		tn *= 0.5;
	}

	return x;
}

static void set_x (double x)
{
	int i;
	x -= floor (x);
	for (i=0; i<MANTISSA_BITS; i++)
	{
		if (x >= 0.5)
		{
			bits[i] = 1;
			x -= 0.5;
		}
		else
		{
			bits[i] = 0;
		}
		x *= 2.0;
	}
}

static void setup(int shift)
{
	int i;
	// srand (331);
	for (i=0; i<shift; i++)
	{
		rand();
	}
	for (i=0; i<MAX_TERMS+MANTISSA_BITS; i++)
	{
		bits[i] = (int) (2.0*rand()/(RAND_MAX+1.0));
	}
}

// basic triangle wave
static double triangle (double x)
{
	double y = 2.0*x;
	if (x <= 0.5) return y;
	return 2.0-y;
}

static void takagi_series_c (double re_w, double im_w, double *prep, double *pimp)
{
	int i;

	double wmag = sqrt(re_w * re_w + im_w *im_w);
	double wcos = re_w / wmag;
	double wsin = im_w / wmag;

	if (wmag > 1.00001) 
	{
		*prep = 0.0;
		*pimp = 0.0;
		return;
	}

	double rn = 1.0;
	double cosn = 1.0;
	double sinn = 0.0;

	double resum = 0.0;
	double imsum = 0.0;
	for (i=0; i<MAX_TERMS; i++)
	{
		double x = get_x (i);
		x = triangle(x);
		resum += rn * cosn * x;
		imsum += rn * sinn * x;

		if (rn < 1.0e-18) break;

		rn *= wmag;
		double tmp = sinn * wcos + cosn * wsin;
		cosn = cosn * wcos - sinn * wsin;
		sinn = tmp;
	}

	*prep = resum;
	*pimp = imsum;
}

// This function sets up a fixed ''x'' value, and graphs the interior of
// the disk (|w| < 1) as a funciton of w.
// 
static double takagi_disk (double re_q, double im_q, int itermax)
{
	if (!is_init) 
	{ 
		is_init = 1; 
		setup(itermax);
		double z = get_x (0);
		printf ("disk for x=%22.18g\n", z);
	}

	double rep, imp;
	takagi_series_c (re_q, im_q, &rep, &imp);
	//return (atan2 (imp, rep)+M_PI) / (2.0*M_PI);
	return sqrt (rep*rep+imp*imp);
	// return rep;
	// return imp;
}

// This function sets up a cylinder slice, with |w| held constant
// but the phase of w varying in the horizontal direction
// In the veritical direction runs the ''x'' of the t_w(x).
// 
static double takagi_cylinder (double x, double y, int itermax)
{
	double w = 0.94;

	static double yprev = -999.0;

	if (y != yprev)
	{
		setup (itermax);
		set_x (y);
		yprev = y;
		double z = get_x (0);
		printf ("disk for x=%22.18g (x=%g y=%g)\n", z, x, y);
	}
 
	double re_q = w * cos (2.0 * M_PI * x);
	double im_q = w * sin (2.0 * M_PI * x);

	double rep, imp;
	takagi_series_c (re_q, im_q, &rep, &imp);
	return (atan2 (imp, rep)+M_PI) / (2.0*M_PI);
	// return sqrt (rep*rep+imp*imp);
	// return rep;
}

// DECL_MAKE_HISTO(takagi_disk);
DECL_MAKE_HISTO(takagi_cylinder);

/* --------------------------- END OF LIFE ------------------------- */
