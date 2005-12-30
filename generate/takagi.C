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

#define MAX_TERMS 1000
#define MANTISSA_BITS 48

static char bits [MAX_TERMS+MANTISSA_BITS];

static void setup(void)
{
	int i;
	for (i=0; i<MAX_TERMS+MANTISSA_BITS; i++)
	{
		bits[i] = (int) (2.0*rand()/(RAND_MAX+1.0));
printf ("duuude its %d\n", bits[i]);
	}
}

// return (2^n x)
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
}

static void takagi_series_c (double re_q, double im_q, double *prep, double *pimp)
{
	//*prep = rep;
	//*pimp = imp;
}

static double takagi_series (double re_q, double im_q, int itermax)
{
	if (!is_init) { is_init = 1; setup(); }

	double rep, imp;
	takagi_series_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	return rep;
}


DECL_MAKE_HISTO(takagi_series);

/* --------------------------- END OF LIFE ------------------------- */
