/*
 * brot-prod.C
 *
 * FUNCTION:
 * Display A kind of Mandelbrot product function.
 *
 * HISTORY:
 * Linas - December 2008
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

static double brot_prod (double re_q, double im_q, int itermax, double param)
{
	int k;

	double rec = re_q;
	double imc = im_q;
	double rez = 0.0;
	double imz = 0.0;
	double rep = 1.0;
	double imp = 0.0;
	for (k=0; k<itermax; k++)
	{
		double tmp = rez*rez - imz*imz;
		imz = 2.0 * rez * imz;
		rez = tmp;
		rez += rec;
		imz += imc;

		tmp = rep * rez - imp * imz;
		imp = rep * imz + imp * rez;
		rep = tmp;

		// if (1.0e35 < rep*rep + imp*imp) break;
		if (1.0e35 < rez*rez + imz*imz) break;
	}

	return (double) k;
}

DECL_MAKE_HEIGHT(brot_prod);

/* --------------------------- END OF LIFE ------------------------- */
