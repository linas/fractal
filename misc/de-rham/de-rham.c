
/* de-rham.C
 *
 * draw de Rham curves by iteration of functions
 *
 * Linas Vepstas may 2005
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"

static double ax = 0.5;
static double ay = sqrt(3.0)/6.0;
static double w = 0.6;

void koch_0 (double *x, double *y)
{
	double re = *x;
	double im = *y;
	*x = ax * re + ay * im;
	*y = ay * re - ax * im;
}

void koch_1 (double *x, double *y)
{
	double re = *x;
	double im = *y;
	*x = (1.0-ax) * re - ay * im;
	*y = -ay * re - (1.0-ax) * im;
	*x += ax;
	*y += ay;
}

void cesaro_0 (double *x, double *y)
{
	double re = *x;
	double im = *y;
	*x = ax * re - ay * im;
	*y = ay * re + ax * im;
}

void cesaro_1 (double *x, double *y)
{
	double re = *x;
	double im = *y;
	*x = (1.0-ax) * re + ay * im;
	*y = -ay * re + (1.0-ax) * im;
	*x += ax;
	*y += ay;
}

void bernoulli_0 (double *x, double *y)
{
	double re = *x;
	double im = *y;
	*x = ax * re - ay * im;
	*y = ay * re + ax * im;
}

void bernoulli_1 (double *x, double *y)
{
	double re = *x;
	double im = *y;
	*x = (1.0-ax) * re + ay * im;
	*y = -ay * re + (1.0-ax) * im;
}

void mink_0 (double *x, double *y)
{
	*x /= *x+1.0;
}

void mink_1 (double *x, double *y)
{
	*x += 1.0;
}

void takagi_0 (double *x, double *y)
{
	double xx = *x;
	double yy = *y;
	*x = xx*0.5;
	*y = xx + w*yy;
}

void takagi_1 (double *x, double *y)
{
	double xx = *x;
	double yy = *y;
	*x = 0.5*xx + 0.5;
	*y = -xx +w *yy +1.0;
}

void fixpt (double val, double *x, double *y)
{
	int i = 0;
	val *= (double) (1<<30);
	unsigned int nt = (int) val;
	for (i=0; i<30; i++)
	{
		if (nt & 0x1) 
		{
			// koch_1 (x,y);
			// cesaro_1 (x,y);
			// mink_1 (x,y);
			// bernoulli_1 (x,y);
			takagi_1 (x,y);
		}
		else
		{
			// koch_0 (x,y);
			// cesaro_0 (x,y);
			// mink_0 (x,y);
			// bernoulli_0 (x,y);
			takagi_0 (x,y);
		}
		nt >>= 1;
	}
} 

main (int argc, char *argv[])
{
	int i;
	int p,q;

	if (4 > argc)
	{
		fprintf (stderr, "Usage: %s <q> <ax> <ay>\n", argv[0]);
		exit (1);
	}
	
	q  = 43;
	q = atoi (argv[1]);
	ax = atof (argv[2]);
	ay = atof (argv[3]);
	w = ay;

	double aa = sqrt (ax*ax+ay*ay);
	double ama = sqrt ((1.0-ax)*(1.0-ax)+ay*ay);
	double c = sqrt (aa*aa+ama*ama);
	printf ("#\n# denom=%d\n#\n", q);
	printf ("# x=%g y=%g |a| = %g |1-a|=%g c=%g\n#\n", ax,ay, aa, ama,c);

	double x = 0.5;
	double y = 0.0;
	for (i=0; i<30; i++) cesaro_0 (&x,&y);
	printf ("# fxpt0 = %g + i %g\n", x ,y );
	cesaro_1 (&x,&y);
	printf ("#\n# 1-fxpt0 = %g + i %g\n", x ,y );

	x = 1.0;
	y = 0.0;
	printf ("#\n# fxpt1 = %g + i %g\n", x ,y );
	cesaro_0 (&x,&y);
	printf ("#\n# 0-fxpt1 = %g + i %g\n", x ,y );
	printf ("#\n");

	for (p=0; p<q; p++) 
	{
		double val = (double) p / (double) q;
		x = 0.5;
		y = 0.0;
		fixpt (val, &x, &y);

		printf ("%d\t%g	%g	%g\n", p,val,x,y);
	}
}
