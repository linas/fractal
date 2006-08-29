/* 
 * affine.c
 *
 * Draw de Rham curves by iteration of affine matrices
 *
 * Linas Vepstas may 2005, august 2006
 */

#include <math.h>
#include <stdio.h>

static double ax = 0.5;
static double ay = sqrt(3.0)/6.0;
static double w = 0.6;

static double d,e,f,g;

void generic_0 (double *x, double *y)
{
	double re = *x;
	double im = *y;
	*x = ax * re + d * im;
	*y = ay * re + e * im;
}

void generic_1 (double *x, double *y)
{
	double re = *x;
	double im = *y;
	*x = ax + (1.0-ax) * re + f * im;
	*y = ay - ay * re + g * im;
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
			generic_1 (x,y);
		}
		else
		{
			generic_0 (x,y);
		}
		nt >>= 1;
	}
} 

main (int argc, char *argv[])
{
	int i;
	int p,q;

	if (8 > argc)
	{
		fprintf (stderr, "Usage: %s <q> <ax> <ay>\n", argv[0]);
		exit (1);
	}
	
	q  = 243;
	q = atoi (argv[1]);
	ax = atof (argv[2]);
	ay = atof (argv[3]);
	d = atof (argv[4]);
	e = atof (argv[5]);
	f = atof (argv[6]);
	g = atof (argv[7]);

	printf ("#\n# d=%.2f\n", d);
	printf ("# e=%.2f\n", e);
	printf ("# f=%.2f\n", f);
	printf ("# g=%.2f\n", g);
	printf ("#\n");
	double x = 0.5;
	double y = 0.0;
#if 0
	double aa = sqrt (ax*ax+ay*ay);
	double ama = sqrt ((1.0-ax)*(1.0-ax)+ay*ay);
	double c = sqrt (aa*aa+ama*ama);
	printf ("#\n# denom=%d\n#\n", q);
	printf ("# x=%g y=%g |a| = %g |1-a|=%g c=%g\n#\n", ax,ay, aa, ama,c);

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
#endif

	for (p=0; p<q; p++) 
	{
		double val = (double) p / (double) q;
		x = 0.5;
		y = 0.0;
		fixpt (val, &x, &y);

		printf ("%d\t%g	%g	%g\n", p,val,x,y);
	}
}
