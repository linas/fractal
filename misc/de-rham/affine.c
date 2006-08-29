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

main (int argc, char *argv[])
{
	int i;
	int p,q;

	if (8 > argc)
	{
		fprintf (stderr, "Usage: %s <q> <ax> <ay>\n", argv[0]);
		exit (1);
	}
	
	double ax, ay, d,e,f,g;
	q  = 243;
	q = atoi (argv[1]);
	ax = atof (argv[2]);
	ay = atof (argv[3]);
	d = atof (argv[4]);
	e = atof (argv[5]);
	f = atof (argv[6]);
	g = atof (argv[7]);

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
	
	printf ("#\n# d=%.2f\n", d);
	printf ("# e=%.2f\n", e);
	printf ("# f=%.2f\n", f);
	printf ("# g=%.2f\n", g);
	printf ("#\n");
	double x = 0.5;
	double y = 0.0;
	for (p=0; p<q; p++) 
	{
		double val = (double) p / (double) q;
		fixpt (val);

#if 0
		printf ("p=%d\tval=%g\n", p,val);
		printf ("%6.3g\t%6.3g\t%6.3g\n", result[0][0], result[0][1], result[0][2]);
		printf ("%6.3g\t%6.3g\t%6.3g\n", result[1][0], result[1][1], result[1][2]);
		printf ("----------------------------\n\n");
#endif
		x = result[0][2];
		y = result[1][2];
		printf ("%d\t%g	%g	%g\n", p,val,x,y);
	}
}
