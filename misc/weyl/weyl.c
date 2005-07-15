/* 
 * weyl.c
 * 
 * Weyl group functions hack
 *
 * Linas Vepstas July 2005
 */

#include <math.h>

double 
olap (double x, int n, int mterms)
{
	double re_acc = 0.0;
	double im_acc = 0.0;

	int m;
	
	int mmax = 4*mterms;

	double re_ex = cos (2.0*M_PI*x);
	double im_ex = sin (2.0*M_PI*x);

	/* add terms 4 at a time (one is zero) */
	for (m=-mmax-4*n; m<mmax-4*n; m+=4)
	{
		double re_emx = cos (2.0*M_PI*(m+2)*x);
		double im_emx = sin (2.0*M_PI*(m+2)*x);

		double re_term = re_ex / ((double)m+1+4*n);
		double im_term = -im_ex / ((double)m+1+4*n);
					 
		re_term -= 2.0 / ((double)m+2+4*n);
		im_term -= 2.0 / ((double)m+2+4*n);
		re_term += re_ex / ((double)m+3+4*n);
		im_term += im_ex / ((double)m+3+4*n);

		re_acc += re_term*re_emx - im_term*im_emx;
		im_acc += im_term*re_emx + re_term*im_emx;
	}

	re_acc /= M_PI;
	im_acc /= M_PI;

	double tmp = re_acc;
	re_acc = im_acc;
	im_acc = -re_acc;

	return re_acc;
}

main ()
{
	int i;

	int nmax = 723;

	for (i=0; i<nmax; i++)
	{
		double x = (double) i / ((double) nmax);

		double y = olap (x,1,5);

		printf ("%d	%g	%g\n", i, x,y);
	}
}
