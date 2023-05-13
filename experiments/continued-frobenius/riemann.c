
/*
 * riemann.c:
 *
 * Plot riemann zeta version of continued fraction frobenius perron operator
 * (for values along imaginary s axis)
 * 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "zetafn.h"

void
init_density_imag (double *re, double *im, int cnt, double ess) 
{
	int i;
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));
		double y = 1.0 / (1.0+x);
		double sq = 1.0/sqrt(x);
		re[i] = cos (ess*log(x)) *sq;
		im[i] = sin (ess*log(x)) *sq;
	}
}

void
init_density_ess (double *re, int cnt, double ess) 
{
	printf ("ess= %g\n", ess);
	ess -= 1.0;
	int i;
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));
		double y = 1.0 / (1.0+x);
		re[i] = pow (x,ess);
	}
}

void
init_density_why (double *re, int cnt, int p) 
{
	int i;
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));
		double ess = p;
		re[i] = pow (1.0-x,ess);
	}
}

long double
ache_mp(int m, int p)
{
	int k;

	long double acc = 0.0L;
	long double sign = 1.0L;
	for (k=0; k<=p; k++)
	{
		long double term = zetam1 (k+m+2);
		term *= binomial (m+k+1,m);
		term *= binomial (p,k);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	return acc;
}

void
init_density_ucfy (double *re, int cnt, int p) 
{
	int i;
	for (i=0; i<cnt; i++) 
	{
		long double x = ((double) (2*i+1)) / ((double) (2*cnt));
		long double y = 1.0-x;
		int m;
		long double ym = 1.0L;
		long double acc = 0.0L;
		for (m=0; m<30; m++)
		{
			acc += ym * ache_mp (m,p);
			ym *= y;
		}
		re[i] = acc;
	}
}


inline double
perron (double x, double *den, int cnt)
{
	/* perform one iteration of frobenius-perron operator. */
   int kc=0,jc=0;
	double val = 0.0;
	int n;
	for (n=1; n<400123; n++)
	{
		double y = 1.0 /(((double) n) + x);

#if STRAIGHT_SAMPLE
		// sample denisty at the nearest sample point
		// this potenitally has some aliasing artifacts
		double dj = y* ((double)cnt);
		// dj -= 0.5; this is wrong!! 
		int j = dj;
		val += den[j] * y * y;
#endif

#define LINEAR_SAMPLE 1
#if LINEAR_SAMPLE
		// sample denisty by interpolating linearly between 
		// two nearest points. This should be pretty damn accurate.
		double dj = y* ((double)cnt) - 0.5;
		int j = dj;
		double frac =  dj -j;
		if (cnt-1 == j) { j--; frac += 1.0; }
		double samp = (1.0-frac) * den[j] + frac * den[j+1];
		val += samp * y * y;
#endif
		if (1==j) kc++;
		// printf ("n=%d j=%d x=%f y=%f val=%f\n", n,j,x,y,val);
		if (0 == j) 
		{
			jc ++;
			if (jc > 4*kc) break;
		}
	}

	return val;
}

void
iterate_perron (double *tgt, double *src, int cnt)
{
	/* one iteration of frobenius-perron, 'exactly', subject to numerics */
	int i;
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));
		tgt[i] = perron (x, src, cnt);
	}
}

double
integrate_density (double *src, int cnt)
{
	/* Riemann integral */
	int i;
	double acc = 0.0;
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));
		acc += x*src[i];
	}
	acc /= (double) cnt;
	return acc;
}

void
do_imag_axis (int cnt)
{
	double *rexs = (double *) malloc (cnt*sizeof(double));
	double *imxs = (double *) malloc (cnt*sizeof(double));
	double *reucfxs = (double *) malloc (cnt*sizeof(double));
	double *imucfxs = (double *) malloc (cnt*sizeof(double));
	
	printf ("#\n");
	printf ("# riemann.dat\n");
	printf ("# reimann zeta drived via continued frac\n");
	printf ("# cnt = %d\n", cnt);
	printf ("#\n");
	printf ("#ess	re	im\n");
	printf ("#\n");
	
	double ess;
	for (ess=7.0; ess<50; ess+=0.2) 
	{
		init_density_imag (rexs, imxs, cnt, ess);
		iterate_perron (reucfxs, rexs, cnt);
		iterate_perron (imucfxs, imxs, cnt);
		double reval = integrate_density (reucfxs, cnt);
		double imval = integrate_density (imucfxs, cnt);
		double re = 0.5*reval - ess*imval;
		double im = 0.5*imval + ess*reval;
		reval = -(4.0-ess*ess) / (4.0+ess*ess);
		imval = -4.0*ess / (4.0+ess*ess);
		re = reval - re;
		im = imval - im;

		printf ("%g	%g	%g\n", ess, re, im);
		fflush (stdout);
	}
}

double
do_real_zeta (double ess, int cnt)
{
	double *src = (double *) malloc (cnt*sizeof(double));
	double *tgt = (double *) malloc (cnt*sizeof(double));
	init_density_ess (src, cnt, ess);
	iterate_perron (tgt, src, cnt);
	double val = integrate_density (tgt, cnt);
	val = -val;
	val += 1.0/(ess-1.0);
	val *= ess;
	return val;
}

void 
do_why (int p, double cnt)
{
	double *src = (double *) malloc (cnt*sizeof(double));
	double *tgt = (double *) malloc (cnt*sizeof(double));
	init_density_why (src,cnt,p);
	iterate_perron (tgt, src, cnt);

	init_density_ucfy (src,cnt,p);

	int i;
	for (i=0; i<cnt; i++)
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));

		printf ("%d\t%g\t%g\t%g\t%g\n", i,x,tgt[i],src[i], tgt[i]-src[i]);
	}

}

void 
do_whyp (int cnt)
{
	double *arr[10];
	int p;
	int i;

	for (p=0; p<10; p++)
	{
		double *src = (double *) malloc (cnt*sizeof(double));
		init_density_ucfy (src,cnt,p);
		arr[p] = src;
	}

	printf ("#\n");
	printf ("# why powers\n");
	printf ("#\n");
	printf ("# U_cf acting on powers of y=1-x\n");
	printf ("# columns are \n");
	printf ("# x, 0,1,2,3,4\n");
	printf ("#\n");
	printf ("#\n");
	for (i=0; i<cnt; i++)
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));

		printf ("%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
			i,x,
			arr[0][i],
		   arr[1][i],
		   arr[2][i],
		   arr[3][i],
		   arr[4][i],
		   arr[5][i],
		   arr[6][i],
		   arr[7][i],
		   arr[8][i],
		   arr[9][i]);
	}

}

int
main (int argc, char * argv[]) 
{
	int cnt = 361;

	if (2 <= argc) 
	{
		cnt = atoi (argv[1]);
	}
	
#ifdef CHECK_BASIC_UCF_EQN
	double ess = 3.0;
	double val = do_real_zeta (ess, cnt);
	printf ("its %12.8g\n", val);
#endif

#ifdef CHECK_ACHE
	do_why (0,cnt);
#endif
	do_whyp (cnt);
	return 0;
}
