
/*
 * frobenius.c:
 *
 * Iterate on a density distribution with the continued-fraction map.
 * Compare results to the formal Frobenius-Perron operator.
 */

#include <math.h>

double
continued_fraction_map (double x) 
{
	double y;
	y = 1/x;
	y -= floor(y);
	return y;
}

double * 
init_density (int cnt) 
{
	double *den = (double *) malloc (cnt * sizeof (double));
	
	int i;
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));
		double y = 1.0 / (1.0+x);
		den[i] = 1.0;
		// den[i] = x;  // some arbitary setup 
		// den[i] = y;  // some arbitary setup 
		// den[i] = 1.0 - y*y;
		// den[i] = y*y;
		// den[i] = -0.75 +1.75*y*y*sqrt(y);
		// den[i] = 1.0-y*y*y*y;
	}
	return den;
}


inline double
perron (double x, double *den, int cnt)
{
	/* perform one iteration of frobenius-perron operator. */
   int kc=0,jc=0;
	double val = 0.0;
	int n;
	for (n=1; n<10123; n++)
	{
		double y = 1.0 /(((double) n) + x);
		double dj = y* ((double)cnt);
		// dj -= 0.5; this is wrong!! 
		int j = dj;

		val += den[j] * y * y;
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

void
renorm_perron (double *tgt, double *src, int cnt)
{
	int i;
	/* remove const and renorm, idea is to capture only the decaying
	 * eigenstates. */

	/* integrate to normalize */
	double snorm = 0.0;
	double tnorm = 0.0;
	for (i=0; i<cnt; i++) 
	{
		snorm += src[i];
		tnorm += tgt[i];
	}
	double lambda = tnorm / snorm;

	// double epsilon = 0.1;
	double epsilon = 0.0;
	
	tnorm = 0;
	for (i=0; i<cnt; i++) 
	{
		tgt[i] -= epsilon * src[i];
		tnorm += tgt[i];
	}
	
	double norm = snorm / tnorm;
	
	for (i=0; i<cnt; i++) 
	{
		tgt[i] *= norm;
	}

#if 0
	double zero = 0.5 * (3.0*tgt[0] - tgt[1]);
	zero = 1.0 / zero;

	for (i=0; i<cnt; i++) 
	{
		tgt[i] *= zero;
	}
#endif

	printf ("# duude lambda = %f\n", lambda); 
}

void
better_renorm_perron (double *tgt, double *src, int cnt)
{
	int i;
	/* remove const and renorm, idea is to capture only 
	 * the first decaying eigenstate. */

	/* integrate to normalize */
	double snorm = 0.0;
	double tnorm = 0.0;
	for (i=0; i<cnt; i++) 
	{
		snorm += src[i];
		tgt[i] = src[i] - tgt[i];
		tnorm += tgt[i];
	}
	double lambda = 1.0 - tnorm / snorm;

	double norm = snorm/tnorm;
	
	for (i=0; i<cnt; i++) 
	{
		tgt[i] *= norm;
	}

#if 0
	double zero = 0.5 * (3.0*tgt[0] - tgt[1]);
	zero = 1.0 / zero;

	for (i=0; i<cnt; i++) 
	{
		tgt[i] *= zero;
	}
#endif

	printf ("# duude better lambda = %f\n", lambda); 
}

void
second_eigen_perron (double *f2, double *f1, double *f0, int cnt)
{
	int i;
	/* try to capture the second eigenstate */

	/* integrate to normalize */
	double s0 = 0.0;
	double s1 = 0.0;
	double s2 = 0.0;
	for (i=0; i<cnt; i++) 
	{
		s0 += f0[i];
		s1 += f1[i];
		s2 += f2[i];
	}
	double oeps = 0.998674;
	double lam1 = 0.3025;
	double A = oeps*oeps-lam1*lam1;
	A /= lam1*oeps*(1.0-lam1*oeps);
	double B = -(1.0+A*lam1)/(lam1*lam1);

	double quad = (s0+A*s1+B*s2)/s0;

	double disc = sqrt(A*A-4.0*B*(1.0-quad));
	double lam2p = 0.5*(-A +disc)/B;
	double lam2n = 0.5*(-A -disc)/B;

	printf ("# duuude A=%f B=%f quad=%f disc=%f l+=%f l-=%f\n", A,B,quad,
						 A*A-4.0*B*(1.0-quad), lam2p,lam2n);
	
	/* subtract the 'zero' */
	for (i=0; i<cnt; i++) 
	{
		f2[i] -= f0[i] + A*f1[i] + B*f2[i];
	}

	double lam2 = lam2p;
	
	double norm = lam2*lam2*(1.0 + A*lam2 + B*lam2);
	norm = 1.0/norm;
	for (i=0; i<cnt; i++) 
	{
		f2[i] *= norm;
	}


#if 0
	double zero = 0.5 * (3.0*tgt[0] - tgt[1]);
	zero = 1.0 / zero;

	for (i=0; i<cnt; i++) 
	{
		tgt[i] *= zero;
	}
#endif

	printf ("# duude 2nd lambda = %f\n", lam2); 
}

void
iterate_density (double *tgt, double *src, int cnt)
{
	/* very approximate, numeric iteration, has lots of jitter & aliasing
	 * and general nastiness & numerical artifacts to it */

	int *r = (int *) malloc (cnt * sizeof(int));
	
	int i;
	for (i=0; i<cnt; i++) 
	{
		r[i] = 0;
		tgt[i] = 0.0;
	}
#define NSAMP 357
	int ns = NSAMP * cnt;
	for (i=0; i<ns; i++) 
	{
		double x = ((double) i) / ((double) ns);
		double y = continued_fraction_map (x);
		int j = (int) (y * ((double) cnt));
		if (j >= cnt) { 
			printf ("oops too big! i=%d x=%f y=%f j=%d\n", i,x,y,j);
			j = cnt-1;
		}
		if (j <0) { 
			printf ("oops too small! i=%d x=%f y=%f j=%d\n", i,x,y,j);
			j = 0;
		}
		int ii = i/NSAMP;
		// tgt[ii] += src[j];   // bunches up, repeats (wrong one)
		// tgt[j] += src[ii];   // spreads out & smooths
		tgt[j] += src[ii];   
		r[j] ++;
// printf ("duuude x=%f y=%f j=%d ii=%d tgt=%f r=%d\n",x,y,j,ii,tgt[j],r[j]);
	}
	for (i=0; i<cnt; i++) 
	{
		// tgt[i] /= (double) r[i];
		tgt[i] /= (double) NSAMP;
// printf ("final i=%d tgt=%f r=%d\n",i,tgt[i],r[i]);
	}
}


void
blur_density (double *den, int cnt, int filt)
{
	/* apply blurry linear filter of half-width 'filt' */
	double *tmp = (double *) malloc (cnt * sizeof (double));
	
	double ren = 1.0 / ((double) filt*filt);
	int i;
	for (i=filt; i<cnt-filt; i++ )
	{
		tmp[i] = ((double) filt) * den[i];
		int j;
		for (j=1; j<filt; j++)
		{
			tmp[i] += ((double) (filt-j)) * (den[i+j] + den[i-j]);
		}
		tmp[i] *= ren;
	}
	
	/* lower boundary */
	for (i=0; i<filt; i++) 
	{
		
		tmp[i] = ((double) filt) * den[i];
		int r = filt;
		int j;
		for (j=1; j<filt; j++)
		{
			int k = filt - j;
			tmp[i] += ((double) k) * den[i+j];
			r += k;
			if (i-j >= 0) 
			{
				tmp[i] += ((double) k) * den[i-j];
				r +=k;
			}
		}
		tmp[i] /= (double) r;
	}
	
	/* upper boundary */
	for (i=cnt-filt; i<cnt; i++) 
	{
		
		tmp[i] = ((double) filt) * den[i];
		int r = filt;
		int j;
		for (j=1; j<filt; j++)
		{
			int k = filt - j;
			tmp[i] += ((double) k) * den[i-j];
			r += k;
			if (i+j < cnt) 
			{
				tmp[i] += ((double) k) * den[i+j];
				r +=k;
			}
		}
		tmp[i] /= (double) r;
	}
	
	for (i=0; i<cnt; i++ )
	{
		den[i] = tmp[i];
	}

	free(tmp);
}

void
bin_density (double *den, int cnt, int bin)
{
	double ren = 1.0 / ((double) bin);
	
	den[0] *= 0.5;
	int i, k;
	for (i=0, k=0; i<cnt; i+= bin, k++) 
	{
		if (k) den[k] = 0.0;
		int j;
		for (j=0; j<bin; j++) 
		{
			den[k] += den[i+j];
		}
		den[k] *= ren;
	}
}

void
print_density (double *den, int cnt)
{
	int i;

	for (i=0; i<cnt; i++) 
	{
		double x = ((double) i) / ((double) cnt);
		printf ("%d	%f	%f\n", i, x, den[i]);
	}
}

main (int argc, char * argv[]) 
{
	int cnt = 361;

	if (2 <= argc) 
	{
		cnt = atoi (argv[1]);
	}
	
#define NITER 56
	double *den[NITER];

	int i;
	for (i=0; i<NITER; i++) 
	{
   	den[i] = init_density (cnt);
	}

#if LAMBDA_ONE
	/* remove common mode */
	iterate_perron (den[1], den[0], cnt);
	for (i=0; i<cnt; i++) 
	{
		/* constant calibrated by truncated convergence */
		den[0][i] -= den[1][i] /0.998674 ;
	}
	iterate_perron (den[1], den[0], cnt);
	better_renorm_perron (den[1], den[0], cnt);
#endif

	int istart = 0;
	for (i=1; i<istart+9; i++) 
	{
		 // blur_density (den[i], cnt, 5);
		// iterate_density (den[i], den[i-1], cnt);
		iterate_perron (den[i], den[i-1], cnt);
		// renorm_perron (den[i], den[i-1], cnt);
		// second_eigen_perron (den[i], den[i-1], den[i-2], cnt);
	}

	// iterate_perron (den[2], den[0], cnt);
	
#if WHATTHE
	printf ("#\n");
	printf ("# frobenius.dat\n");
	printf ("# cnt = %d\n", cnt);
	printf ("#\n");
	printf ("#i	x	init	num	frob\n");
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) i) / ((double) cnt);
		printf ("%d	%f	%f	%f	%f\n", i, x, den[0][i],
				den[1][i],
				den[2][i]);
	}
#endif 
	
#if 1
	printf ("#\n");
	printf ("# frobenius.dat\n");
	printf ("# blur = 10\n");
	printf ("# cnt = %d\n", cnt);
	printf ("# istart = %d\n", istart);
	printf ("#\n");
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));
		printf ("%d	%f	%f	%f	%f	%f	%f	%f	%f	%f\n", i, x, 
				den[istart+0][i],
				den[istart+1][i],
				den[istart+2][i],
				den[istart+3][i],
				den[istart+4][i],
				den[istart+5][i],
				den[istart+6][i],
				den[istart+7][i],
				den[istart+8][i]
				);
	}
#endif 
}
