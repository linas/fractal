
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
		double x = ((double) i) / ((double) cnt);
		den[i] = x;  // some arbitary setup 
	}
	return den;
}


inline double
perron (double x, double *den, int cnt)
{
	/* perform one iteration of frobenius-perron operator. */

	double val = 0.0;
	int n;
	for (n=1; n<100; n++)
	{
		double y = 1.0 /(((double) n) + x);
		double dj = y* ((double)cnt);
		// dj +=0.5;
		int j = dj;

		val += den[j] * y * y;
		// printf ("n=%d x=%f y=%f val=%f\n", n,x,y,val);
		if (0 == j) break;
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
#define NSAMP 7357
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
	
#define NITER 8
	double *den[NITER];

	int i;
	for (i=0; i<NITER; i++) 
	{
   	den[i] = init_density (cnt);
	}

	for (i=1; i<NITER; i++) 
	// for (i=1; i<2; i++) 
	{
		 // blur_density (den[i], cnt, 5);
		// iterate_density (den[i], den[i-1], cnt);
		iterate_perron (den[i], den[i-1], cnt);
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
	printf ("#\n");
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) i) / ((double) cnt);
		printf ("%d	%f	%f	%f	%f	%f	%f	%f	%f	%f\n", i, x, den[0][i],
				den[1][i],
				den[2][i],
				den[3][i],
				den[4][i],
				den[5][i],
				den[6][i],
				den[7][i],
				den[8][i]
				);
	}
#endif 
}
