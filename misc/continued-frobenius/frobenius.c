
/*
 * frobenius.c:
 *
 * Iterate on a density distribution with the continued-fraction map.
 * Compare results to the formal Frobenius-Perron operator.
 * Currently best used to obtain numeric approx for the first
 * decaying eigenmode.
 * 
 * A good fit to this first decaying eignemode is provided by 
 * den[i] = -0.75 +1.75*y*y*sqrt(y); with eigenvalue of -0.3025
 * where y = 1/(1+x)
 *
 * Linas Vepstas  <linas@linas.org> Dec 2003
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "zetafn.h"

double
continued_fraction_map (double x) 
{
	double y;
	y = 1/x;
	y -= floor(y);
	return y;
}

double eigenvector_numeric_candidate (double x)
{
	double v[30];

#if OLD_GRAPHICAL_FIT
v[0] =    -0.75 + 1.75/ pow (2.0, 2.5);
v[1] =     0.386699021;
v[2] =    0.3383616433;
v[3] =    0.2537712325;
v[4] =    0.1744677223;
v[5] =    0.1134040195;
v[6] =    0.0708775122;
v[7] =   0.04303277527;
v[8] =   0.02555071031;
v[9] =   0.01490458102;
v[10] =  0.008570134085;
v[11] =  0.004869394366;
v[12] =  0.002739034331;
v[13] =  0.001527538377;
v[14] = 0.0008456016015;
v[15] = 0.0004650808808;
v[16] = 0.0002543411067;
v[17] = 0.0001383914845;
v[18] = 7.496205412e-05;
v[19] = 4.044005551e-05;
v[20] = 2.173652984e-05;
v[21] = 1.164456956e-05;
v[22] =  6.21925874e-06;
v[23] = 3.312431286e-06;
v[24] =  1.75972912e-06;
v[25] = 9.326564338e-07;
v[26] = 4.932317679e-07;
v[27] = 2.603167664e-07;
v[28] = 1.371311537e-07;
v[29] = 7.211207222e-08;
#endif


v[0]=-0.570182;
v[1]=0.414614;
v[2]=0.51254;
v[3]=0.385943;
v[4]=0.246797;
v[5]=0.145288;
v[6]=0.0813875;
v[7]=0.0441307;
v[8]=0.0233975;
v[9]=0.0122086;
v[10]=0.00629728;
v[11]=0.00322096;
v[12]=0.00163734;
v[13]=0.00082857;
v[14]=0.000417907;
v[15]=0.00021027;
v[16]=0.000105611;
v[17]=5.29768e-05;
v[18]=2.65497e-05;
v[19]=1.32967e-05;
v[20]=6.65611e-06;
v[21]=3.3308e-06;
v[22]=1.66637e-06;
v[23]=8.33526e-07;
v[24]=4.16882e-07;
v[25]=2.08483e-07;
v[26]=1.04256e-07;
v[27]=5.21332e-08;
v[28]=2.60684e-08;
v[29]=1.30348e-08;

	int k;
	double acc = 0.0;
	
	for (k=0; k<30; k++)
	{
		acc += v[k];
	}
	
	for (k=0; k<30; k++)
	{
		v[k] /= acc;
	}
	
	double xn = 1.0;
	acc = 0.0;
	x = 1.0-x;
	for (k=0; k<30; k++)
	{
		acc += xn * v[k];
		xn *= x;
	}
	return acc;
}

long double v[30];

void init_funky (void)
{
	int k;
	long double acc = 0.0;
	
	// create basic v_k
	long double sign;
	sign = 1.0L;
	long double vlen = 0.0L;
	for (k=0; k<30; k++)
	{
		v[k] = (2.0/3.0)*(k+1)*zetam1(k+2);
		long double ns=1.0L;
		int n;
		for (n=k; n<60; n++)
		{
			v[k] -= ns*(n+1)*zetam1(n+2);
			ns = -ns;
		}
		vlen += v[k]*v[k];
		
		sign = -sign;
	}
	printf ("# v-zero = %Lg\n", v[0]);
	// norm to unit length
	vlen = 1.0L/sqrtl(vlen);
	printf ("#  renorm=%Lg\n", vlen);
	for (k=0; k<30; k++)
	{
		v[k] *= vlen;
		v[k] -= 1.007181703 * sqrtl(3.0L) *pow (0.5,k+1);
	}
}
	
double funky (double x)
{
	int k;
	double xn = 1.0;
	double acc = 0.0;
	x = 1.0-x;
	for (k=0; k<30; k++)
	{
		acc += xn * v[k];
		xn *= x;
	}
	return acc;
}

long double sk[30];
void
init_skunky(void)
{
	int k,n;
	for (k=0; k<30; k++)
	{
		sk[k] = 0.0;
		long double sn = 1.0L;
		if (k%2==1) sn = -1.0L;
		for (n=k; n<70; n++)
		{
			sk[k] += sn*(n+1)*zetam1(n+2);
			sn = -sn;
		}
	}
}

double skunky (double x)
{
	long double acc = (2.0L/3.0L);
	acc *= trigamma (1+x);
	
	int k;
	double xn = 1.0;
	long double ny = x-1.0;
	for (k=0; k<30; k++)
	{
		acc -= xn * sk[k];
		xn *= ny;
	}

	acc *= 23.85215416;
	acc -= 1.007181703 * sqrtl(3.0L) / (1.0+x);

	return acc;
}

double wunky (double x)
{
	long double acc = (2.0L/3.0L);
	acc *= trigamma (1+x);
	
	int k;
	long double ny = x-1.0;
	double xn = ny;
	long double sn = 1.0;
	for (k=0; k<60; k++)
	{
		long double term = 1.0L - xn;
		term /= (1.0-ny);
		term *= (k+1) * zetam1(k+2);
		acc -= sn * term;
		xn *= ny;
		sn = -sn;
	}

	acc *= 23.85215416;
	acc -= 1.007181703 * sqrtl(3.0L) / (1.0+x);

	return acc;
}

double clunky (double x)
{
	long double acc = (2.0L/3.0L);
	acc *= trigamma (1.0+x);

	acc -= (zetam1(2)-0.25) / (2.0-x);
	
	int k;
	long double y = 1.0-x;
	double xn = 1.0;
	for (k=0; k<60; k++)
	{
		long double term = k+1;
		term *= zetam1(k+2);
		term *= xn;
		term *= y / (1.0+y);
		acc -=term;
		xn *= y;
	}

	acc *= 23.85215416;
	acc -= 1.007181703 * sqrtl(3.0L) / (1.0+x);

	return acc;
}

double junky (double x)
{
	long double acc = (2.0L/3.0L);
	acc *= trigamma (1.0+x);

	acc -= (zetam1(2)-0.25) / (2.0-x);

	long double y = 1.0-x;
	acc -= y*trigamma(2.0-y)/ (1.0+y);
	
	acc *= 23.85215416;
	acc -= 1.007181703 * sqrtl(3.0L) / (1.0+x);

	return acc;
}


double dunky (double x)
{
	double val;
	val = 0.25 - zetam1(2);
	val /= (2.0-x);

	double tmp = 2.0/3.0;
	tmp -= (1.0-x)/(2.0-x);
	val += trigamma (1.0+x) * tmp;

	// renorm and ortho
	val *= 23.85215416;
	val -= 1.007181703 * sqrtl(3.0L) / (1.0+x);

	return val;
}

double * 
init_density (int cnt) 
{
	double *den = (double *) malloc (cnt * sizeof (double));
	
	init_funky();
	init_skunky();
	int i;
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));
		double y = 1.0-x;
		double z = 1.0 / (1.0+x);
		// den[i] = 1.0;
		// den[i] = 3.0*z - 2.0;    // orthogonal to zeroth eigenvalue
		// den[i] = 4.0*y - 3.0*z;     // orthogonal to zeroth eigenvalue
		// den[i] = 1.6*y*y - 0.6*z;     // orthogonal to zeroth eigenvalue
		// den[i] = (8.0**y*y - 3.0*z)/5.0;     // orthogonal to zeroth eigenvalue
		// den[i] = (16.0*y*y*y - 3.0*z)/13.0;     // orthogonal to zeroth eigenvalue
		// den[i] = (32.0*y*y*y*y - 3.0*z)/29.0;     // orthogonal to zeroth eigenvalue
		// den[i] = (64.0*y*y*y*y*y - 3.0*z)/61.0;     // orthogonal to zeroth eigenvalue
		// den[i] = x;  // some arbitary setup 
		// den[i] = z;  // zeroth eigenvector
		// den[i] = 1.0 - z*z;
		// den[i] = z*z;
		// den[i] = -0.75 +1.75*z*z*sqrt(z);
		// den[i] = 1.0-z*z*z*z;
		// den[i] = sin (-14.1 *log(x)) / sqrt(x);
		// den[i] = funky(x);
		// den[i] = skunky(x);
		den[i] = dunky(x);
	}
	return den;
}

void
set_ortho_density (double *den, int cnt, int n) 
{
	double num = 2.0;
	int i;
	for (i=0; i<n; i++)
	{
		num *= 2.0;
	}
	
	for (i=0; i<cnt; i++) 
	{
		double x = ((double) (2*i+1)) / ((double) (2*cnt));
		double y = 1.0-x;
		double z = 1.0 / (1.0+x);
		y = pow (y, n);
		den[i] = (num*y - 3.0*z)/(num-3.0);     // orthogonal to zeroth eigenvalue
	}
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

#if 1
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

#if 1
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
derivative_of_density (double *tgt, double *src, int cnt)
{
	int i;
	for (i=0; i<cnt-1; i++) 
	{
		tgt[i] = (src[i+1]-src[i])*((double) cnt);
	}
	tgt[cnt-1] = tgt[cnt-2];
}

void
rescale_density (double *tgt, double *src, int cnt, double scale)
{
	int i;
	for (i=0; i<cnt; i++) 
	{
		tgt[i] = scale * src[i];
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

int
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

#if GET_FIRST_DECAYING_EIGENVECTOR
	/* remove common mode */
	iterate_perron (den[1], den[0], cnt);
	for (i=0; i<cnt; i++) 
	{
		/* constant calibrated by truncated convergence */
		den[0][i] -= den[1][i] /0.998674 ;
	}
	iterate_perron (den[1], den[0], cnt);
	better_renorm_perron (den[1], den[0], cnt);
	
	int istart = 0;
	for (i=2; i<istart+9; i++) 
	{
		iterate_perron (den[i], den[i-1], cnt);
		// renorm_perron (den[i], den[i-1], cnt);
		better_renorm_perron (den[i], den[i-1], cnt);
	}
#endif

#ifdef JUST_PLAIN_ITERATE
	int istart = 0;
	for (i=1; i<istart+9; i++) 
	{
		 // blur_density (den[i], cnt, 5);
		// iterate_density (den[i], den[i-1], cnt);
		iterate_perron (den[i], den[i-1], cnt);
		// renorm_perron (den[i], den[i-1], cnt);
		// better_renorm_perron (den[i], den[i-1], cnt);
		// second_eigen_perron (den[i], den[i-1], den[i-2], cnt);
	}
#endif
	
#define ITERATE_AND_RESCALE
#ifdef ITERATE_AND_RESCALE
	int istart = 0;
	double acc = 1.0;
	for (i=1; i<istart+9; i++) 
	{
		iterate_perron (den[i], den[i-1], cnt);
		double scale;
		// scale = 1.0/den[i][0];   // non-extrapolation
		//scale = 1.5 * den[i][0] - 0.5 * den[i][1];   // linear extrapolation
		scale = 1.875 * den[i][0] - 1.25 * den[i][1] + 0.375*den[i][2];  // quad
		// scale = 1.875 * den[i][cnt-1] - 1.25 * den[i][cnt-2] + 0.375*den[i][cnt-3];  // quad
		// double pscale = 1.875 * den[i-1][cnt-1] - 1.25 * den[i-1][cnt-2] + 0.375*den[i-1][cnt-3];  // quad
		// scale /= pscale;
		acc *= scale;
		printf ("# i=%d scale=%15.12f acc=%f\n", i, scale, acc);
		
		double correction = 1.0L - 1.0L/(1.202*(double)cnt);
		// scale *= correction;
		scale = 1.0/scale;
		rescale_density (den[i], den[i], cnt, scale);
		fflush (stdout);
	}
#endif


#if SHOW_ORTHO
	int istart = 0;
	for (i=0; i<istart+9; i++) 
	{
		set_ortho_density (den[i], cnt, i);
	}
#endif

#if DERIV
	for (i=0; i<istart+9; i++) 
	{
		derivative_of_density (den[i], den[i+1], cnt);
	}
#endif
	
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
		printf ("%d	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f\n", i, x, 
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
	return 0;
}
