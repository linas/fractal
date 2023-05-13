
/*
 * eicheck.c
 *
 * validate eigenvector candidates
 */
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_sf_zeta.h>

#include "zetafn.h"
#include "ache.h"


long double zeta_times_k_alt (int n)
{
	int k;
	long double sign = -1.0L;
	long double acc = 0.0L;
	for (k=n; k<n+60; k++)
	{
		long double term = zetam1 (k+2);
		term *= (long double) (k+1);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	return acc;
}

long double zeta_times_k_alt_p (int n)
{
	int k;
	long double acc = zeta_times_k_alt (n);
	acc +=(2.0L/3.0L) * ((long double) (n+1)) * zetam1(n+2);
	return acc;
}

long double zalt (int n)
{
	int k;
	long double sign = 1.0L;
	long double acc = 0.0L;
	if (n!=0)
	{
	for (k=1; k<n; k++)
	{
		long double term = zetam1 (k+2);
		term *= (long double) (k+1);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	acc -= 0.25L;
	acc *= sign;
	acc += (2.0L/3.0L) * ((long double) (n+1)) * zetam1(n+2);
	}
	else
	{
		acc = 0.25L - zetam1(2)/3.0L;
	}
			  
	return acc;
}

long double whacky (int n)
{
	int k;
	long double acc = 0.0L;
	long double sign = 1.0L;
	for (k=n; k<n+30; k++)
	{
		acc += sign * (k+1)*zetam1(k+2)*pow (0.5,k);
		sign = -sign;
	}

	acc -= 0.6666 * (n+1)*zetam1(n+2);
	return acc;
	
}


int
main (int argc, char *argv[])

{
	int k,m,n,p;

	n=0;
	if (argc==2) {
		n = atoi(argv[1]);
	}

#define WSZ 40
	long double zero[WSZ];
	long double acc = sqrtl(3.0L) / 2.0L;
   for (k=0; k<WSZ; k++)
   {
		zero[k] = acc;
		acc *= 0.5;
	}

	long double para[WSZ];
	long double vec[WSZ];
	long double len = 0.0;
	long double sign = 1.0L;
   for (k=0; k<WSZ; k++)
   {
		// this one works pretty good deviation=0.003
		// vec[k] = -0.5*(k+1)*zetam1(k+2) + (k+2)*zetam1(k+3);
		
		// this one is good to dev=0.0006
		// vec[k] = zeta_times_k_alt_p (k);
		// zalt is an alternate expression for above
		// vec[k]= zalt (k);
		
		vec[k]= whacky (k);
		printf ("duuude initial vec[%d] = %Lg\n", k, vec[k]);

		len += vec[k]*vec[k];
		sign = -sign;
	}

	len = 1.0L/sqrtl(len);
   for (k=0; k<WSZ; k++)
   {
		vec[k] *= len;
		para[k] = -0.3035 * vec[k];
	}

	long double dot = 0.0;
   for (k=0; k<WSZ; k++)
   {
		dot += vec[k] *zero[k];
	}

	long double vplen = 0.0L;
   for (k=0; k<WSZ; k++)
   {
		vec[k] -= dot*zero[k];
		vplen += vec[k]*vec[k];
	}
	vplen = sqrtl (vplen);
	printf (" v dot zero = %Lg\n", dot);
	printf ("length of v_perp = %Lg\n", vplen);

	long double hdot = 0.0;
	long double pdot = 0.0;
	long double hvlen = 0.0L;
	long double hv[WSZ];
	long double hvlast[WSZ];
	for (m=0; m<WSZ; m++)
	{
		long double acc = 0.0;
		long double term;
		for (k=0; k<WSZ; k++)
		{
			term = ache_mp (m, k);
			term *= vec[k];
			acc += term;
		}
		hvlast[m] = term /acc;
		hv[m] = acc;
		hvlen += hv[m]*hv[m];
		hdot += hv[m] *zero[m];
		
		para[m] = hv[m] - para[m];
		pdot += para[m]*zero[m];
		// printf ("duude m=%d parallel = %Lg\n", m, para[m]);
	}
	hvlen = sqrtl (hvlen);
	printf ("length of w=Hv_perp is %Lg  and wlen/v_perp_len=%Lg \n", hvlen, hvlen/vplen);
	printf (" w dot zero = %Lg  and w.zero/wlen=%Lg\n", hdot, hdot/hvlen);

	printf ("pdot = %Lg\n", pdot);

	long double wplen = 0.0;
   for (k=0; k<WSZ; k++)
   {
		hv[k] -= hdot*zero[k];
		para[k] -= pdot*zero[k];
		// printf ("duude m=%d perp of paraallel = %Lg\n", k, para[k]);
		wplen += hv[k]*hv[k];
	}
	wplen = sqrtl(wplen);
	printf ("length of w_perp=%Lg, wplen/vplen=%Lg\n", wplen, wplen/vplen);

	// printf ("dotty=%Lg %Lg\n", dot, hdot/dot);
	//
	// ------------ try again
#ifdef RE_PARALLEL
	long double adj = pdot + dot ;
	adj /= -1.3035;
	printf ("duude adj=%20.15Lg\n", adj);

	len = 0.0;
   for (k=0; k<WSZ; k++)
   {
		vec[k] = zeta_times_k_alt_p (k);
		len += vec[k]*vec[k];
	}
	len = 1.0L/sqrtl(len);
   for (k=0; k<WSZ; k++)
   {
		vec[k] *= len;
	}

	len = 0.0;
   for (k=0; k<WSZ; k++)
   {
		vec[k] += adj*zero[k];
		len += vec[k]*vec[k];
	}
	len = 1.0L/sqrtl(len);
   for (k=0; k<WSZ; k++)
   {
		vec[k] *= len;
	}

	hvlen = 0.0;
	for (m=0; m<WSZ; m++)
	{
		long double acc = 0.0;
		long double term;
		for (k=0; k<WSZ; k++)
		{
			term = ache_mp (m, k);
			term *= vec[k];
			acc += term;
		}
		hv[m] = acc;
		hvlen += hv[m]*hv[m];
	}
	hvlen = sqrtl (hvlen);
	printf ("duude hvlen = %Lg\n", hvlen);

#endif

	long double avg = 0.0;
	for (m=0; m<WSZ; m++)
	{
		avg += hv[m] / vec[m];
	}
	avg /= WSZ;
	
	long double dev = 0.0;
	for (m=0; m<WSZ; m++)
	{
		long double term = avg - hv[m] / vec[m];
		dev += term*term;
	}
	dev /= WSZ;
	dev = sqrt (dev);
	
	printf ("avg=%Lg dev = %Lg\n", avg, dev);

	for (m=0; m<WSZ; m++)
	{
		acc = hv[m];
		long double val = vec[m];
		long double lambda = acc/val;
		printf ("%d v=%Lg  Hv=%Lg  hvlast=%Lg    lam=%Lg  d=%Lg\n", 
							 m, val, acc, hvlast[m], lambda, lambda-avg);
	}
	
	return 0;
}
