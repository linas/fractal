
/*
 * galois.c
 *
 * solution to galois-theory inspired schroedinger equation
 * January 2006 -- Linas Vepstas
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Return a_n given higer an's and energy */
double 
zero_ansatz (int n, double a8, double a6, double a4, double a2, double e)
{
	double a = a8 *(n*n+15*n+56);
	a += a6*(e -8*n*n -8*13*n -16*21-14-0.75);
	a += a4*(16*n*n + 16*11*n +16*28 + 5 -16*e);
	a += a2*(16*e - 36);

	a /= 16.0;
	return a;
}

double 
gauss_ansatz (int n, double b0, double b2, double b4, double e)
{
	double b = 0.0;

	if (4<=n)
	{
		n += 4;
		b = (2*n-17.0/4.0 - e)*b0;
		b += b2*(-n*n + 4*n -23.0/8.0 + 0.5*e);
		b -= b4*(e/16.0 -0.5*n*n +(3.0*n)/8.0 -31.0/64.0);
		b /= (n*n+3*n+2)/16.0;
	}
	else if (4==n)
	{
		return -((7*16-1-4*e)*b2 + (56+32*e)*b0)/48.0;
	}
	else if (2==n)
	{
		return (31.0/8.0 - 0.5*e)*b0;
	}
	else if (n=0)
	{
		return b0;
	}
	return b;
}

void make_coeffs (double *ar, int nlast, double energy)
{
	int n;
	int nstart = 0;
	if (nlast%2)
	{
		nstart = 1;
	}
	ar[nstart] = 1.0;
	
	n=nstart;
	ar[n+2] = gauss_ansatz (n+2, ar[n], 0.0, 0.0, energy);
	n += 2;
	ar[n+2] = gauss_ansatz (n+2, ar[n-2], ar[n], 0.0, energy);
	for (n=nstart+4; n<=nlast; n+=2)
	{
		// ar[n] = zero_ansatz (n, ar[n+8], ar[n+6], ar[n+4], ar[n+2], energy);
		ar[n+2] = gauss_ansatz (n+2, ar[n-4], ar[n-2], ar[n], energy);
	}

	for (n=nstart; n<nlast; n+=2)
	{
		printf ("%d	%g	%g\n", n, ar[n], ar[n+2]/(ar[n]));
	}
}

/* return value of wave function at x*/
double psi_wf (double x, double * ar, int nlast)
{
	int n;
	double xn = 1.0;
	int nstart = 0;
	if (nlast%2)
	{
		xn = x;
		nstart = 1;
	}
	double psi = 0.0;
	for (n=nstart; n<nlast; n+=2)
	{
		psi += ar[n] * xn;
		// printf ("%d	ar=%g	xn=%g \t term=%g \tpsi=%g\n", n, ar[n], xn, xn*ar[n], psi);
		xn *= x*x;
	}
	// psi *= x*x-0.25;

	psi *= exp (-0.5*x*x);
	return psi;
}

main(int argc, char* argv[])
{
	int n;
#define NMAX 1000
	double ar[NMAX];

	if (argc != 2)
	{
		fprintf (stderr, "Usage: %s <energy>\n", argv[0]);
		exit (1); 
	}
	double energy = atof (argv[1]);

	int nlast = 120;

	make_coeffs (ar, nlast, energy);
	psi_wf (5.0, ar, nlast);

// #define WAVE_FUNCTION
#ifdef WAVE_FUNCTION
	/* wave function */
	double x;

	for (x=0.0; x<7.0; x +=0.4)
	{
		double psi = psi_wf (x, ar, nlast);
		printf ("%g	%g\n", x, psi);
	}
#endif

// #define ASYMPTOTIC
#ifdef ASYMPTOTIC
	double x = 5.0;
	for (energy=0.0; energy>-50; energy-=0.1)
	{
		make_coeffs (ar, nlast, energy);
		double psi = psi_wf (x, ar, nlast);
		printf ("%g	%g\n", energy, psi);
	}
#endif
}

