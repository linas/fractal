
/* 
 * valuate.c
 *
 * Evaluate the frobenius perron operator on explicit expressions
 * for the eigenfunctions.   Bzzzt broken right now.
 *
 * Linas December 2003
 */

#include <stdio.h>

/* return psi'(1+x)  */
inline double
psi_prime(double x)
{
	int n;
	double acc = 0.0;
	for (n=1; n<10123; n++)
	{
		x += 1.0;
		acc += 1.0 / (x*x);
	}	
	return acc;
}

/* same as above, more accuracy */
inline double
better_psi_prime(double x)
{
	int n;
	double acc = 0.0;
	double en = 0.0;
	double npx = x;
	for (n=1; n<3123; n++)
	{
		npx += 1.0;
		en += 1.0;
		acc -= x*(x+2.0*en) / (en*en*npx*npx);
	}	

	acc += 1.644934066848226;
	return acc;
}

/* return value for first eigenvector */
inline double
eigenvect (double x)
{
	// return  psi_prime(x);
	// return  better_psi_prime(x)/ 1.644934066848226;
	return 1.0/(1.0+x);
}

/* appy frobenius perron on eigenvect at x */
double
frob (double x)
{
	double acc = 0.0;
	int n;
	for (n=1; n<20123; n++)
	{
		x += 1.0;
		double y = 1.0 / x;
		acc += y* y* eigenvect (y);
		// if (0 == n%100) { printf ("."); fflush (stdout); }
	}
	return acc;
}


double meric (double x)
{
	double acc = 0.0;

	acc = 0.954966;
  	acc += -0.274879 *x;
  	acc += 0.112234 *x*x;
	acc += -0.029104 *x*x*x;
	acc += 0.008093 *x*x*x*x;
	acc += -0.002180 *x*x*x*x*x;
  	acc += 0.000585 *x*x*x*x*x*x;
	acc += -0.000156 *x*x*x*x*x*x*x;
	acc += 0.000041 *x*x*x*x*x*x*x*x;
	acc += -0.000011 *x*x*x*x*x*x*x*x*x;
	
	return acc;
}
  
main ()
{
	int i;

	printf ("# \n");
	printf ("# eigenvactor evaluation \n");
	printf ("# \n");
	printf ("# i	x	eig	frob(eig) f/e\n");
	printf ("# \n");
	fflush (stdout);

#define NPTS 30
	for (i=0; i<NPTS; i++)
	{
		double x = ((double) i) / ((double) NPTS);
		double e = eigenvect (x);
		// double e = meric (x);
		double f = frob(x);
		printf ("%d	%f	%f	%f	%f\n", i,x,e,f, f/e);
		fflush (stdout);
	}
}
