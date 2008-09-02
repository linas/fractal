
/*
 * zeta-fourier.c
 *
 * Fourier transforms of zeta-function-zero based propator.
 *
 * Linas September 2008
 */

#include <math.h>
#include <stdio.h>

double *zeros;

void read_zeros (int nzeros)
{
	int i;
	zeros = (double *) malloc (nzeros * sizeof(double));
	FILE * fh = fopen ("zeros1");
	for (i=0; i<nzeros; i++)
	{
		char buff[100];
		fgets(buff, 100, fh);
		zeros[i] = atof(buff);
printf ("its %g\n", zeros[i]);
	}
	fclose (fh);
}

void coeffs(int nmax)
{
	int n;

	double *bn = (double *) malloc (nmax * sizeof(double));
	for (n=0; n<nmax; n++)
	{
		bn[n] = 0.0;
		for (i=0; i<nzeros; i++)
		{
			b[n] += 1.0 / (n*n - zeros[i]*zeros[i]);
		}

		b[n] *= 2.0*n;

		printf ("%d\t%8.6g\n", n, b[n]);
	}
}


main()
{

	read_zeros(30);
	coeffs();
}
