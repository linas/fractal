
/*
 * zeta-fourier.c
 *
 * Fourier transforms of zeta-function-zero based propator.
 *
 * Linas September 2008
 */

double *zeros;

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

	coeffs();
}
