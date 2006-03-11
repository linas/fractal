
/*
 * hermite.C
 *
 * Green's function for the quantum harmonic oscillator
 * K(x,y)= sum_n H_n(x) H_N(y) / w_n   where w_n = n+1/2
 *
 * Linas Vepstas March 2006
 */

/* compute hermite polynomial (physicists norm) up to order
 * nmax, using simple recursive algo, and place results in array "vals"
 */
void hermite_array (double x, int nmax, double *H_n)
{
	H_n[0] = 0.0;
	H_n[1] = 2.0*x;

	int n;
	for (n=1; n<nmax-1; n++)
	{
		H_n[n+1] = 2.0*(x*H_n[n] - n*H_n[n-1]);
	}
}
 
