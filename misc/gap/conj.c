
/*
 * conj.c
 *
 * Conjugate map explorer
 *
 * December 2004
 */

double a = 1.0;

double 
phi (double x)
{
	double y;

	y = 0.5 + 0.5*a*(a+1)*(1-2*x)/(x*x-x-a*(a+1));

	return y;
}

double 
phiinv (double y)
{
	double rad;

	rad = 4*a*a*(a+1)*(a+1)/((2*y-1)*(2*y-1));
	rad += 1 + 4*a*(a+1);
	rad = 0.5*sqrt(rad);

	double b = 0.5 - a*(a+1)/(2*y-1);

	if (y<0.5) return b-rad;
	return b+rad;
}
	
int
main()
{
	int i;

	int imax=23;
	for (i=0; i<imax; i++) 
	{
		double x = ((double) i)/((double) imax);

		double y = phiinv (x);
		double xx = phi (y);

		printf ("%d %g  %g  %g\n", i, x, y, xx);
	}
	return 0;
}
