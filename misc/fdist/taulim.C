/*
 * taulim.C
 *
 * Explore limiting behaviour of the GKW eigenfuncs at zero.
 *
 * March 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double triangle(double x)
{
	x -= floor(x);
	if (x < 0.5) return 2.0*x;
	return 2.0-2.0*x;
}

double sum(int k, int n, int lmax)
{
	int l;

	double tk = pow(2.0, k-n);
	double acc = 0.0;
#define LMAX 10123
	for (l=0; l<lmax; l++)
	{
		double tlp1 = 2*l+1;

		double alpha = 1.0;
		acc += alpha * triangle (tlp1*tk) / tlp1;
	}

	acc *= n+1;
	return acc;
}

main(int argc, char * argv[])
{
	int n;
	int k = 2;

	k = atoi(argv[1]);
	int lmax = atoi(argv[2]);

	for (n=1; n<40; n++)
	{
		double y = sum(k,n,lmax);

		printf ("its %d sum=%g\n", n, y);

	}

}
