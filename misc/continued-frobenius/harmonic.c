
/* Harmonic number sums
 *
 * Per erroneous ??? bits in mathowrld
 *
 * August 2005 -- linas
 */

#include <math.h>

double harmonic (int n)
{
	int k;

	if (n<=0) return 0.0;
	
	double acc = 0.0;
	for (k=1; k<=n; k++)
	{
		acc += 1.0 / k;
	}

	return acc;
}

double harmonic_z (double z)
{
	int i;

	double zn=1.0;
	double fac = 1.0;
	double acc = 0.0;
	
	for (i=0; i<30; i++)
	{
		acc += harmonic (i) * zn / fac;
		fac *= i+1;
		zn *= z;
	}

	return acc;
}

double gosper_z (double z)
{
	int i;

	double zn = -z;
	double fac = 1.0;
	double acc = 0.0;
	
	for (i=1; i<30; i++)
	{
		acc += zn / (i*fac);
		fac *= i+1;
		zn *= -z;
	}

	acc *= exp (z);

	return acc;
}


main () {
	int i;

	for (i=0; i<20; i++)
	{
		double z = 0.05*i;

		double h = hamonic_z (z);
		double g = gosper_z (z);
		printf ("its z=%g h=%g  g=%g\n", z,h,g);
	}

}
