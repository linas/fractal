/* 
 * basic dedekind eta function sum 
*/

#include <math.h>

double eta (double ex)
{
	double prod = 1.0;
	double xn = ex;
	int n;
	for (n=1; ; n++)
	{
		double term = 1.0 - xn;
		prod *= term;
		xn *= ex;
		if (1.0e-16 > xn) break;
	}
	prod = prod*prod;
	prod = prod*prod;

	prod *= pow (ex, 1.0/6.0);
	return prod;
}

main ()
{
	double tee =0.5;
	for (tee=0.5; ; tee*= 0.5)
	{
		double s = eta (tee);
		printf ("its %g %g \n", tee, s);
	}
}
