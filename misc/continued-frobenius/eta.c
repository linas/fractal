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

double sum (double ex)
{
	int n;
	double sum = 0.0;
	for (n=1; ; n++)
	{
		double term = ex + (double) n;
		term = 1.0/ term;
		term = term*term* eta(term);
		sum += term;
	}
	return sum;
}

main ()
{
	double ex;
	for (ex=0.95; ex>0.0; ex -= 0.05)
	{
		double e = eta (ex);
		double s = sum (ex);
		printf ("its x=%g eta=%g sum=%g rat=%g\n", ex, e, s, s/e);
	}
}
