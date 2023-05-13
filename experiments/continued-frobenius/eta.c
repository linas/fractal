/* 
 * basic dedekind eta function sum 
*/

#include <math.h>

long double eta (long double ex)
{
	long double prod = 1.0L;
	long double xn = ex;
	int n;
	for (n=1; ; n++)
	{
		long double term = 1.0L - xn;
		prod *= term;
		xn *= ex;
		if (1.0e-20 > xn) break;
	}
	prod = prod*prod;
	prod = prod*prod;

	prod *= powl (ex, 1.0L/6.0L);
	return prod;
}

long double sum (long double ex)
{
	int n;
	long double sum = 0.0;
	for (n=1; ; n++)
	{
		long double term = ex + (long double) n;
		term = 1.0L/ term;
		term = term*term* eta(term);
		if (n>50 && 1.0e-12 > term) break;
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
