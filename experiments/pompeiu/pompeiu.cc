
/*
 * One very specific Pompieu function, just for grins.
 *
 * Linas Vepstas September 2015
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double gee(double x)
{
	double sum= 0.0;
	for (int lev=0; lev<14; lev++)
	{
		unsigned long deno = 1 << lev;
		double indo = 1.0 / (double) deno;
		for (unsigned long num = 1; num<deno; num+=2)
		{
			double qi = num * indo;
			double term = indo * pow(fabs(x-qi), 1.0L/3.0L);
			if (x < qi) term = -term;
			sum += term;
		}
	}
	return sum;
}

int main()
{
	for (double x=0.0; x<1.0; x += 0.00111234)
	{
		printf("%g	%g\n", x, gee(x));
	}
}
