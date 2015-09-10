
/**
 * takagi, but for a simple mobius
 *
 * Sept 2015
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double mob_tent(double x)
{
	x -= floor(x);
	if (0.5 < x) return (1.0 - x) / x;
	return x / (1.0 - x);
}

double mob_tak(double x, double alpha)
{
	double sum = 0.0;
	double ap = 1.0;
	double tp = 1.0;
	for (int k=0; k<860; k++)
	{
		sum += ap * mob_tent(tp * x);
		ap *= alpha;
		tp *= 2.0;
		if (ap < 1.0e-16) break;
	}
	return sum;
}

int main(int argc, char * argv[])
{
	double alpha = atof(argv[1]);

	printf ("#\n# alpha=%g\n#\n", alpha);
	for (double x=0.0; x< 1.0; x += 0.002)
	// for (double x=0.49; x< 0.51; x += 0.00002)
	{
		double y = mob_tak(x, alpha);
		printf("%g\t%g\n", x, y);
	}
}

