
/*
 * Recurrance in the transer operator.
 *
 * Dec 2017 Linas Vepstas
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double xiter (double x, double K, int lvl)

	if (K < x) return 0.0;
	if (lvl < 0)
	{
		if (x < 0.5*K) return 1.0;
		return -1.0;
	}
	double sum = 0.0;
	double otk = 0.5 / K;
	double xtk = otk*x;
	sum += xiter(xtk);
	sum += xiter(xtk+0.5);
	sum *= otk;
	return sum;
}

int main(int argc, char* argv[])
{
	double K = atof(argv[1]);
#define NPTS 600
#define NREC 15

	for (int i=0; i< NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		printf ("%d	%g", i, x);
		for (int r=0; r<NREC; r++)
		{
			double y = xiter(x, K, r);
			printf("	%g", y);
		}
		printf ("\n");
	}
}
