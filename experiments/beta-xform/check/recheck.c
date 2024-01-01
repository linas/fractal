/*
 * recheck.c
 * Recheck old (negative) results.
 *
 * January 2024
 */

#include <stdio.h>
#include <stdlib.h>

double base(double x)
{
	return 1.0;
}

double dense(double beta, double x)
{
	double tit = 0.5*beta;
	double obn = 1.0;
	double sum = 0.0;
	for (int i=0; i<24; i++)
	{
		if (x < tit)
		{
			sum += obn * base (x/beta);
		}
		if (0.5 < tit) tit -= 0.5;
		tit *= beta;
		obn /= beta;
		// printf("i=%d x=%g tit=%g obn=%g sum=%g\n", i, x, tit, obn, sum);
	}
	return sum;
}

int main(int argc, char* argv[])
{
	double beta = atof(argv[1]);
	int npts = 1000;
	printf("#\n# beta=%g\n#\n", beta);
	for (int i=0; i< npts; i++)
	{
		double x = (((double) i) + 0.5) / ((double) npts);
		double f = dense(beta, x);
		printf("%d	%g	%g\n", i, x, f);
	}
}
