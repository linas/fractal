/*
 * recheck.c
 * Recheck old (negative) results.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double base(double x)
{
	// return 1.0;
	return 0.5-x;
	// return x -0.5;
	// return sin(2.0* M_PI * x);
	// if (x < 0.8*0.8) return 0.0;
// printf("yo x=%g\n", x);
	// if (x < 0.3333) return 0.0;
	// return 1.0;
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
			sum += obn * base (x);
		}
		if (0.5 < tit) tit -= 0.5;
		tit *= beta;
		obn /= beta;
		// printf("i=%d x=%g tit=%g obn=%g sum=%g\n", i, x, tit, obn, sum);
	}
	return sum;
}

double ell(double beta, double x)
{
	double lo = x / beta;
	double hi = lo + 0.5;

	double ellie = dense(beta, lo) + dense(beta, hi);
	ellie /= beta;
	return ellie;
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
		double elf = ell(beta, x);
		double elg = ell(beta, beta*x);
		printf("%d	%g	%g	%g	%g\n", i, x, f, elf, elg);
	}
}
