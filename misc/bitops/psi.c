/*
 * psi.c
 *
 * Attempt at orthonormal fns.
 * January 2018
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double psi_0(double x, double K)
{
	if (x<K) return 1.0/sqrt(K);
	return 0.0;
}

double psi_1(double x, double K)
{
	double step = 2.0*K*K-K;
	if (x<step) return 1.0/sqrt(2.0*step);
	if (x<K) return -1.0/sqrt(2.0*(K-step));
	return 0.0;
}

/* Transfer operator applied to function */
double xfer(double x, double K, double (*fun)(double, double))
{
	if (K<x) return 0.0;
	double res = x / (2.0*K);
	double elf = fun(res, K) + fun(0.5+res, K); 
	elf /= 2.0*K;
	return elf;
}

/* One or the other part of the transfer operator applied to function */
double part(double x, double K, double (*fun)(double, double), int which)
{
	if (K<x) return 0.0;
	double res = x / (2.0*K);
	double elf = fun(0.5*which+res, K); 
	elf /= 2.0*K;
	return elf;
}

#define MAXN 50
double midpoints[MAXN];
int mid_sequence[MAXN];

void find_midpoints(double K)
{
	/* First, find the midpoints; the numbering here is off-by-one:
	 * the "first" midpoint is in midpoints[2].
	 */
	midpoints[0] = 0.0;
	midpoints[1] = K;

	double m = K;
	for (int i=2; i< MAXN; i++)
	{
		if (m <= 0.5) m = 2.0*K*m;
		else m = 2.0*K*m - K;
		midpoints[i] = m;
	}

	/* Now sort them in sequential order. Use a super-stupid sort algo. */
	mid_sequence[0] = 0;

	for (int j=1; j< MAXN; j++)
	{
		double m = mid_sequence[j-1];
		double lst = K;
		int idx = MAXN;
		for (int i=1; i< MAXN; i++)
		{
			if (m < midpoints[i] && midpoints[i] < lst) { lst = midpoints[i]; idx = i; }
		}
		mid_sequence[j] = idx;

		printf("%d	%d	%g\n", j, idx, lst);
	}
}


int main(int argc, char* argv[])
{
	if (argc < 1)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	printf("#\n# K=%g\n#\n", K);

#if 0
#define NPTS 201
	double s = 0.0;
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double (*f)(double, double) = psi_1;
		double ef = f(x, K);
		double p0 = part(x, K, f, 0);
		double p1 = part(x, K, f, 1);
		double y = xfer(x, K, f);
		s += y / ((double) NPTS);
		printf("%d	%g	%g %g	%g	%g	%g\n", i, x, ef, p0, p1, y, s);
	}
#endif

	find_midpoints(K);

}
