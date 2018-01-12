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
int lower_sequence[MAXN];
int upper_sequence[MAXN];

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
		double mid = midpoints[mid_sequence[j-1]];
		double lub = K;
		int idx = MAXN;
		for (int i=1; i< MAXN; i++)
		{
			if (mid < midpoints[i] && midpoints[i] < lub) {
				lub = midpoints[i]; idx = i;
			}
		}
		mid_sequence[j] = idx;

		printf("%d	%d	%g\n", j, idx, lub);
	}

	printf("#\n#------------------------\n#\n");

	/* Compute the pointers to the lower and the upper ends */
	lower_sequence[0] = 0;
	lower_sequence[1] = 0;
	upper_sequence[0] = 1;
	upper_sequence[1] = 1;
	for (int j=2; j< MAXN; j++)
	{
		double low = 0.0;
		double hi = K;
		int lidx = 0;
		int uidx = 1;
		for (int i=1; i< j; i++)
		{
			if (low < midpoints[i] && midpoints[i] < midpoints[j])
			{
				low = midpoints[i];
				lidx = i;
			}
			if (midpoints[i] < hi && midpoints[j] < midpoints[i])
			{
				hi = midpoints[i];
				uidx = i;
			}
		}
		lower_sequence[j] = lidx;
		upper_sequence[j] = uidx;
		printf("%d	%g	%g	%g\n", j, low, midpoints[j], hi);
	}	
}

/* Return the n'th wave function at the value of x. */
double psi_n(double x, double K, int n)
{
	/* Get the lower, middle and upper bounds */
	n++; /* Off-by-one! */
	double lower = midpoints[lower_sequence[n]];
	if (x < lower) return 0.0;
	double upper = midpoints[upper_sequence[n]];
	if (upper < x) return 0.0;
	double middle = midpoints[n];
	if (x < middle)
		return 1.0 / sqrt(2.0 * (middle - lower));

	return 1.0 / sqrt(2.0 * (upper - middle));
}

/* Verify orthogonality */
double psi_prod(int m, int n)
{
	if (m == 0 && n == 0) return 1.0;
	m++; n++;

	double mlo = midpoints[lower_sequence[m]];
	double nhi = midpoints[upper_sequence[n]];
	if (nhi <= mlo) return 0.0;

	double mhi = midpoints[upper_sequence[m]];
	double nlo = midpoints[lower_sequence[n]];
	if (mhi <= nlo) return 0.0;

	double mc = midpoints[m];
	if (mlo <= nlo && nhi <= mc) return 0.0;
	if (mc <= nlo && nhi <= mhi) return 0.0;

	double nc = midpoints[n];
	if (nlo <= mlo && mhi <= nc) return 0.0;
	if (nc <= mlo && mhi <= nhi) return 0.0;

	return 1.0;
}

void verify_ortho(void)
{
	for (int j=0; j< MAXN-1; j++)
	{
		for (int i=0; i< MAXN-1; i++)
		{
			double d = psi_prod(i, j);
			if (i != j && 0.0 < d) printf("Error: off-diag %d %d\n", i, j);
			if (i == j && d < 1.0) printf("Error: not diag %d %g\n", i, d);
		}
	}
	printf("Done verifying orthogonality\n");
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
	verify_ortho();
}
