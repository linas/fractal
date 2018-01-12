/*
 * psi.c
 *
 * Orthonormal Hessenberg matrix functions.
 * Construct wave functions which put the downshift operator
 * into Hessenberg form.
 *
 * January 2018
 */
#include <math.h>
#include <stdbool.h>
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

//		printf("%d	%d	%g\n", j, idx, lub);
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
//		printf("%d	%g	%g	%g\n", j, low, midpoints[j], hi);
	}	
}

/* Return the n'th wave function at the value of x. */
double psi_n(double x, double K, int n)
{
	// printf("psi ask for %d x=%g K=%g\n", n, x, K);
	if (0 == n)
	{
		if (K < x) return 0.0;
		return 1.0 / sqrt(K);
	}

	/* Get the lower, middle and upper bounds */
	n++; /* Off-by-one! */
	double lower = midpoints[lower_sequence[n]];
	if (x < lower) return 0.0;
	double upper = midpoints[upper_sequence[n]];
	if (upper < x) return 0.0;

	double middle = midpoints[n];
	double norm = 1.0 / (middle - lower);
	norm += 1.0 / (upper - middle);
	norm = 1.0 / sqrt(norm);

	if (x < middle)
		return norm / (middle - lower);

	return -norm / (upper - middle);
}

/* Compute the matrix elements of the unit operator. These are
 * used to verify orthogonality. Viz no mistakes were made.
 */
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

/* Verify orthogonality */
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

int fapp = 0;
double psi(double x, double K)
{
	return psi_n(x, K, fapp);
}


/* Return matrix entry for the transfer operator in the hessenberg
 * basis.  Brute-force integration. */
double hess_brute(double K, int m, int n)
{
	fapp = n;
// #define IPTS 5311201
#define IPTS 811201
	double s = 0.0;
	for (int i=0; i< IPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) IPTS);
		x *= K;
		double Lfn = xfer(x, K, psi);
		// double Lfn = psi_n(x, K, n);
		double fm = psi_n(x, K, m);
		s += fm*Lfn;
		// printf("wtf %d %d x=%g l=%g r=%g s=%g\n", m, n, x, fm, Lfn, s*K/((double) IPTS));
	}
	s *= K / (double) IPTS;
	return s;
}

/* Return matrix entry for the transfer operator in the hessenberg
 * basis. Using interval math, and thus "perfectly" accurate. */
double hess(double K, int m, int n)
{
	// printf("ask for %d %d\n", m, n);
	if (m == 0 && n == 0) return 1.0; // XXX
	m++; n++;

	// left wavelet
	double mlo = midpoints[lower_sequence[m]];
	double mce = midpoints[m];
	double mhi = midpoints[upper_sequence[m]];

	// right wavelet
	double nlo = midpoints[lower_sequence[n]];
	double nce = midpoints[n];
	double nhi = midpoints[upper_sequence[n]];

	// printf("left= %g %g %g right = %g %g %g\n", mlo, mce, mhi, nlo, nce, nhi);
	// right wavelet in bottom part of xfer oper
	double bnlo = 2.0 * K * nlo;
	double bnce = 2.0 * K * nce;
	double bnhi = 2.0 * K * nhi;

	// right wavelet in top part of xfer oper
	double tnlo = bnlo - K;
	double tnce = bnce - K;
	double tnhi = bnhi - K;

	// printf("bot= %g %g %g top = %g %g %g\n", bnlo, bnce, bnhi, tnlo, tnce, tnhi);
	// Whether or not to skip computing these branches
	bool dobot = true;
	bool dotop = true;
	if (K < bnlo) dobot = false;    // if (nlo < 0.5)
	if (tnhi < 0.0) dotop = false;  // if (nhi < 0.5)
	// printf("dobot=%d dotop=%d\n", dobot, dotop);

	// Clamp to boundaries
	if (K < bnce) bnce = K;
	if (K < bnhi) bnhi = K;

	if (tnlo < 0.0) tnlo = 0.0;
	if (tnce < 0.0) tnce = 0.0;

	// printf("bot= %g %g %g top = %g %g %g\n", bnlo, bnce, bnhi, tnlo, tnce, tnhi);
	if (dobot)
	{
		// If the wavelets don't intersect, then nothing to do
		if (bnhi <= mlo) goto punt;
		if (mhi <= bnlo) goto punt;
		if (mlo <= bnlo && bnhi <= mce) goto punt;
		if (mce <= bnlo && bnhi <= mhi) goto punt;
	}

punt:
	if (dotop)
	{
		// If the wavelets don't intersect, then nothing to do
		if (tnhi <= mlo) return 0.0;
		if (mhi <= tnlo) return 0.0;
		if (mlo <= tnlo && tnhi <= mce) return 0.0;
		if (mce <= tnlo && tnhi <= mhi) return 0.0;
	}

	// If we are here, then there is some non-trivial intersection.
	// This will take some effort to get right.
	// printf("Start hard work\n");

	// Get the wavelet values. Undo the off-by-one
	m--; n--;
	double vmlo = psi_n(0.5*(mlo+mce), K, m);
	double vmhi = psi_n(0.5*(mhi+mce), K, m);

	double vnlo = psi_n(0.5*(nlo+nce), K, n);
	double vnhi = psi_n(0.5*(nhi+nce), K, n);
	// printf("values left= %g %g  right=%g %g\n", vmlo, vmhi, vnlo, vnhi);

	// Accumulate values here.
	double acc = 0.0;
	if (dobot)
	{
		double xa = fmax(bnlo, mlo);
		double xb = fmin(bnce, mce);
		double xc = fmax(bnce, mce);
		double xd = fmin(bnhi, mhi);

		acc += vmlo * vnlo * (xb-xa);
		if (bnce < mce) acc += vmlo * vnhi * (xc - xb);
		else acc += vmhi * vnlo * (xc - xb);
		acc += vmhi * vnhi * (xd-xc);
	}

	if (dotop)
	{
		double xa = fmax(tnlo, mlo);
		double xb = fmin(tnce, mce);
		double xc = fmax(tnce, mce);
		double xd = fmin(tnhi, mhi);

		acc += vmlo * vnlo * (xb-xa);
		if (tnce < mce) acc += vmlo * vnhi * (xc - xb);
		else acc += vmhi * vnlo * (xc - xb);
		acc += vmhi * vnhi * (xd-xc);
	}

	acc /= 2.0 * K;
	return acc;
}

void show_melts(double K)
{
	int mxi = MAXN-1;
	mxi = 5;
// hess(K, 4, 3);
// return;
	for (int i=0; i< mxi; i++)
	{
		int js = i-1;
		if (js < 0) js = 0;
		for (int j=js; j< mxi; j++)
		{
			double g = hess_brute(K, i, j);
			double h = hess(K, i, j);
#if 0
			if (i==j && g < 0.99) printf("Error diag: %d %g", j, g);
			if (i !=j && 5.0e-5 < g) printf("Error off-diag: %d %d %g\n", i, j, g);
#endif
if (g < 2.0e-5) g = 0;
			printf("%d	%d	%g	%g\n", i, j, g, h);
#if 0
if (2e-5 < g) printf("------------------- %g %g %g and %g %g %g\n",
midpoints[lower_sequence[i+1]],
midpoints[i+1],
midpoints[upper_sequence[i+1]],
midpoints[lower_sequence[j+1]],
midpoints[j+1],
midpoints[upper_sequence[j+1]]);
#endif

		}
		printf("# ===========\n");
	}
}

#ifndef NOMAIN
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
	show_melts(K);
}
#endif
