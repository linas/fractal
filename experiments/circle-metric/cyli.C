/*
 * cyli.C
 * Fast and direty hack for the tongue region
 *
 * November 2023
 */

#include <math.h>
#include <stdio.h>

double pert(double x, double K)
{
	return x - K * sin(2.0 * M_PI * x); 
}

// Terrible stupid root finder; brute force subdivision
double root (double y, double K, double xlo, double ylo, double xhi, double yhi)
{
#define EPS 1.0e-15
	double xmid = 0.5 * (xlo+xhi);
	if (abs(xhi-xlo) < EPS) return xmid;

	double ymid = pert(xmid, K);
	if (y<ymid) return root(y, K, xlo, ylo, xmid, ymid);

	return root(y, K, xmid, ymid, xhi, yhi);
}

// Inverse of pert.
double alpha (double y, double K)
{
	return root(y, K, 0.0, 0.0, 1.0, 1.0);
}

double cinv(double y, double K, double omega)
{
	if (y < omega) return alpha(y-omega+1.0, K);
	return alpha(y-omega, K);
}

double pprime(double x, double K)
{
	return 1.0 - 2.0 * M_PI * K * cos(2.0 * M_PI * x);
}

double jaco(double y, double K, double omega)
{
	return pprime(cinv(y, K, omega), K);
}

int main(int argc, char* argv[])
{
	double omega = atof(argv[1]);
	double K = atof(argv[2]);

#define NPTS 500
	for (int i=0; i< NPTS; i++)
	{
		double y = (i+0.5) / NPTS;
		double jc = jaco(y, K, omega);
		printf("%d	%f	%f\n", i, y, jc);
	}
}
