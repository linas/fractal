/*
 * cyli.C
 * Fast and direty hack for the tongue region
 *
 * November 2023
 */

#include <math.h>
#include <stdio.h>

// Perturbed diagonal
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
double alpha(double y, double K)
{
	return root(y, K, 0.0, 0.0, 1.0, 1.0);
}

// Inverse of standard curcle map.
double cinv(double y, double omega, double K)
{
	if (y < omega) return alpha(y-omega+1.0, K);
	return alpha(y-omega, K);
}

// First derivative of pert.
double pprime(double x, double K)
{
	return 1.0 - 2.0 * M_PI * K * cos(2.0 * M_PI * x);
}

// Jacobian of circle map.
double jaco(double y, double omega, double K)
{
	return pprime(cinv(y, omega, K), K);
}

// print jacobian
void dump_jaco(double omega, double K)
{
#define JNPTS 500
	for (int i=0; i< JNPTS; i++)
	{
		double y = (i+0.5) / JNPTS;
		double jc = jaco(y, omega, K);
		printf("%d	%f	%f\n", i, y, jc);
	}
}

// Compute transfer operator applied to fun
double transfer(double y, double (*fun)(double, double, double),
                double omega, double K)
{
	double ci = cinv(y, omega, K);
	double jc = jaco(y, omega, K);
	double tr = fun(ci, omega, K) / jc;
	return tr;
}

// Iterate transfer operator on fun n times.
double triter(double x, int n, double (*fun)(double, double, double),
              double omega, double K)
{
	double ci = cinv(x, omega, K);
	double jc = jaco(x, omega, K);
	double rho = fun(ci, omega, K);
	double y = rho / jc;
	x = ci;
	for (int i=1; i<n; i++)
	{
		double ci = cinv(x, omega, K);
		double jc = jaco(x, omega, K);
		y = y / jc;
		x = ci;
	}
	return y;
}

// Compute average of of transfer operator applied to
// fun one, two, three, ... times. The value returned turns out to
// be equal to the invariant measure.
double taverage(double x, int n, double (*fun)(double, double, double),
                double omega, double K)
{
	double ci = cinv(x, omega, K);
	double jc = jaco(x, omega, K);
	double rho = fun(ci, omega, K);
	double y = rho / jc;
	x = ci;
	double sum = y;
	for (int i=1; i<n; i++)
	{
		double ci = cinv(x, omega, K);
		double jc = jaco(x, omega, K);
		y = y / jc;
		x = ci;

		sum += y;
	}
	return sum / n;
}

double foo(double x, double omega, double K)
{
	return 1.0;
}

double f1(double x, double omega, double K)
{
	return transfer(x, foo, omega, K);
}

double f2(double x, double omega, double K)
{
	return transfer(x, f1, omega, K);
}

double f3(double x, double omega, double K)
{
	return transfer(x, f2, omega, K);
}

double f4(double x, double omega, double K)
{
	return transfer(x, f3, omega, K);
}

double rho(double x, double omega, double K)
{
	return taverage(x, 100, foo, omega, K);
}

double lrho(double x, double omega, double K)
{
	return transfer(x, rho, omega, K);
}

void dump_transfer(double omega, double K, double (*fun)(double, double, double))
{
#define TNPTS 500
	for (int i=0; i< TNPTS; i++)
	{
		double y = (i+0.5) / TNPTS;
		double jc = jaco(y, omega, K);
#if 0
		// Manual recursion
		double t1 = f1(y, omega, K);
		double t2 = f2(y, omega, K);
		double t3 = f3(y, omega, K);
		double t4 = f4(y, omega, K);

		// Double-check loop recursion
		double tr = f4(y, omega, K);
		double tc = triter(y, 4, foo, omega, K);
		printf("yoo %g\n", tr-tc);
#endif

		double t1 = triter(y, 1, foo, omega, K);
		double t2 = triter(y, 100, foo, omega, K);
		double t3 = rho(y, omega, K);
		double t4 = lrho(y, omega, K);
		printf("%d	%f	%f	%f	%f	%f	%f\n", i, y, jc, t1, t2, t3, t4);
	}
}

int main(int argc, char* argv[])
{
	double omega = atof(argv[1]);
	double K = atof(argv[2]);

	// dump_jaco(omega, K);
	dump_transfer(omega, K, foo);
}
