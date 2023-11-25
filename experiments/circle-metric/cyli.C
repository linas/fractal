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

// Compute transfer operator applied to fun.
// It's applied exactly once.
double transfer(double y, double (*fun)(double, double, double),
                double omega, double K)
{
	double ci = cinv(y, omega, K);
	double jc = jaco(y, omega, K);
	double tr = fun(ci, omega, K) / jc;
printf("enter trans with y=%f ci=%f jc=%f tr=%f\n", y, ci, jc, tr);
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
printf("enter triter first y=%f ci=%f jc=%f tr=%f\n", x, ci, jc, y);
	x = ci;
	for (int i=1; i<n; i++)
	{
		double ci = cinv(x, omega, K);
		double jc = jaco(x, omega, K);
		y = y / jc;
printf("enter triter next y=%f ci=%f jc=%f tr=%f\n", x, ci, jc, y);
		x = ci;
	}
	return y;
}

double unit(double x, double omega, double K)
{
	return 1.0;
}

// Compute average of the transfer operator applied to
// fun one, two, three, ... times. Since the transfer operator
// preserves the measure, then when it's applied to anything that
// integrates to 1.0, the returned average will be the invariant
// measure.
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

// Integral of this is exactly zero.
double slope(double x, double omega, double K)
{
	// return x-0.5;
	// return (x>0.5) ? 2.0*(1.0-x) : - 2.0*x;
	// return sin(2.0*M_PI*(x-omega));
	// return sin(2.0*M_PI*(x+omega));
	return cos(2.0*M_PI*(x-omega));
}

// Failed attempt to compute the shift. This won't work because it's not
// actually a shift. Hmm. Let me think ...
double tshift(double x, double shift, int n,
              double omega, double K)
{
	double ci = cinv(x, omega, K);
	double jc = jaco(x, omega, K);
	double rho = slope(ci, omega, K);
	double y = rho / jc;
	x = ci;
	double sum = y;
	double shn = shift;
	double norm = 1.0;
	for (int i=1; i<n; i++)
	{
		double ci = cinv(x, omega, K);
		double jc = jaco(x, omega, K);
		y = y / jc;
		x = ci;

		sum += shn * y;
		norm += shn;
		shn *= shift;
	}
	return sum / norm;
}

double f1(double x, double omega, double K)
{
	// return transfer(x, unit, omega, K);
	return transfer(x, slope, omega, K);
}

// Apply xfer operator exactly once to the previous result.
// For fun=unit, this is identical to triter(x, 2, fun, omega, K);
// and likewise for f3 and f4. But for other funs it is not.
// So WTF am I doing wrong?
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

double mu(double x, double omega, double K)
{
	return taverage(x, 100, unit, omega, K);
}

double trans_mu(double x, double omega, double K)
{
	return transfer(x, mu, omega, K);
}

double half(double x, double omega, double K)
{
	return tshift(x, 30, 0.5, omega, K);
}

double trans_half(double x, double omega, double K)
{
	return transfer(x, half, omega, K);
}

double avslope(double x, double omega, double K)
{
	return taverage(x, 100, slope, omega, K);
}

double trans_slope(double x, double omega, double K)
{
	return transfer(x, avslope, omega, K);
}

void dump_transfer(double omega, double K, double (*fun)(double, double, double))
{
#define TNPTS 500
	for (int i=0; i< TNPTS; i++)
	{
		double y = (i+0.5) / TNPTS;
		double jc = jaco(y, omega, K);

		// Manual recursion
		double t1 = f1(y, omega, K);
		double t2 = f2(y, omega, K);
		double t3 = f3(y, omega, K);
		double t4 = f4(y, omega, K);

#if 0
		// Double-check loop recursion
		double tr = f4(y, omega, K);
		double tc = triter(y, 4, unit, omega, K);
		printf("yoo %g\n", tr-tc);
#endif
		printf("%d	%f	%f	%f	%f	%f	%f\n", i, y, jc, t1, t2, t3, t4);
	}
}

void dump_invariant(double omega, double K)
{
#define INPTS 500
	for (int i=0; i< INPTS; i++)
	{
		double y = (i+0.5) / INPTS;
		double jc = jaco(y, omega, K);

		double t1 = triter(y, 1, unit, omega, K);
#define NAVG 100
		double t2 = triter(y, NAVG, unit, omega, K);
		double t3 = mu(y, omega, K);
		double t4 = trans_mu(y, omega, K);
		printf("%d	%f	%f	%f	%f	%f	%f\n", i, y, jc, t1, t2, t3, t4);
	}
}

double avslo(double x, double omega, double K)
{
	return taverage(x, 50, slope, omega, K);
}

void dump_shift(double omega, double K, double shift)
{
#define SNPTS 500
	for (int i=0; i< SNPTS; i++)
	{
		double y = (i+0.5) / SNPTS;

#if 0
		double t1 = tshift(y, 30, 0.5, omega, K);
		double t2 = trans_half(y, omega, K);
		double t3 = taverage(y, 100, slope, omega, K);
		double t4 = trans_slope(y, omega, K);
#endif
#if 0
		double t1 = taverage(y, 10, slope, omega, K);
		double t2 = taverage(y, 20, slope, omega, K);
		double t3 = taverage(y, 30, slope, omega, K);
		double t4 = taverage(y, 40, slope, omega, K);
#endif
		double t1 = taverage(y, 50, slope, omega, K);
		double t2 = triter(y, 50, slope, omega, K);
		double t3 = triter(y, 51, slope, omega, K);
		double t4 = transfer(y, avslo, omega, K);
		printf("%d	%f	%f	%f	%f	%f\n", i, y, t1, t2, t3, t4);
	}
}

void dump_debug(double omega, double K, double (*fun)(double, double, double))
{
#define TNPTS 500
	for (int i=0; i< TNPTS; i++)
	{
		double y = (i+0.5) / TNPTS;

		// Manual recursion
		double t1 = f1(y, omega, K);
		double t2 = f2(y, omega, K);
		double t3 = triter(y, 2, fun, omega, K);
		// double t3 = f3(y, omega, K);
		double t4 = triter(y, 3, fun, omega, K);

#if 0
		// Double-check loop recursion
		double tr = f4(y, omega, K);
		double tc = triter(y, 4, fun, omega, K);
		printf("yoo %g\n", tr-tc);
#endif
		printf("%d	%f	%f	%f	%f	%f\n", i, y, t1, t2, t3, t4);
	}
}


int main(int argc, char* argv[])
{
	double omega = atof(argv[1]);
	double K = atof(argv[2]);

	// dump_jaco(omega, K);
	// dump_transfer(omega, K, unit);
	// dump_invariant(omega, K);
	// dump_shift(omega, K, 0.5);
	// dump_debug(omega, K, unit);
	// dump_debug(omega, K, slope);

	double y = 0.555;
printf("call f1-----\n");
	f1(y, omega, K);
printf("call f2------\n");
	f2(y, omega, K);
printf("call trit------\n");
	triter(y, 2, slope, omega, K);
}
