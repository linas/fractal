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
	return tr;
}

// Iterate transfer operator on fun n times.
double triter(double x, int n, double (*fun)(double, double, double),
              double omega, double K)
{
	double jacobian = 1.0;
	for (int i=0; i<n; i++)
	{
		double ci = cinv(x, omega, K);
		double jc = pprime(ci, K);
		jacobian *= jc;
		x = ci;
	}
	// Use the x coming out of the loop above. It's cinv iterated.
	double rho = fun(x, omega, K);
	double y = rho / jacobian;
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
	double sum = 0.0;
	double jacobian = 1.0;
	for (int i=0; i<n; i++)
	{
		double ci = cinv(x, omega, K);
		double jc = pprime(ci, K);
		jacobian *= jc;
		x = ci;

		double rho = fun(x, omega, K);
		double y = rho / jacobian;
		sum += y;
	}
	return sum / n;
}

// Integral of this is exactly zero.
double slope(double x, double omega, double K)
{
return -1.0;
	// return (x>0.5) ? 2.0:0.0;
	// return x-0.5;
	// return (x>0.5) ? 1.0:-1.0;
	// return (x>0.5) ? 2.0*(1.0-x) : - 2.0*x;
	// return sin(2.0*M_PI*(x-omega));
	// return sin(2.0*M_PI*(x+omega));
	// return cos(2.0*M_PI*(x-omega));
}

// Failed attempt to compute the shift. This won't work because it's not
// actually a shift. Hmm. Let me think ...
double tshift(double x, int n, double shift,
              double (*fun)(double, double, double),
              double omega, double K)
{
	double sum = 0.0;
	double shn = 1.0;
	double norm = 0.0;
	double jacobian = 1.0;
	for (int i=0; i<n; i++)
	{
		double ci = cinv(x, omega, K);
		double jc = pprime(ci, K);
		jacobian *= jc;
		x = ci;

		double rho = fun(x, omega, K);
		double y = rho / jacobian;

		sum += shn * y;
		norm += shn;
		shn *= shift;
	}
	return sum / norm;
}

double f1(double x, double omega, double K)
{
	return transfer(x, unit, omega, K);
	// return transfer(x, slope, omega, K);
}

// Apply xfer operator exactly once to the previous result.
// This is identical to triter(x, 2, fun, omega, K);
// and likewise for f3 and f4. Used only for testing.
double f2(double x, double omega, double K)
{
	double rv = transfer(x, f1, omega, K);
	return rv;
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
#define NAVG 200
	return taverage(x, NAVG, unit, omega, K);
}

double trans_mu(double x, double omega, double K)
{
	return transfer(x, mu, omega, K);
}

double half(double x, double omega, double K)
{
	return tshift(x, 30, 0.5, slope, omega, K);
}

double trans_half(double x, double omega, double K)
{
	return transfer(x, half, omega, K);
}

double avslope(double x, double omega, double K)
{
	return taverage(x, 30, slope, omega, K);
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
#define NSTEP 100
		double t2 = triter(y, NSTEP, unit, omega, K);
		double t3 = mu(y, omega, K);
		double t4 = trans_mu(y, omega, K);
		printf("%d	%f	%f	%f	%f	%f	%f\n", i, y, jc, t1, t2, t3, t4);
	}
}

// Well, this doesn't do what was vaguely hoped for:
// build shift states that work. Because, well, they're
// not built the right way.
void dump_shift(double omega, double K, double shift)
{
#define SNPTS 500
	for (int i=0; i< SNPTS; i++)
	{
		double y = (i+0.5) / SNPTS;

		double t1 = taverage(y, 20, slope, omega, K);
		double t2 = taverage(y, 40, slope, omega, K);
		double t3 = taverage(y, 80, slope, omega, K);
		double t4 = taverage(y, 160, slope, omega, K);
#if 0
		double t1 = tshift(y, 30, 0.5, slope, omega, K);
		double t2 = trans_half(y, omega, K);
		double t3 = taverage(y, 130, slope, omega, K);
		double t4 = transfer(y, avslope, omega, K);
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
	dump_shift(omega, K, 0.5);
}
