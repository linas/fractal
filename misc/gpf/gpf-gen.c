/*
 * Generating functions for greatest prime factors.
 *
 * April 2016
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gpf.h>
#include "gpf-gen-bignum.h"

/*
 * Ordinary generating function for the greatest common factor.
 */
double gpf_ordinary(double x)
{
	double sum = 0;
	double xn = x;

	if (x < 1.0e-16) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x;
		if (n*xn < 1.0e-16*sum) break;
	}

	return sum;
}

/*
 * Exponential generating function for the greatest common factor.
 */
double gpf_exponential(double x)
{
	double sum = 0;
	double xn = x;

	if (x < 1.0e-16) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x / ((double) n);

		if (n*xn < 1.0e-16*sum) break;
	}

	return sum;
}

double complex gpf_cpx_exponential(double complex z)
{
	double complex sum = 0;
	double complex zn = z;

	if (cabs(z) < 1.0e-16) return z;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * zn;
		zn *= z / ((double) n);

		if (n*cabs(zn) < 1.0e-16*cabs(sum)) break;
	}

	return sum;
}

double complex gpf_cpx_exponential_rt(double r, double theta)
{
	double complex z =  r * (cos(theta) + I *sin(theta));
	return gpf_cpx_exponential(z);
}

double gpf_bignum(double r, double theta)
{
	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);

	cpx_set_d(z, r*cos(theta), r*sin(theta));

	// theta *= 2.0 * M_PI;
	// cpx_gpf_exponential(sum, z, 200);
	cpx_gpf_poch_rising(sum, z, 40);

	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);
	return rv;
}

double complex gpf_cpx_bignum(double r, double theta)
{
	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);

	// theta *= 2.0 * M_PI;
	cpx_set_d(z, r*cos(theta), r*sin(theta));

	// cpx_gpf_exponential(sum, z, 60);
	cpx_gpf_poch_rising(sum, z, 60);

	double complex rv = cpx_get_re(sum) + I * cpx_get_im(sum);
	return rv;
}

int zero_count(double radius)
{
	int count = 0;
	double delta = 0.25 / radius;
	double prev = 0.0;
	for (double theta = 0.0; theta < 2.0*M_PI; theta += delta)
	{
		double complex egz = gpf_cpx_bignum(radius, theta);
		double phase = atan2(cimag(egz), creal(egz));

		if ((0.0 < prev) && (phase < 0.0))
		{
			count ++;
			if ((prev < 2.0) || (-2.0 < phase))
			{
				printf("# step not fine enough at theta=%g\n", theta);
				fprintf(stderr, "Big fail, step not fine enough at theta=%g\n", theta);
				fprintf(stderr, "delta=%g r=%g prev=%g ph=%g\n", delta, radius, prev, phase);
				// exit(1);
			}
		}
		prev = phase;
	}

	return count;
}

int main(int argc, char* argv[])
{
#ifdef ORD
	for (double x=0.0; x< 1.0; x+= 0.002)
	{
		double y = gpf_ordinary(x);
		double z = gpf_exponential(x);
		printf("%g\t%g\t%g\n", x, y, z);
	}
#endif

// #define EXPO
#ifdef EXPO
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <r>\n", argv[0]);
		exit(1);
	}
	double dom = atof(argv[1]);
	printf("#\n# Max = %g\n#\n", dom);
	for (double x=0.0; x< dom; x+= 0.01*dom)
	{
		// double y = gpf_exponential(x);
		// double z = y * exp(-x);
		// printf("%g\t%g\t%g\n", x, y, z);
		double r = x;
		double y = gpf_bignum(r, M_PI);
		// double z = y * log(r) / (r*r);
		// printf("%g\t%20.18g\t%20.18g\t%20.18g\n", x, r, y, z);
		printf("r=%g g=%g\n", r, y);
		fflush(stdout);
	}
#endif
#ifdef RATIONALS
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <r>\n", argv[0]);
		exit(1);
	}
	double dom = atof(argv[1]);
	printf("#\n# Max = %g\n#\n", dom);
	for (double x=0.0; x< dom; x+= 0.001*dom)
	{
		double w0 = gpf_bignum_exponential(x, 0.0);
		double w1_2 = gpf_bignum_exponential(x, 1.0/2.0);
		double w1_3 = gpf_bignum_exponential(x, 1.0/3.0);
		double w1_4 = gpf_bignum_exponential(x, 1.0/4.0);
		double w1_5 = gpf_bignum_exponential(x, 1.0/5.0);
		double w1_6 = gpf_bignum_exponential(x, 1.0/6.0);
		printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\n", x, w0, w1_2, w1_3, w1_4, w1_5, w1_6);
		fflush(stdout);
	}
#endif
#ifdef QUADRATIC_IRRATIONALS
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <r>\n", argv[0]);
		exit(1);
	}
	double dom = atof(argv[1]);
	printf("#\n# Max = %g\n#\n", dom);
	for (double x=0.0; x< dom; x+= 0.0002*dom)
	{
		double complex w0   = gpf_cpx_bignum_exponential(x, sqrt(2.0)/2.0);
		double complex w1_2 = gpf_cpx_bignum_exponential(x, sqrt(2.0)/3.0);
		double complex w1_3 = gpf_cpx_bignum_exponential(x, sqrt(3.0)/2.0);
		double complex w1_4 = gpf_cpx_bignum_exponential(x, sqrt(3.0)/3.0);
		double complex w1_5 = gpf_cpx_bignum_exponential(x, sqrt(5.0)/2.0);
		double complex w1_6 = gpf_cpx_bignum_exponential(x, sqrt(5.0)/3.0);
		printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", x,
			creal(w0), cimag(w0),
			creal(w1_2), cimag(w1_2),
			creal(w1_3), cimag(w1_3),
			creal(w1_4), cimag(w1_4),
			creal(w1_5), cimag(w1_5),
			creal(w1_6), cimag(w1_6));
		fflush(stdout);
	}
#endif

#ifdef TRANSCENDY
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <r>\n", argv[0]);
		exit(1);
	}

	double phi = 0.5 * sqrt(5.0) - 0.5;
	double dom = atof(argv[1]);
	printf("#\n# Max = %g\n#\n", dom);
	for (double x=0.0; x< dom; x+= 0.0002*dom)
	{
		double tp = 2.0 * M_PI;
		double complex w0   = gpf_cpx_bignum_exponential(x, 1.0/tp);
		double complex w1_2 = gpf_cpx_bignum_exponential(x, (1.0/3.0)/tp);
		double complex w1_3 = gpf_cpx_bignum_exponential(x, (2.0/3.0)/tp);
		double complex w1_4 = gpf_cpx_bignum_exponential(x, 0.75/tp);
		double complex w1_5 = gpf_cpx_bignum_exponential(x, phi/tp);
		double complex w1_6 = gpf_cpx_bignum_exponential(x, (5.0/7.0)/tp);
		printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", x,
			creal(w0), cimag(w0),
			creal(w1_2), cimag(w1_2),
			creal(w1_3), cimag(w1_3),
			creal(w1_4), cimag(w1_4),
			creal(w1_5), cimag(w1_5),
			creal(w1_6), cimag(w1_6));
		fflush(stdout);
	}
#endif

#ifdef ZERO_RAYS
	// for (double r=0.87; r< 0.89; r+= 0.001*0.02)
	// for (double r=3.23; r< 3.25; r+= 0.001*0.02)
	// for (double r=3.74; r< 3.75; r+= 0.001*0.02)
	// for (double r=19.95; r< 19.96; r+= 0.001*0.02)
	for (double t=0.83; t< 0.84; t+= 0.001*0.02)
	{
		// double w = gpf_bignum_exponential(r, 0.5 * 0.83295289206477);
		// double w = gpf_bignum_exponential(r, 0.5 * 0.42458721923649);
		// double w = gpf_bignum_exponential(r, 0.5 * 0.59817818048564);
		// double w = gpf_bignum_exponential(r, 0.5 * 0.768238424116);
		// printf("%g\t%g\n", r, w);
		double w = gpf_bignum_exponential(0.88022349562601, 0.5*t);
		printf("%g\t%g\n", t, w);
		fflush(stdout);
	}
#endif
#ifdef FOURIER_ANALYSIS
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <r>\n", argv[0]);
		exit(1);
	}
	double dom = atof(argv[1]);
#define NPTS 2000
	double doms[NPTS+1];
	double vals[NPTS+1];
	int i=0;
	for (double x=10.0; i < NPTS; x+= dom/NPTS)
	{
		doms[i] = x;
		double r = x*x;
		double y = gpf_bignum_exponential(r, 0.0);
		double z = (y * log(r) / (r*r)) - 1.75;
		vals[i] = z;

		// printf("# %g\t%g\n", x, z);
		i++;
		if (i%100==0)
		{
			printf("# Done with %d points to sqrt(r)=%g\n", i, x);
			fflush(stdout);
		}
	}

	// Quick n dirty fourier analysis.
	printf("#\n# fourier %d points to %g\n#\n", NPTS, dom);
	for (double f=0.0; f < 8.0; f+= 0.003)
	{
		double samp = 0.0;
		double camp = 0.0;
		for (i=0; i<NPTS; i++)
		{
			double si = sin(f*doms[i]);
			double co = cos(f*doms[i]);
			samp += si * vals[i];
			camp += co * vals[i];
		}
		samp /= (double) NPTS;
		camp /= (double) NPTS;
		printf("%20.18g\t%20.18g\t%20.18g\n", f, samp, camp);
		fflush(stdout);
	}
#endif

#define ZERO_COUNT
#ifdef ZERO_COUNT
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <rlo> <rhi>\n", argv[0]);
		exit(1);
	}
	double rlo = atof(argv[1]);
	double rhi = atof(argv[2]);
	for (double r=rlo; r<= rhi; r+= 1)
	{
		int count = zero_count(r);
		printf("%g\t%d\n", r, count);
		fflush(stdout);
	}
#endif
}
