/*
 * force.c
 * Blindly hunt for eigenfunctions of the transfer operator.
 * This is done by iterating, and looking for fixed points that
 * might accidentally appear in the ocean. A handful are found
 * and pinned down.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double invar(double beta, double x)
{
	double midpnt = 0.5*beta;
	double obn = 1.0;
	double sum = 0.0;
	double norm = 0.0;
	for (int i=0; i<1000; i++)
	{
		if (x < midpnt) sum += obn;
		norm += midpnt*obn;

		if (0.5 < midpnt) midpnt -= 0.5;
		midpnt *= beta;
		obn /= beta;
		if (obn < 1e-15) break;
	}
	return sum / norm;
}

/* Some magic numbers!
 *
 * #define NHIST 27717  with "SAW" intialization and with beta=1.61803
 * converges in a stable way to the first decaying golden eigenfunc.
 * Other intializations are NOT stable and do NOT converge!
 * Other numbers of bins do NOT result in a stable convergence!
 * Even the explicit initialization of the eigenfunc is NOT  stable!
 * Other beta values are NOT stable, including beta=1.618034
 *
 * What is causing (in)stability?
 *
 * Other magic numbers:
 * Using the SAW or GOLD_ONE initialization:
 * #define NHIST 28658 (fibo 22) + 1 with beta=1.618034 has ringing pattern
 * #define NHIST any fib +1 seems to work, with ringing.
 */

// #define NHIST 360361
// #define NHIST 160361
// #define NHIST 27721
// #define NHIST 27719
// #define NHIST 27717
// #define NHIST 28657 // (fibo 22)
#define NHIST 28658
// #define NHIST 196418 // (fib 26)
// #define NHIST 196419

// #define NHIST 18561 // n=2 fib plus one
// #define NHIST 25282    // n=3 fib plus one

double histo[NHIST];
double histn[NHIST];

int enx (double x)
{
	int n = floor (NHIST * x);
	if (0 > n) n=0;
	if (NHIST <= n) n=NHIST-1;
	return n;
}

double dense(double x)
{
	if (1.0 <= x) return 0.0;
	return histo[enx(x)];
}

double ell(double beta, double y)
{
	if (0.5 * beta < y) return 0.0;

	double xlo = y / beta;
	double xhi = xlo + 0.5;

	double dlo = dense(xlo);
	double dhi = dense(xhi);

	double ellie = dlo + dhi;
	ellie /= beta;
	return ellie;
}

// Kernel generator, used to build the kernel.
// Can be anything at all.
double kerg(double x, double beta)
{
	double a = 0.25*(beta-1.0);
	double m = 11.481425;  // Overall scale factor, doesn't matter.
	// return m*(x-a);
	return m*(x-a)*(x-a);
}

// This function is in the kernel of the xfer function, for
// any possible kerg() function above. If and only if; there
// are no others.
double saw(double x, double beta)
{
	x -= floor(x);  // mod 1
	if (0.5*beta < x) return 0.0;
	if (x < 0.5*(beta-1.0)) return kerg(x, beta);
	if (0.5 < x) return -kerg(x-0.5, beta);
	return 0.0;
}

double blancmange(double x, int l, double w, double beta)
{
	if (0.5*beta < x) return 0.0;

	double xn = x;
	double wn = 1.0;
	double sum = 0.0;
	double tlp = 2*l+1;
	for (int i=0; i<1000; i++)
	// for (int i=0; i<1; i++)
	{
		sum += wn * saw(tlp * xn, beta);
		wn *= w;
		if (0.5 < xn) xn -= 0.5;
		xn *= beta;
		if (wn < 1e-15) break;
	}

	// Wow. This almost works ...
	// sum *= invar(beta,x);
	// sum *= 2.0 + invar(beta,x);

	return sum;
}

void blanc_setup(double beta, int l, double w)
{
	for (int i=0; i<NHIST; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = blancmange(x, l, w, beta);
	}
}

void setup(double beta)
{
// #define SAW
#ifdef SAW
	for (int i=0; i<NHIST; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = x-0.5;
	}
#endif

// #define STEP
#ifdef STEP
	int midp = enx(0.25*beta);
	for (int i=0; i<midp; i++) histn[i] = 1.0;
	for (int i=midp; i<2*midp; i++) histn[i] = -1.0;
	for (int i=2*midp; i<NHIST; i++) histn[i] = 0.0;
#endif

// #define GOLD_ONE
#ifdef GOLD_ONE
	double sum2 = 0.0;
	int half = enx(0.5);
	for (int i=0; i<half; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = beta*x-0.5;
		sum2 += histn[i];
	}

	double sum1 = 0.0;
	int endp = enx(0.5*beta);
	for (int i=half; i<endp; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = x-0.5;
		sum1 += histn[i];
	}

	for (int i=endp; i<NHIST; i++) histn[i] = 0.0;

	printf("# gold one sum_2 = %g sum_1 = %g\n", sum2, sum1);
	sum2 /= (double) half;
	sum1 /= (double) endp-half;
	printf("# Interval-weighted sum_2 = %g sum_1 = %g ratio s1/s2= %g\n",
		sum2, sum1, sum1/sum2);
#endif

// #define QUAD_ONE
#ifdef QUAD_ONE
	// Works great, now that we found all the algebra mistakes.
	double lambda = 1.0 / (beta*beta);
	double a = -1.0;
	double b = 0.125 * beta;
	printf("# quadratic lambda = %g a=%g b=%g\n", lambda, a, b);

	double sum2 = 0.0;
	int half = enx(0.5);
	for (int i=0; i<half; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = (beta*beta*x*x +a*beta*x + b)*beta*lambda;
		sum2 += histn[i];
	}
	double sum1 = 0.0;
	int endp = enx(0.5*beta);
	for (int i=half; i<endp; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = x*x +a*x + b;
		sum1 += histn[i];
	}
	for (int i=endp; i<NHIST; i++) histn[i] = 0.0;

	printf("# qudratic sum_2 = %g sum_1 = %g\n", sum2, sum1);
	sum2 /= (double) half;
	sum1 /= (double) endp-half;
	printf("# Interval-weighted sum_2 = %g sum_1 = %g ratio s1/s2= %g\n",
		sum2, sum1, sum1/sum2);
#endif

// #define GOLD_TWO
#ifdef GOLD_TWO
	int m1 = enx(0.5*beta*(beta-1.0));
	for (int i=0; i<m1; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = beta*beta*x-0.5;
	}
	int half = NHIST/2;
	for (int i=m1; i<half; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = beta*x-0.5;
	}
	int endp = enx(0.5*beta);
	for (int i=half; i<endp; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = x-0.5;
	}
	for (int i=endp; i<NHIST; i++) histn[i] = 0.0;
#endif

// #define GOLD_THREE
#ifdef GOLD_THREE
	// Bingo! This works for index=3

	double em1 = 0.5*beta*(beta-1.0);
	// Either defintiion of C works.
	// double C = 0.5 + 0.25 / beta;
	double C = em1 - 0.25 / beta;

	// ival 3 -- correct by inference
	int half = enx(0.5);
	for (int i=0; i<half; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = beta*x - C;
	}

	// ival 2
	int m1 = enx(em1);
	for (int i=half; i<m1; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		// Either version works
		// histn[i] = ((beta + 1.0)/ beta) *x +0.5 - 2.0*C;
		histn[i] = ((beta + 1.0)/ beta) *x - em1;
	}

	// ival 1  --
	int endp = enx(0.5*beta);
	for (int i=m1; i<endp; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = x - C;
	}
	for (int i=endp; i<NHIST; i++) histn[i] = 0.0;
#endif

// #define GOLD_FOUR
#ifdef GOLD_FOUR
	// Works great!
	int m1 = enx(0.5*beta*(beta-1.0));
	for (int i=0; i<m1; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = beta*beta*beta*x-0.5;
	}
	int nit = enx(0.5/beta);
	for (int i=m1; i<nit; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = beta*beta*x-0.5;
	}
	int half = enx(0.5);
	for (int i=nit; i<half; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = beta*x-0.5;
	}
	int endp = enx(0.5*beta);
	for (int i=half; i<endp; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NHIST);
		histn[i] = x-0.5;
	}
	for (int i=endp; i<NHIST; i++) histn[i] = 0.0;
#endif
}

double normalize(double beta)
{
	// Remove constant
	int m0 = enx(0.5*beta);
	double delta = 1.0 / ((double) m0);

	double sum = 0.0;
	for (int i=0; i<m0; i++) sum += histn[i];
	for (int i=0; i<m0; i++) histn[i] -= sum * delta;

	// Normalize
	double norm = 0.0;
	for (int i=0; i<m0; i++) norm += fabs(histn[i]);
	norm *= delta;

	// Hack to avoid divide by zero.
	if (norm < 0.001) norm = 1.0;
	for (int i=0; i<m0; i++) histo[i] = histn[i] / norm;

	printf("# renorm, sum=%g norm = %g\n", sum * delta, norm);
	fprintf(stderr, "Renorm -> sum=%g norm = %g\n", sum * delta, norm);
	return norm;
}

double step(double beta)
{
	double delta =1.0 / ((double) NHIST);

	// Run step
	for (int i=0; i<NHIST; i++)
	{
		double x = (((double) i) + 0.5) * delta;
		histn[i] = ell(beta, x);
	}

	return normalize(beta);
}

#define NCAP 10
double capt[NCAP][NHIST];

void capture(int n)
{
	for (int i=0; i<NHIST; i++) capt[n][i] = histo[i];
}

int main(int argc, char* argv[])
{
	if (argc != 5)
	{
		// fprintf(stderr, "Usage: %s beta nsteps\n", argv[0]);
		fprintf(stderr, "Usage: %s beta nsteps ll w\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int nsteps = atoi(argv[2]);
	int ll = atoi(argv[3]);
	double w = atof(argv[4]);

	// setup(beta);
	blanc_setup(beta, ll, w);
	printf("#\n# beta=%g NHIST=%d\n", beta, NHIST);
	printf("# blanc w=%g l=%d eig=%g\n", w, ll, 2*w/beta);
	fprintf(stderr, "Beta=%g blanc w=%g l=%d eig=%g\n", beta, w, ll, 2*w/beta);

	normalize(beta);

// #define SHOW_NORMS
#ifdef SHOW_NORMS
	for (int i=0; i< nsteps; i++)
	{
		double lam = step(beta);
		printf("%d	%g\n", i, lam);
	}
#endif

#define SHOW_DENS
#ifdef SHOW_DENS
	for (int i=0; i< nsteps; i++)
		step(beta);

	// double clam[NCAP];

	fprintf(stderr, "----- nstep=%d\n", nsteps);
	for (int i=0; i< NCAP; i++)
	{
		capture(i);
		// clam[i] =
		step(beta);
	}

	for (int j=0; j< NHIST; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NHIST);
		printf("%d	%g", j, x);
		// printf("	%g", invar(beta, x));

		for (int i=0; i< NCAP-1; i++)
		{
			printf("	%g", capt[i][j]);
			// printf("	%g", capt[i+1][j] - capt[i][j]);
			//if (1.0e-2 < fabs(capt[i][j]))
			// printf("	%g", (capt[i+1][j] - capt[i][j]) / capt[i][j]);
		}
		printf("\n");
	}
#endif

// #define CAPTURE_DENS
#ifdef CAPTURE_DENS
	// Record only the almost-eigenfunctions
	int ncap = 0;
	for (int i=0; i< nsteps; i++)
	{
		double lam = step(beta);
		if (fabs(lam - 1.0/beta) < 0.004)
		{
			capture(ncap);
			ncap++;
			if (NCAP <= ncap) break;
		}
	}

	printf("#\n# beta=%g ncaptures = %d\n#\n", beta, ncap);

	for (int j=0; j< NHIST; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NHIST);
		printf("%d	%g", j, x);

		for (int i=0; i< ncap; i++)
		{
			printf("	%g", capt[i][j]);
		}
		printf("\n");
	}
#endif
}
