/*
 * force.c
 * Recheck old (negative) results.
 * They're still negative.
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
	return sum;
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
// #define NHIST 28658
// #define NHIST 196418 // (fib 26)
// #define NHIST 196419

// #define NHIST 18561 // n=2 fib plus one
#define NHIST 25282    // n=3 fib plus one

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

#define QUAD_ONE
#ifdef QUAD_ONE
	// Broken
	// Should have worked, but I guess there's an algebra error somewhere
	double lambda = 0.5 * (beta + sqrt(beta + 5.0)) / (2.0*beta + 1.0);
	double a = -beta;
	double b = 0.25 * (1.0 - 2.0*beta) / (lambda*lambda*beta*beta - lambda*beta -1.0);
	printf("# quadratic lambda = %g a=%g b=%g\n", lambda, a, b);

double b2=beta*beta;
double b4=b2*b2;
double b6=b4*b2;
double b8=b4*b4;
double lalt = (b4+sqrt(b8+4*b6)) / (2*b6);
printf("lalt=%g\n", lalt);

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

printf("f1 at m0 should be %g\n", 0.25*beta*beta + 0.5*beta*a +b);
printf("f1 at m0 founde %g\n", histn[endp-1]);

printf("f2 at m1 should be %g\n", lambda*beta*(0.25*beta*beta + 0.5*beta*a +b));
printf("f2 at m1 found %g\n", histn[half-1]);

double enf2 = lambda*beta*(0.25*beta*beta + 0.5*beta*a +b);
enf2 *= lambda*beta;
printf("eiegf2 at m1 %g\n", enf2);
double xlo=0.5/beta;
double lef = lambda*beta*(beta*beta*xlo*xlo + beta*a*xlo +b);
printf("lef at m1 %g vs array %g\n", lef, histn[enx(xlo)]);
double xhi=xlo + 0.5;
double rig = (xhi*xhi + a*xhi +b);
printf("rig at m1 %g vs aray %g\n", rig, histn[enx(xhi)-1]);
printf("sum lef+rig at m1 %g\n", lef+rig);
printf("-------\n");

printf("quad terms eig %g vs %g + %g = %g \n", lambda*lambda*b4, lambda*b2, 1.0/b2,
lambda*b2 + 1.0/b2);

printf("line terms eig %g vs %g + %g = %g \n", 
lambda*lambda*b2*beta*a, lambda*beta*a, 1.0/beta + a/beta,
lambda*beta*a + 1.0/beta + a/beta);

printf("cone =%g  vs %g + %g = %g\n", lambda*lambda*b2*b,
lambda*beta*b, 0.25 + 0.5*a+b,
lambda*beta*b + 0.25 + 0.5*a+b
);

printf("-------\n");
double acc = 0.0;
for (int i=half; i<endp; i++)
{
	double x = (((double) i) + 0.5) / ((double) NHIST);
	int j = enx(x);
	double f1 = histn[j];
	f1 *= lambda*beta;
	int k = enx(x/beta);
	double f2 = histn[k];
	double ad = fabs(f1-f2);
	acc += ad;
	if (5e-5 < ad)
		printf("%d %g expect %g got %g diff %g\n", i, x, f1, f2, f1-f2);
}
printf("ovr %d total acc=%g\n", endp-half, acc);

#if 1
acc = 0.0;
for (int i=0; i<half; i++)
{
	double x = (((double) i) + 0.5) / ((double) NHIST);
	double ef2 = lambda*beta*histn[enx(x)];
	double lef = histn[enx(x/beta)];
	double rig = histn[enx(0.5+x/beta)];
	double ad = fabs(ef2 - (lef+rig));
	acc += ad;
	// if (5e-5 < ad)
		// printf("%d %g egn %g got %g + %g diff %g ska %g\n", i, x, ef2, lef, rig, ad, ad/x);

   // double egn = lambda*lambda*b2*( b2*x*x +a*beta*x + b);
	// printf("at %d %g expect egn=%g got %g delt=%g\n", i, x, egn, ef2, egn-ef2);

	double ylo = lambda*beta*( beta*x*x +a*x + b);
	printf("at %d %g expect ylo=%g got %g delt=%g\n", i, x, ylo, lef, ylo-lef);

	// double xhi = x/beta+0.5;
	// double yhi = xhi*xhi + a*x/beta + 0.5*a + b;
	// printf("at %d %g expect yhi=%g got %g delt=%g\n", i, x, yhi, rig, yhi-rig);
}
printf("yepovr %d total acc=%g\n", half, acc);
#endif

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
	for (int i=0; i<m0; i++) histo[i] = histn[i] / norm;

	printf("# renorm, sum=%g norm = %g\n", sum * delta, norm);
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
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s beta nsteps\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int nsteps = atoi(argv[2]);

	setup(beta);
	normalize(beta);

// #define SHOW_NORMS
#ifdef SHOW_NORMS
	printf("#\n# beta=%g NHIST=%d\n#\n", beta, NHIST);
	for (int i=0; i< nsteps; i++)
	{
		double lam = step(beta);
		printf("%d	%g\n", i, lam);
	}
#endif

#define SHOW_DENS
#ifdef SHOW_DENS
	printf("#\n# beta=%g\n#\n", beta);

	for (int i=0; i< nsteps; i++)
		step(beta);

	for (int i=0; i< NCAP; i++)
	{
		capture(i);
		step(beta);
	}

	for (int j=0; j< NHIST; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NHIST);
		printf("%d	%g", j, x);

		for (int i=0; i< NCAP; i++)
		{
			printf("	%g", capt[i][j]);
		}
		printf("\n");
	}
#endif

// #define CAPTURE_DENS
#ifdef CAPTURE_DENS
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
