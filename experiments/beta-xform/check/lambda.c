/*
 * lambda.c
 *
 * Version of unstack.c with lambda inserted into it.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "unref.c"

// ==============================================================

int mode = 0;

// #define QUADRATIC_GIVES_LINEAR
#ifdef QUADRATIC_GIVES_LINEAR
// This coefficient combo used on the quadratic seems to give
// another ergodic function sequence when iterated. But the
// characteristic eigennorm is 1/beta and the functions are
// linear-like; the parabola is washed out.
double c1 = 1.0;
double c0 = 0.9008;
#endif

// This works, but the linear slopes up! The other sloped down!
// double c1 = 0.0;
// double c0 = -1.578;

// This one is flat flat but eigen is still 1/beta. Its not going
// to be parabolic.
double c1 = 0.83;
double c0 = 0.4795;

// Arbitrary function
double nu(double x)
{
	if (x < 0.0) fprintf(stderr, "Error nu fail neg %g\n", x);
	if (1.0 < x) fprintf(stderr, "Error nu fail pos %g\n", x);

#if MODAL
	if (0 == mode) return 1.0;
	if (1 == mode) return x;
	return gp_invar(1.6, x);
#endif

	// return 1.0;
	// return x-0.5;
	// return x - 0.5*0.8262;     // appropriate for beta=1.6
	// return x - 0.5*0.82616;
	// return x - 0.5*0.826154;
	// return x - 0.5*0.8261542;
	// return x - 0.5*0.82615419;
	// return x - 0.5*0.826154195;

	// Bernoulli poly B_2
	// The result is senstive to this being B_2.
	// Being able to integrate to exactly zero is important.
	// return x*x - x  + 1.0 / 6.0;
	// return x*x - x  + 0.16666;
	return x*x - c1*x  + c0 / 6.0;

	// Bernoulli poly B_3
	// return x*x*x - 1.5*x*x  + 0.5*x;

	// Bernoulli poly B_4
	// return x*x*x*x - 2.0*x*x*x  + x*x - 1.0/30.0;
}

// ==============================================================
#include "unlambda.c"

int main(int argc, char* argv[])
{
#ifdef COMPARE_CONST_AND_INVAR
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s beta n\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int n = atoi(argv[2]);

	int imax = 814;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);

		double y = gp_invar(beta, x);
		printf("%d	%f	%f", i, x, y);

		double blam = 1.0;
		mode = 0;
		double z1 = nul_n(beta, blam, x, n);
		mode = 22;
		double z2 = nul_n(beta, blam, x, n);
		printf("	%f	%f	%f\n", z1, z2, z2/z1);
		fflush(stdout);
	}
#endif

#define DO_GRAPH
#ifdef DO_GRAPH
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s beta lambda n\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	double lambda = atof(argv[2]);
	int n = atoi(argv[3]);

	// c1 = atof(argv[4]);
	// c0 = atof(argv[5]);

	double blam = beta * lambda;

#define NIT 6
	double l1norm[NIT];
	double l2norm[NIT];
	double l2orth[NIT];
	double dotzero[NIT];
	double dotlast[NIT];
	double dotorth[NIT];
	double dotinv[NIT];
	for (int j=0; j<NIT; j++)
	{
		l1norm[j] = 0.0;
		l2norm[j] = 0.0;
		l2orth[j] = 0.0;
		dotzero[j] = 0.0;
		dotlast[j] = 0.0;
		dotorth[j] = 0.0;
		dotinv[j] = 0.0;
	}

	double scale = lambda * beta;
	// scale = lambda;
	double scan = pow(scale, n);

	printf("#\n# beta=%f n=%d scale=%f\n#\n", beta, n, scale);
	fflush(stdout);

	int imax = 814;
	double delta = 1.0 / ((double) imax);
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);

		double yinvrnt = gp_invar(beta, x);
		printf("%d	%f	%f", i, x, yinvrnt);

		double plm = scan;
		double ylast = yinvrnt;
		for (int j=0; j<NIT; j++)
		{
			double y = nul_n(beta, blam, x, n+j);
			y *= plm;
			plm *= scale;
			printf("	%f", y);

			l1norm[j] += fabs(y) * delta;
			l2norm[j] += y*y * delta;
			l2orth[j] += (y*y/yinvrnt) * delta;
			dotzero[j] += y * delta;
			dotlast[j] += y*ylast * delta;
			dotorth[j] += (y*ylast/yinvrnt) * delta;
			dotinv[j] += y*yinvrnt * delta;
			ylast = y;
		}
		printf("\n");
		fflush(stdout);
	}

	printf("#\n");
	printf("# l1norm= ");
	for (int j=0; j<NIT; j++) printf(" %g", l1norm[j]);
	printf("\n# l2norm= ");
	for (int j=0; j<NIT; j++) printf(" %g", l2norm[j]);
	printf("\n# l2orth= ");
	for (int j=0; j<NIT; j++) printf(" %g", l2orth[j]);
	printf("\n# dotzero= ");
	for (int j=0; j<NIT; j++) printf(" %g", dotzero[j]);
	printf("\n# dotlast= ");
	for (int j=0; j<NIT; j++) printf(" %g", dotlast[j]);
	printf("\n# dotorth= ");
	for (int j=0; j<NIT; j++) printf(" %g", dotorth[j]);
	printf("\n# dotinv= ");
	for (int j=0; j<NIT; j++) printf(" %g", dotinv[j]);
	printf("\n#\n");

	printf("# cos= ");
	for (int j=0; j<NIT; j++)
	{
		double chk = beta * dotorth[j] / l2orth[j];
		printf(" %g", chk);
	}
	printf("\n#\n");
#endif

// #define BISECT
#ifdef BISECT
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s beta n\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int n = atoi(argv[2]);
	c0 = atof(argv[3]);

	double blam = beta;
	double x = 0.01;
	double yinvrnt = gp_invar(beta, x);
	double scale = beta;
	double scan = pow(scale, n);
	double y = nul_n(beta, blam, x, n);
	y *= scan;
	y /= yinvrnt;

	printf("%g n=%d rat= %g\n", c0, n, y);
#endif
}

// ==============================================================
