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
	return x - 0.5*0.826154195;

	// Bernoulli poly B_2
	// The result is senstive to this being B_2.
	// Being able to integrate to exactly zero is important.
	// return x*x - x  + 1.0 / 6.0;
	// return x*x - x  + 0.16666;

	// Bernoulli poly B_3
	// return x*x*x - 1.5*x*x  + 0.5*x;

	// Bernoulli poly B_4
	// return x*x*x*x - 2.0*x*x*x  + x*x - 1.0/30.0;
}

// ==============================================================
#include "unlambda.c"

int main(int argc, char* argv[])
{
#ifdef COMPARSE_CONST_AND_INVAR
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

	double blam = beta * lambda;

#define NIT 6
	double l1norm[NIT];
	double l2norm[NIT];
	double dotzero[NIT];
	double dotlast[NIT];
	double dotinv[NIT];
	for (int j=0; j<NIT; j++)
	{
		l1norm[j] = 0.0;
		l2norm[j] = 0.0;
		dotzero[j] = 0.0;
		dotlast[j] = 0.0;
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

		double yinv = gp_invar(beta, x);
		printf("%d	%f	%f", i, x, yinv);

		double plm = scan;
		double ylast = yinv;
		for (int j=0; j<NIT; j++)
		{
			double y = nul_n(beta, blam, x, n+j);
			y *= plm;
			plm *= scale;
			printf("	%f", y);

			l1norm[j] += fabs(y) * delta;
			l2norm[j] += y*y * delta;
			dotzero[j] += y * delta;
			dotlast[j] += y*ylast * delta;
			dotinv[j] += y*yinv * delta;
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
	printf("\n# dotzero= ");
	for (int j=0; j<NIT; j++) printf(" %g", dotzero[j]);
	printf("\n# dotlast= ");
	for (int j=0; j<NIT; j++) printf(" %g", dotlast[j]);
	printf("\n# dotinv= ");
	for (int j=0; j<NIT; j++) printf(" %g", dotinv[j]);
	printf("\n#\n");
#endif

}

// ==============================================================
