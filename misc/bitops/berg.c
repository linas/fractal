/* 
 * berg.c
 *
 * Crazy bergman polynomials
 *
 * January 2018
 */

#define NOMAIN
#include "psi.c"
#include "psibig.c"


int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K dim\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	int dim = atoi(argv[2]);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K);

	double lam = 1.0;
	for (int i=0; i<dim; i++)
	{
		double h = hess(K, i+1, i);
		lam /= h;
		printf("%d	%g	%g\n", i, h, lam);
	}
}

 
