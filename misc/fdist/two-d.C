/*
 * Verify calculations of two-dimensional rep of dyadic group
 * and specifically verify eigenvalues of dyadic GKW
 *
 * March 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double bee(double x)
{
	x -= floor(x);
	if (x<0.5) return 0.0;
	return 1.0;
}

double beew(double w, double x)
{
	int n;
	double wn = 1.0;
	double tn = 1.0;
	int nmax = floor(-25.0 / log(w));
	double acc = 0.0;
	for (n=0; n<nmax; n++)
	{
		double term = wn * bee(tn*x);
		acc +=term;
		wn *= w;
		tn *= 2.0;
	}
	return acc;
}

double beewl(double w, int l, double x)
{
	double tlp1 = 2.0*l+1.0;
	return beew(w, tlp1*x);
}

int main( int argc, char * argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s w l\n", argv[0]);
		exit (1);
	}
	double w = atof(argv[1]);
	int l = atoi(argv[2]);

	printf("#\n# w=%g l=%d\n#\n", w, l);

	int i;
	int imax = 300;
	for (i=0; i<imax; i++)
	{
		double x = ((double) i) / ((double) imax);
		double bwl = beewl(w,l,x);
		printf("%d	%g	%g\n", i, x, bwl);
	}

	return 0;
}
