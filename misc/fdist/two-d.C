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
	int nmax = floor(-25.0 / log(fabs(w)));
	double acc = 0.0;
	for (n=0; n<nmax; n++)
	{
		double term = wn * bee(tn*x);
		acc +=term;
		wn *= w;
		tn *= 2.0;
	}
	return (1.0-w)*acc;
}

double beewl(double w, int l, double x)
{
	double tlp1 = 2.0*l+1.0;
	return beew(w, tlp1*x);
}

double Lcbee(double w, int l, double x)
{
	int n;
	double acc = 0.0;
	double tn = 0.5;
	for (n=1; n<40; n++)
	{
		double term = beewl(w,l, (2.0-x)*tn);
		acc += tn * term;
		tn *= 0.5;
	}
	return acc;
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
	int imax = 2400;
	for (i=0; i<imax; i++)
	{
		double x = ((double) i) / ((double) imax);
		double bwl = beewl(w,l,x);
		// bwl *= 1.0-w;

		double lcbwl = Lcbee(w,l,x);

		double gbwl = 0.0;
		switch(l)
		{
		case 0:
			// the l=0 case -- works
			gbwl = (1.0 -w * bwl) / (2.0-w);
			break;

		case 1:
			// the l=1 case
			if (x < 1.0/3.0)
				gbwl = (0.25*(5.0-w) -w * bwl) / (2.0-w);
			else if (x < 2.0/3.0)
				gbwl = (0.25*(1.0+3.0*w) -w * bwl) / (2.0-w);
			else
				gbwl = (0.5*(3.0-w) -w * bwl) / (2.0-w);
			break;

		case 2:
			// the  l=2 case 
			if (x < 1.0/5.0)
				gbwl = (0.125*(9.0-w) -w * bwl) / (2.0-w);
			else if (x < 2.0/5.0)
				gbwl = (0.125*(1.0+7.0*w) -w * bwl) / (2.0-w);
			else if (x < 3.0/5.0)
				gbwl = (0.25*(7.0-3.0*w) -w * bwl) / (2.0-w);
			else if (x < 4.0/5.0)
				gbwl = (0.25*(3.0+w) -w * bwl) / (2.0-w);
			else 
				gbwl = (0.25*(5.0-w) -w * bwl) / (2.0-w);
			break;

		case 3:
			// the l=3 case -- for first 1/7th
			gbwl = (0.125*(11.0-3.0*w) -w * bwl) / (2.0-w);
			break;

		case 4:
			// the l=4 case -- for first 1/9th
			gbwl = (0.0625*(17.0-w) -w * bwl) / (2.0-w);
			break;

		case 5:
			// the l=5 case -- for first 1/11th
			gbwl = (0.0625*(21.0-5.0*w) -w * bwl) / (2.0-w);
			break;

		case 6:
			// the l=6 case -- for first 1/13th
			gbwl = (0.0625*(19.0-3.0*w) -w * bwl) / (2.0-w);
			break;

		case 7:
			// the l=7 case -- for first 1/15th
			gbwl = (0.0625*(23.0-7.0*w) -w * bwl) / (2.0-w);
			break;

		case 8:
			// the l=8 case -- for first 1/15th
			gbwl = (0.03125*(33.0-w) -w * bwl) / (2.0-w);
			break;
		}

		double diff = lcbwl - gbwl;

		printf("%d	%8.5g	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, bwl, lcbwl, gbwl, diff);
	}

	return 0;
}
