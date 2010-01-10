
/*
 * graph of hurwitz zeta
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_zeta.h>

// Egeinvector:
/**
#
# right 1'th eigenvector[0]=0.570509 (normalized=1)
# right 1'th eigenvector[1]=-0.416245 (normalized=-0.729603)
# right 1'th eigenvector[2]=-0.511962 (normalized=-0.897377)
# right 1'th eigenvector[3]=-0.38492 (normalized=-0.674695)
# right 1'th eigenvector[4]=-0.24617 (normalized=-0.431492)
# right 1'th eigenvector[5]=-0.145082 (normalized=-0.254304)
# right 1'th eigenvector[6]=-0.0814122 (normalized=-0.142701)
# right 1'th eigenvector[7]=-0.0442329 (normalized=-0.0775324)
# right 1'th eigenvector[8]=-0.0235004 (normalized=-0.0411919)
# right 1'th eigenvector[9]=-0.0122866 (normalized=-0.0215362)
# right 1'th eigenvector[10]=-0.00634884 (normalized=-0.0111284)
# right 1'th eigenvector[11]=-0.00325233 (normalized=-0.00570075)
# right 1'th eigenvector[12]=-0.00165539 (normalized=-0.0029016)
# right 1'th eigenvector[13]=-0.000838545 (normalized=-0.00146982)
# right 1'th eigenvector[14]=-0.000423261 (normalized=-0.0007419)
# right 1'th eigenvector[15]=-0.000213082 (normalized=-0.000373494)
# right 1'th eigenvector[16]=-0.000107064 (normalized=-0.000187664)
# right 1'th eigenvector[17]=-5.37183e-05 (normalized=-9.41585e-05)
# right 1'th eigenvector[18]=-2.69248e-05 (normalized=-4.71944e-05)
# right 1'th eigenvector[19]=-1.34853e-05 (normalized=-2.36373e-05)
# right 1'th eigenvector[20]=-6.75051e-06 (normalized=-1.18324e-05)
**/

/** first non-trivial Right-hand side eigenvector */
double revec(double x)
{
	int n;
	double y, yn, sum, norm;
	double eigenvector[21];
	eigenvector[0]=0.570509;
	eigenvector[1]=-0.416245;
	eigenvector[2]=-0.511962;
	eigenvector[3]=-0.38492;
	eigenvector[4]=-0.24617;
	eigenvector[5]=-0.145082;
	eigenvector[6]=-0.0814122;
	eigenvector[7]=-0.0442329;
	eigenvector[8]=-0.0235004;
	eigenvector[9]=-0.0122866;
	eigenvector[10]=-0.00634884;
	eigenvector[11]=-0.00325233;
	eigenvector[12]=-0.00165539;
	eigenvector[13]=-0.000838545;
	eigenvector[14]=-0.000423261;
	eigenvector[15]=-0.000213082;
	eigenvector[16]=-0.000107064;
	eigenvector[17]=-5.37183e-05;
	eigenvector[18]=-2.69248e-05;
	eigenvector[19]=-1.34853e-05;
	eigenvector[20]=-6.75051e-06;

	y = 1.0 - x;

	sum = 0.0;
	norm = 0.0;
	yn = 1.0;
	for (n=0; n<21; n++)
	{
		sum += eigenvector[n] * yn;
		norm += eigenvector[n];
		yn *= y;
	}

	return sum/norm;
}



main (int argc, char * argv[])
{
	int i;
	double s = 2.0;

	s = atof (argv[1]);

	int npts = 300;
	double n2 = gsl_sf_hzeta(2.0, 1.0);
	double n25 = gsl_sf_hzeta(2.5, 1.0);
	double n3 = gsl_sf_hzeta(3.0, 1.0);
	double n4 = gsl_sf_hzeta(4.0, 1.0);

	// double norm = 1.0 / (n2 - s * n3);
	// double norm = 1.0 / (n2 - s * n4);
	double norm = 1.0 / (n2 - s);
	// double norm = 1.0 / (n25 - s);
	for (i=0; i< npts; i++)
	{
		double x = ((double) i) / ((double) npts);
		double y2 = gsl_sf_hzeta(2.0, 1.0+x);
		// double y25 = gsl_sf_hzeta(2.5, 1.0+x);
		// double y3 = gsl_sf_hzeta(3.0, 1.0+x);
		// double y4 = gsl_sf_hzeta(4.0, 1.0+x);

		// This one almost works for s=1
		// double y = norm * (y2 - s * y3);

		// almost works for s=2.5
		// double y = norm * (y2 - s * y4);

		// fair for s=1
		// double y = norm * (y2 - s);

		// fair for s=1.44
		// double y = norm * (y2 - s/(1.0+x));

		// quite good for s=1.2
		// but a little too flat
		// double y = norm * (y2 - s/sqrt(1.0+x));

		// really very good for s=1.28
		// just a shade too flat, though
		// double y = norm * (y2 - s/pow(1.0+x, 0.666666666));

		// Wow, really very good, for s=1.32
		// maybe a shade too curved .. hard to say ...
		double y = norm * (y2 - s/pow(1.0+x, 0.75));

		// quite good for s=0.66
		// but just a little too flat.
		// double y = norm * (y25 - s);

      // really quite good for s=0.8
		// but just a little too curved ...
		// double y = norm * (y25 - s/sqrt(1.0+x));

		double r = revec(x);
		printf("%d	%g	%g	%g\n", i, x, y, r);
	}
}
