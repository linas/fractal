
/* 
 * hypergeometric series
 */

#include <complex.h>
#include <math.h>

// double hyperg_2F1(double a, double b, double c, double re_z, double im_z);

complex double 
hyperg_2F1(double a, double b, double c, complex double z)
{
	int n;

	complex double acc = 0.0;
	double an = 1.0;
	double bn = 1.0;
	double cn = 1.0;
	double nf = 1.0;
	complex double zn = 1.0;
	n = 0;
	
	while (1)
	{
		double rat = an / cn;
		rat *= bn /nf;
		complex double term = zn * rat;

		if (n >167) break;
		if (isnan (creal(term)) || isnan (cimag (term))) break;
		acc += term;

		double ta = cabs (term);
		double sa = cabs (acc);
		if (ta < 1.0e-15*sa) break;
		// printf ("duude t= %g +I %g acc=%g +I %g\n", 
		//    creal (term), cimag (term), creal (acc), cimag (acc));

		n++;
		an *= a;
		bn *= b;
		cn *= c;
		nf *= n;
		zn *= z;

		// printf ("n=%d a=%g b=%g c=%g n!=%g\n", n, an, bn, cn, nf);

		a += 1.0;
		b += 1.0;
		c += 1.0;
	}

	return acc;
}


main () 
{
	int i;
	double complex j;

	double a = 1.0/12.0;
	double b = 1.0/12.0;
	double c = 2.0/3.0;

	j = 0.67+0.67I;

	for (i=0; i<20; i++)
	{
		j = 0.9 * cexp (2.0*M_PI*i*0.05);

		double complex f = hyperg_2F1 (a,b,c,j);

		printf ("j=%d t= %g +I %g\n", j, creal (f), cimag (f));
	}

}
