
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
	
	while (1)
	{
		complex double term = an*bn*zn / (cn*nf);
		acc += term;

		double ta = cabs (term);
		double sa = cabs (acc);
		if (ta < 1.0e-15*sa) break;

		an *= a;
		bn *= b;
		cn *= c;
		nf *= n;
		zn *= z;

		a += 1.0;
		b += 1.0;
		c += 1.0;
		n++;
	}

	return acc;
}


main () 
{

	double complex j;

	double a = 1.0/12.0;
	double b = 1.0/12.0;
	double c = 2.0/3.0;

	j = 0.1;

	double complex f = hyperg_2F1 (a,b,c,j);

	printf (" duude result = %g +I %g\n", creal (f), cimag (f));

}
