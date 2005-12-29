
/*
 * trinary.c
 *
 * try to come up with a 3-adic type continued fraction
 */

#include <complex.h>
#include <math.h>

main ()
{
	int i;
	double complex gen, elt;

	gen = cexp (2.0*M_PI*I/3.0);

	printf ("generator is %g +i%g \n", creal (gen), cimag(gen));

	double x = 0.125425;

	double complex vt = x;
	double complex mult;

	elt = gen;
	for (i=0; i<26; i++)
	{
		vt = cpow (vt, gen);
		mult = vt / elt;
		double m = creal (mult);
		m = cabs(mult);
		m = floor (m);
		// if (0.0 > m) m+=1.0;
		// m = trunc (m);
		int mm = m;
		printf (" its %d %g +i%g  and %g +i%g\n",  mm,
							 creal(vt), cimag(vt),
							 creal(mult), cimag(mult));
		vt -= m *elt;
		
		elt *= gen;
	}
}
