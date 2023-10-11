
/*
 * mathcheck.c
 *
 * do some simple floating-point validation
 * Linas Dec 2003
 */

#include <math.h>
#include <complex.h>

main ()
{
	long double x;
	x = 1.012345678901234567890123456789012345L;

	printf ("long double = %40.35Lg\n", x);

	double complex c = 1.012345678901234567890123456789012345;
	c += 3.012345678901234567890123456789012345I;

	printf ("complex = %25.20g +i %25.20g\n", creal(c), cimag(c));
	
	long double complex cl = 1.012345678901234567890123456789012345L;
	cl += 3.012345678901234567890123456789012345LI;

	printf ("long complex = %40.35Lg +i %40.35Lg\n", creall(cl), cimagl(cl));
	
	x = 0.0123456789L;
	x *= 1.0e-10L;
	x += 0.0123456789L;
	x *= 1.0e-10L;
	x += 0.0123456789L;
	x *= 1.0e-10L;
	x += 0.0123456789L;

	printf ("long double = %40.35Lg\n", x);
	x -= 0.0123456789L;
	x *= 1.0e10L;

	printf ("long double = %40.35Lg\n", x);
	x -= 0.0123456789L;
	x *= 1.0e10L;
			  
	printf ("long double = %40.35Lg\n", x);
	x -= 0.0123456789L;
	x *= 1.0e10L;

	printf ("long double = %40.35Lg\n", x);
	x -= 0.0123456789L;
	x *= 1.0e10L;
	printf ("long double = %40.35Lg\n", x);
	
}
