

/* quickie test case 


*/

#include <math.h>
#include <values.h>
#include <stdio.h>
#include "bessel.h"

#define EPS 5.0e-13
#define J0_EPS EPS
int
main ()
{
	int i;
	double z, b[600];
	double s, c, u, err;

#define ZRANGE 1000.0
#define NRANGE 10

	printf ("Max double value=%g\n", MAXDOUBLE);

	/* perform broad tests */
	for (z=-ZRANGE; z<ZRANGE; z+=1.0)
	{
		bessel (NRANGE, z, b);

		s = sin(z)/z;
		err = fabs ((s-b[0])/s);
		if (J0_EPS < err)
		{
			printf ("\nBad sinc Test: z=%g me=%g sys=%g err=%g\n", 
				z, b[0], s, err);
		}

		/* compare to j1 */
		s = (s-cos(z))/z;
		err = fabs ((s-b[1])/s);
		if (EPS < err)
		{
			printf ("\nBad j1 Test: z=%g me=%g sys=%g err=%g\n", 
				z, b[1], s, err);
		}

		/* compare to sine and cosine */
		s = z * b[0];
		c = b[0] - z * b[1];

		u = s*s+c*c;
		err = u - 1.0;
		if (J0_EPS < err)
		{
			printf ("\nBad sin/cos Test: z=%g err=%g\n", 
				z, err);
		}
		


#ifdef HAVE_SPH
These etest are broken because j0, j1, jn are cylindrical bessels, not
sphereical bessels.

		s = j0 (z);
		err = fabs ((s-b[0])/s);
		if (EPS < err)
		{
			printf ("Bad j0 Test: z=%g me=%g sys=%g err=%g\n", 
				z, b[0], s, err);
		}

		s = j1 (z);
		err = fabs ((s-b[1])/s);
		if (EPS < err)
		{
			printf ("Bad j1 Test: z=%g me=%g sys=%g err=%g\n", 
				z, b[1], s, err);
		}

		for (i=0; i<NRANGE; i++) 
		{
			s = jn (i,z);
			err = fabs ((s-b[i])/s);
			if (EPS < err)
			{
				printf ("Bad Test: i=%d z=%g me=%g sys=%g err=%g\n", 
					i, z, b[i], s, err);
			}
		}
#endif
		i = z;
		if (0 == i%10) { printf ("."); fflush (stdout); }
	}
	printf ("\n");

	return 0;
}
