

/* quickie test case */

#include <math.h>
#include <stdio.h>
#include "bessel.h"

#define EPS 1.0e-15
int
main ()
{
	int i;
	double z, b[600];

#define ZRANGE 10.0
#define NRANGE 100

	/* perform broad tests */
	for (z=-ZRANGE; z<ZRANGE; z+=1.0)
	{
		bessel (NRANGE, z, b);

		s = j0 (z);
		err = fabs ((s-b[0])/s);
		if (EPS < err)
		{
			printf ("Bad j0 Test: i=%d z=%g me=%g sys=%g err=%g\n", 
				i, z, b[0], s, err);
		}

		s = j1 (z);
		err = fabs ((s-b[1])/s);
		if (EPS < err)
		{
			printf ("Bad j1 Test: i=%d z=%g me=%g sys=%g err=%g\n", 
				i, z, b[1], s, err);
		}

		for (i=0; i<NRANGE; i++) 
		{
			double s, err;
			s = jn (i,z);
			err = fabs ((s-b[i])/s);
			if (EPS < err)
			{
				printf ("Bad Test: i=%d z=%g me=%g sys=%g err=%g\n", 
					i, z, b[i], s, err);
			}
		}
		i = z;
		if (0 == i%10) { printf ("."); fflush (stdout); }
	}
	printf ("\n");

	return 0;
}
