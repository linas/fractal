

/* quickie test case */

#include <stdio.h>
#include "bessel.h"

#define EPS 1.0e-15
main ()
{
	int i;
	double z, b[600];


	for (z=-1000.0; z<1000.0; z+=1.0)
	{
		bessel (500, z, b);
		for (i=0; i<500; i++) 
		{
			double s, err;
			s = jn (i,z);
			err = fabs ((s-b[i])/s);
			if (EPS < err)
			{
				printf ("Bad Test: i=%d z=%g me=%g sys=%g\n", i, z, b[i], s);
			}
		}
		i = z;
		if (0 == i%10) { printf ("."); fflush (stdout); }
	}
	printf ("\n");
}
