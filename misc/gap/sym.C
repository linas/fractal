
/* sym.C
 *
 * Explore symmetries
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"

main (int argc, char *argv[])
{
	double z;
	int i;
	ContinuedFraction f;
	f.SetEvenize();

	int m = -2;
	int n = 5;
	f.SetRatio (m,n);
	double fact=f.ToFarey();
	printf ("# norm=%g\n", fact);
	

	int nmax = 23;
	for (i=0; i<nmax; i++)
	{

		int p = i;
		int q = nmax;
		double x = ((double) p)/ ((double) q);
		
		f.SetRatio (p,q);
		
		double y = f.ToFarey();
		
		// f.SetRatio (p+q,p+2*q);
		// z = 4.0* f.ToFarey() -2.0;
		
		// f.SetRatio (p+2*q,p+3*q);
		// z = 8.0* f.ToFarey()-6.0;
		
		// f.SetRatio (p+3*q,p+4*q);
		// z = 16.0* f.ToFarey()-14.0;
		
		// f.SetRatio (q,(n+1)*q-p);
		// z = f.ToFarey();
		// z *= (double) (1<<n);
		// z -= 1.0;

#if SURPRISE
		f.SetRatio ((m+1)*p+q, (2*m+1)*p+2*q);
		z = f.ToFarey();
		// y /= (1<<(2+m));
		// y +=0.5;
		z -= 0.5;
		z *= (1<<(2+m));
#endif
		
		f.SetRatio (m*p+q, (2*m+1)*p+2*q);
		z = f.ToFarey();
		z -= 0.5;
		z *= (1<<(2+m));
		z = -z;
		
		printf("%5d	%8.6g	%8.6g	%8.6g	%8.6g\n", i,x,y,z,y-z);

	}
}

