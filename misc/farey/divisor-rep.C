
/* divisor-rep.C
 *
 * Linas Vepstas January 2011
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"
#include "Prime.h"



main (int argc, char * argv[])
{
	int deno = atoi (argv[1]);
	
	struct Prime * pr = CreatePrime ();

	int i;
	for (i=1; i<deno; i++)
	{
		int num  = i;
		ContinuedFraction cf;
		cf.SetRatio (num,deno);
		double x = cf.ToFarey ();
		int nt = cf.GetNumTerms();

		double ber = 1.0;
		int j;
		for (j=1; j<=nt; j++)
		{
			int p = cf.GetTerm(j);
			int f = GetPrime (pr, j-1);

			while (p) { ber *= f; p--;  }
		}
		ber = log(ber);
		printf ("%d	%g	%g\n", i, x, ber);
	}
}
