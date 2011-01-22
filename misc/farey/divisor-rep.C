
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
	int wtf = GetPrime (pr, deno);
printf("wtf = %d\n", wtf);

	int i;
	for (i=1; i<deno; i++)
	{
		int num  = i;
		ContinuedFraction cf;
		cf.SetRatio (num,deno);
		double x = cf.ToFarey ();
		int nt = cf.GetNumTerms();
		int j;
		for (j=1; j<=nt; j++)
		{
			int p = cf.GetTerm(j);
			printf("ahh %d %d \n", j, p);
		}
		printf ("its %g\n", x);
	}
}
