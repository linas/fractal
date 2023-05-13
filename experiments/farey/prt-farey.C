
/* prt-farey.C
 *
 * Linas Vepstas February 2004
 */

#include <stdlib.h>

#include "Farey.h"

main (int argc, char * argv[])
{
	int num = atoi (argv[1]);
	int deno = atoi (argv[2]);
	
	ContinuedFraction cf;
	cf.SetRatio (num,deno);
	double x = cf.ToFarey ();
	printf ("its %g\n", x);
	// cf.Print ();
}
