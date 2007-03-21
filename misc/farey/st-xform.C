
/* st-xform.C
 *
 * Linas Vepstas March 2007
 */

#include <stdlib.h>

#include "Farey.h"

main (int argc, char * argv[])
{
	int npts = atoi (argv[1]);
	
	ContinuedFraction cf, st;

	int n;
	for (n=1; n<npts; n++)
	{
		int num = n;
		int deno = npts;
		cf.SetRatio (num, deno);

		int intpart = cf.GetTerm (0);
		st.SetTerm (0, intpart);
		
		int r=1, t=1;
		int a = cf.GetTerm (t);
		a--;
		while (34 > r)
		{
			while (0 < a)
			{
				st.SetTerm (r, 2);
				r++;
				a--;
			}
			t++;
			a = cf.GetTerm (t);
			st.SetTerm (r, a);
			r++;
			t++;
			a = cf.GetTerm (t);
		}
		double y = st.ToReal();
		double x = ((double) num)/ ((double) deno);
		printf ("%d\t%g\t%g\n", n, x, y);
	}
}
