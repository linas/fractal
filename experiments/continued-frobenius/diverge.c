
/* diverge.c
 *
 * validate some simple divergent sums
 * Linas Dec 2003
 */

#include <math.h>


main ()
{
	double t = 1.0;

	for (t=1.0; t>1e-6; t *= 0.5)
	{
		int k;

		double acc = 0.0;
		double s = 1.0;
		for (k=0; k<100000; k++)
		{
			double kay = (double) k;
			double reg = exp (-t*t*kay*kay);
			if (reg < 1e-16) break;
			acc += s* (kay+6.0)*(kay+5.0)*(kay+4.0)*(kay+3.0)*(kay+2.0) * reg;
			// acc += s* (kay+5.0)*(kay+4.0)*(kay+3.0)*(kay+2.0) * reg;
			// acc += s* (kay+4.0)*(kay+3.0)*(kay+2.0) * reg;
			// acc += s* (kay+3.0)*(kay+2.0) * reg;
			// acc += s* reg;
			s = -s;
		}
		printf ("t=%f acc = %f\n", t, acc);
	}
}
