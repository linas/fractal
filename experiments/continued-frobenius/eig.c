
#include <math.h>

main ()
{
	double lam[6];

	lam[0] = 1.0;
	lam[1] = 0.30366;
	lam[2] = 0.10086;
	lam[3] = 0.0354;
	lam[4] = 0.01265;
	lam[5] = 0.00445;

	int k;

	double pr = 0.0;
	for (k=0; k<6; k++)
	{
		double ll = log (lam[k]);
		printf ("its %d %g %g %g\n", k, ll, ll-pr, lam[k-1]/lam[k]);
		pr = ll;
	}
}
