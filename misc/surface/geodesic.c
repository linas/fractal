
/*
 * geodesic.c
 *
 * Graphs of the lengths and energies of geodesics on the riemann 
 * surface/fundamental domain. For now,the symmetric orbits only.
 */

#include <math.h>

int get_n_of_rho (double rho)
{
	rho *= rho;
	rho = sqrt (rho-0.75);
	rho -= 0.5;
	rho = floor (rho);
	return (int) rho;
}

double geo_length (double rho)
{
	int n = get_n_of_rho (rho);

	double nu = rho*rho + n*(n-2);
	double de = -nu;
	double cr = 2.0*rho*(n+1);
	nu += cr;
	de += cr;

	double len = log (nu/de);
	return len;
}


main () 
{
	int i;

	int n = 20;
	double rho = 1.0;
	delta = 0.1;
	rho += delta;
	for (i=0; i<n; i++)
	{
		double len = geo_length (rho);

		printf ("%d	%g	%g\n", i, rho, len);
		rho += delta;
	}
}
