/*
 * potts.c
 *
 * Potts model
 *
 * Linas September 2005
 *
 * State is encoded as a p-adic string.
 *
 * Actually, 2-adic in the current implementation
 */

/* Nearest neighbor interaction */
double nearest_neighbor (double s)
{
	double s0,s1;

	s0 = -1.0;
	if (s> 0.5) {
		s0 = 1.0;
		s -= 0.5;
	}
	s *= 2.0;
	
	s1 = -1.0;
	if (s> 0.5) s1 = 1.0;
	
	return 0.2 * s0*s1;
}

double energy (double (*func)(double), double s, int n, double temperature)
{
	int i;

	double en = 0.0;
	for (i=0; i<n; i++)
	{
		en += func (s);
		
		if (s> 0.5) s -= 0.5;
		s *= 2.0;
	}
}

