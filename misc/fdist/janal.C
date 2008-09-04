/*
 * janal.C
 * 
 * Fourier transform of
 * Distribution of the Farey Numbers on the unit interval
 * From hypothesized first principles.
 *
 * Linas October 2004
 * Linas Sept 2008
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void coeffs (int freq_max)
{
	int k, n, m;

	for (m=0; m< freq_max; m++)
	{
		double jm = 1.0;
		double tp = 1.0;
		int tn = 1.0;
		for(n=0; n<10; n++)
		{
			double term = 0.0;
			for(k=0; k<tn; k++)
			{
				term += cos(M_PI*m*(2*k+1) * tp);
			}
printf ("duude m=%d n=%d term=%g\n", m, n, term);
			jm += tp * term;
			tn *= 2;
			tp *= 0.5;
		}

		printf ("%d	%8.6g\n", m, jm);
		fflush (stdout);
	}

}

main(int argc, char *argv[])
{
	int i;

	if (argc < 2)
	{
		fprintf (stderr, "Usage: %s <freq>\n", argv[0]);
		exit (1);
	}
	int max_freq = atoi (argv[1]);

	coeffs (max_freq);
}

