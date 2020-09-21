/*
 * sidorov.c
 *
 * Understand bit-sequence mappings.  The expander and compresor
 * functions. This time with alternative expansions, per sidorov.
 *
 * Linas Vepstas Dec 2017; Sept 2020
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// Comp[ute the m from the sidorov paper. This is one more than the
// run of zeros.  It depends only on beta.  Note that:
// m=1 for K < 0.25 (1+sqrt(5)) = 0.80902
// m=2 up until about 0.878
// m=3 up until about 0.967
int emrun(double K)
{
	double beta = 2.0*K;
	double loga = (beta - 1.0) / (2.0-beta);
	loga = log(loga) / log(beta);
	loga = floor(loga) + 1.0;
	return (int) loga;
}

// Note that beta = 2*K
double sdr(double y, double K, int em)
{
	// Generate the beta expansion bits in a greedy fashion.
	char grebits[50];
	double greedy[50];
	for (int i=0; i<50; i++)
	{
		greedy[i] = y;
		if (0.5 <= y)
		{
			y -= 0.5;
			grebits[i] = 1;
		}
		else grebits[i] = 0;
		y *= 2.0*K;
	}

	// Search for em runs.
	// The Sidorov paper has an error, the m is off by one.
	char lobits[50];
	for (int i=0; i<50-em; i++)
	{
		lobits[i] = grebits[i];
		if (1 == grebits[i])
		{
			bool found = true;
			for (int j=1; j<=em; j++)
			{
				if (1 == grebits[i+j])
				{
					found = false;
					break;
				}
			}
			if (found)
			{
				lobits[i] = 0;
				y = greedy[i];
				printf("# got one at i=%d y=%g next=%g\n", i, y, y*2*K);
				y *= 2.0*K;

				for (int j=i+1; j<50; j++)
				{
					if (0.5 <= y)
					{
						y -= 0.5;
						lobits[j] = 1;
					}
					else lobits[j] = 0;
					y *= 2.0*K;
				}
				break;
			}
		}
	}

	double Jay = K;

	// Reconstruct both sequences;
	double hiacc = 1.0e-30;
	double loacc = 1.0e-30;
	for (int i=0; i<50; i++)
	{
		hiacc *= 1.0 / (2.0*Jay);
		loacc *= 1.0 / (2.0*Jay);
		if (grebits[50-i-1])
		{
			hiacc += 0.5;
		}
		if (lobits[50-i-1])
		{
			loacc += 0.5;
		}
	}

	printf("# hi=");
	for (int i=0; i<50; i++) printf("%d", grebits[i]);
	printf("\n");
	printf("# lo=");
	for (int i=0; i<50; i++) printf("%d", lobits[i]);
	printf("\n");

	return hiacc-loacc;
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);

	int em = emrun(Kay);
	printf("#\n# K=%g m=%d\n#\n", Kay, em);

#if 1
	int npts = 313;
	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		double z = sdr (x, Kay, em);
		printf("%d	%g	%g\n", i, x, z);
	}
#endif
}
