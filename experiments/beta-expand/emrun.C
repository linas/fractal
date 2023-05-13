/*
 * emrun.C
 *
 * Linas Vepstas Sept 2020
 */

#include <math.h>

// Compute the m from the Sidorov paper. This is the length of the
// run of zeros we need to see, before exploring an alternate branch.
// It depends only on beta. Note that the Sidorov paper incorrectly
// states that the run of zeros is one less than m. This is false.
// The required run  of zeros has to be at least m.
// Note that:
// m=1 for K < 0.25 (1+sqrt(5)) = 0.8090169943749475
// m=2 up until about 0.878
// m=3 up until about 0.967
int emrun(double K)
{
	double beta = 2.0*K;
	double gold = 0.5 * (1.0 + sqrt(5));
	if (beta <= gold) return 1;

	double loga = (beta - 1.0) / (2.0-beta);
	loga = log(loga) / log(beta);
	loga = floor(loga) + 1.0;
	return (int) loga;
}

// ================================================================
