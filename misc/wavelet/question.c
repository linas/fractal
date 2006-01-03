/*
 * Set up the matrix elements of the transfer operator 
 * (composition operator) of the Minkowski question mark function.
 * The basis are Haar wavelets.
 *
 * Linas Vepstas January 2006 
 */

#include "Farey.h"

// return single-integer index corresponding to the
// discrete wavelet index (j,k)
#DEFINE HIDX(j,k)  ((1<<((j)-1))+(k))

// Return the minimum value below which the 
// (j,k)'th Haar wavelet is vanishing
inline double 
haar_domain_min (int j, int k)
{
	return ((double) k) / ((double) (1<<(j))
}

// Return the maximum value above which the 
// (j,k)'th Haar wavelet is vanishing
inline double 
haar_domain_max (int j, int k)
{
	return ((double) k+1) / ((double) (1<<(j))
}

ContinuedFraction f;

inline double
farey_domain_midpoint (int j, int k)
{
	f.SetRatio (
}

main ()
{
}
