/*
 * FUNCTION:
 * return greatest prime factor (larged prime divisor)
 *
 * HISTORY:
 * April 2016 -- linas
 */

#include "gpf.h"
#include "prime.h"

/* ------------------------------------------------------------ */
/** Return the greatest prime factor.
 */
unsigned long gpf (unsigned long n)
{
	for (unsigned int nth = 1; ; nth++)
	{
		unsigned long p = get_nth_prime(nth);
		while (n % p == 0)
		{
			n /= p;
		}
		if (n < p*p) return n;
	}

	return 0;
}
