
/* 
 * FUNCTION:
 * return greatest common factor (greatest common divisor)
 *
 * HISTORY:
 * Oct 2004 -- linas
 */

#include "gcf.h"

/* ------------------------------------------------------------ */
/* Return the greatest common factor, 32-bit in accurate */
int
gcf32 (int nume, int denom)
{
	int t;
	t = nume % denom;
	nume = denom;
	denom = t;

	/* Euclids algorithm for obtaining the gcf */
	while (0 != denom)
	{
		t = nume % denom;
		nume = denom;
		denom = t;
	}

	/* nume now holds the GCD (Greatest Common Divisor) */
	return nume;
}

