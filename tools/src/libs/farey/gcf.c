/*
 * FUNCTION:
 * return greatest common factor (greatest common divisor)
 *
 * HISTORY:
 * Oct 2004 -- linas
 */

#include "gcf.h"

/* ------------------------------------------------------------ */
/**
 * Return the greatest common factor, 64-bit accurate.
 * Implements Euclid's algorithm.
 */
unsigned long
gcf64 (unsigned long nume, unsigned long denom)
{
	unsigned long t;
	t = nume % denom;
	nume = denom;
	denom = t;

	/* Euclid's algorithm for obtaining the gcf */
	while (0 != denom)
	{
		t = nume % denom;
		nume = denom;
		denom = t;
	}

	/* nume now holds the GCD (Greatest Common Divisor) */
	return nume;
}

unsigned long
lcm64 (unsigned long a, unsigned long b)
{
	return a * (b / gcf64(a,b));
}
