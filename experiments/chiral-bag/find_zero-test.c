
/* quick and dirty test */

#include <math.h>
#include <stdio.h>

#include "find_zero.h"

int main ()
{
	double p;

	p = FindZero (sin, 0.0, 15, 0.001, 6.283, 100);

	printf ("found root w/ error=%g\n", p-M_PI);

	return 0;
}
