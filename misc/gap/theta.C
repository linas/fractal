
/* 
 * theta.C
 *
 * Generate the line segmments of the hyperbolic curves.
 */
#include <math.h>
#include <stdio.h>

double theta (double x, double *next)
{
	if (0.5 > x) { *next = 0.5; return (0.5*x); }
	if (0.75 > x) { *next = 0.75; return (x-0.25); }
	*next = 1.0;
	return (2.0*x -1.0);
}

void
recurse (int cnt, double (*fn)(double, double*), double pos, double x)
{
	double next;
	double y = fn (x, &next);
	if (1 == cnt)
	{
		printf ("%8.6g	%8.6g\n", pos, y);
		if (1.0 == y) return;
		return;
	}

printf ("duude cnt=%d pos=%g x=%g y=%g next=%g\n", cnt, pos, x, y, next);

	if (y < next)
	{
		recurse (cnt-1, fn, pos, y);
	}

	recurse (cnt, fn, pos, next);
	return;
}

void draw (void)
{
	recurse (3, theta, 0.0, 0.0);
}

main ()
{

	draw();
}
