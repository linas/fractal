
/* 
 * theta.C
 *
 * Generate the line segmments of the hyperbolic curves.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double theta (double x)
{
	if (0.5 > x) { return (0.5*x); }
	if (0.75 > x) { return (x-0.25); }
	return (2.0*x -1.0);
}

void
recurse (int cnt, double (*fn)(double), double pos, double val)
{
	double next;
	double y = fn (val);
	if (1 == cnt)
	{
		printf ("%8.6g	%8.6g	%8.6g\n", pos, y, (1.0-pos)*y);
		return;
	}
	recurse (cnt-1, fn, pos, y);
}

void draw (int nrec)
{
	int i;
	recurse (nrec, theta, 0.0, 0.0);
	recurse (nrec, theta, 0.5, 0.5);
	recurse (nrec, theta, 0.75, 0.75);
	double delt = 0.125;
	double x = 0.75;
	for (i=1; i<nrec; i++) {
		x += delt;
		recurse (nrec, theta, x,x);
		delt *= 0.5;
	}
	recurse (nrec, theta, 1.0, 1.0);
}

main (int argc, char*argv[])
{
	int nrec=atoi(argv[1]);

	draw(nrec);
}
