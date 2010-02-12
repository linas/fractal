/* 
 * counting.c
 * Crude quick-n-dirty prime-couting attempt.
 * 
 * Linas Vepstas February 2010
 */

#include <stdio.h>
#include "question.h"

main()
{
	int nmax, k;
	int logstep = 10;
	int nmax = 100;

	nsteps = 1 << logstep;

	double acc = 0.0;
	double delta = 1.0 / ((double) nsteps);
	for (k=nsteps; k <= nmax*nsteps; k++)
	{
		double x = k * delta;
		double term  = question_mark(nsteps, k);
		acc += delta*term;
		printf("%d	%f %g	%g\n", k, x, term, acc);
	}
}

