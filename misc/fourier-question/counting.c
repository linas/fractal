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
	int k;
	int logstep = 10;
	int prtstep = logstep - 3;
	int nmax = 1000;

	int nsteps = 1 << logstep;
	int psteps = 1 << prtstep;

	double acc = 0.0;
	double delta = 1.0 / ((double) nsteps);
	for (k=1; k <= nsteps; k++)
	{
		double x = k * delta;
		double term  = 1.0 / question_inverse(x);
		acc += delta*term;
		if (k%psteps == 0) printf("%d	%f %g	%g\n", k, x, term, acc);
	}
	for (k=nsteps+1; k <= nmax*nsteps; k++)
	{
		double x = k * delta;
		double ox = 1.0 / x;
		double term  = question_inverse(ox);
		acc += delta*term;
		if (k%psteps == 0) printf("%d	%f %g	%g\n", k, x, term, acc);
	}
}

