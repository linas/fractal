/* 
 * counting.c
 * Crude quick-n-dirty prime-couting attempt.
 * Results are fantastic! Wow!
 * 
 * Linas Vepstas February 2010
 */

#include <math.h>
#include <stdio.h>
#include "question.h"

main()
{
	double x, term;
	int k;
	// int logstep = 10;
	int logstep = 7;
	// int prtstep = logstep - 3;
	int prtstep = logstep+5;
	int nmax = 10*1000*1000;

	int nsteps = 1 << logstep;
	int psteps = 1 << prtstep;

	double lg2 = 1.0 / log(2.0);

	double acc = 0.0;
	double li = 0.0;
	double delta = 1.0 / ((double) nsteps);
	for (k=1; k < nsteps; k++)
	{
		x = k * delta;
		term  = lg2 * 1.0  / question_inverse(x);
		term = 0.0;
		acc += delta*term;
		li += delta / log(x);
		if (k%psteps == 0) printf("%d	%f %g	%g	%g\n", k, x, term, acc, li);
	}
	acc = -1.0;
	for (k=nsteps+1; k <= nmax*nsteps; k++)
	{
		x = k * delta;
		double ox = 1.0 / x;
		term  = 1.0 / question_inverse(ox);
		term -= 1.0;
		term = lg2 / term;
		acc += delta*term;
		li += delta / log(x);
		if (k%psteps == 0) printf("%d	%f %g	%g	%g\n", k, x, term, acc, li);
	}

	/* Print the last one no matter what */
	printf("%d	%f %g	%g	%g\n", k, x, term, acc, li);
}

