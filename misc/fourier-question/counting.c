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
	long double x, term;
	long long k;
	// int logstep = 10;
	int logstep = 12;
	// int prtstep = logstep - 3;
	int prtstep = logstep+5;
	long long nmax = 10*1000*1000;

	long long nsteps = 1LL << logstep;
	long long psteps = 1LL << prtstep;

	long double lg2 = 1.0L / logl(2.0L);

	long double acc = 0.0L;
	long double li = 0.0L;
	long double delta = 1.0L / ((long double) nsteps);
#if 0
	for (k=1; k < nsteps; k++)
	{
		x = k * delta;
		term  = lg2 * 1.0L  / question_inverse(x);
		term = 0.0L;
		acc += delta*term;
		li += delta / logl(x);
		if (k%psteps == 0) printf("%Ld	%Lf	%Lg	%Lg	%Lg\n", k, x, term, acc, li);
	}
#endif 
	acc = -1.0L;
	for (k=nsteps+1LL; k <= nmax*nsteps; k++)
	{
		x = k * delta;
		long double ox = 1.0L / x;
		term  = 1.0L / question_inverse(ox);
		term -= 1.0L;
		term = lg2 / term;
		acc += delta*term;
		// li += delta / logl(x);
		if (k%psteps == 0)
		{
			printf("%Ld	%Lf	%Lg	%Lg	%Lg\n", k, x, term, acc, li);
			fflush (stdout);
		}
	}

	/* Print the last one no matter what */
	printf("%Ld	%Lf	%Lg	%Lg	%Lg\n", k, x, term, acc, li);
}

