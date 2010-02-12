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
	int prtstep = logstep;
	// int prtstep = logstep+5;
	// long long nmax = 10*1000*1000;
	long long nmax = 100*1000;

	long long nsteps = 1LL << logstep;
	long long psteps = 1LL << prtstep;

	long double lg2 = 1.0L / logl(2.0L);

	long double rqm1 = 0.0L;
	long double qacc = 0.0L;
	long double rqp1 = 0.0L;
	long double li = 0.0L;
	long double delta = 1.0L / ((long double) nsteps);

#if 1
	for (k=1; k < nsteps; k++)
	{
		x = k * delta;
		// term  = lg2 * 1.0L  / question_inverse(x);
		// term = 0.0L;
		// acc += delta*term;
		li += delta / logl(x);
		// if (k%psteps == 0) printf("%Ld	%Lf	%Lg	%Lg	%Lg\n", k, x, term, acc, li);
	}
#endif 


	printf("#\n# Integrals of inverse question mark\n#\n");
	printf("#\n# columns: k, x, li, rqm1, qi, rqp1\n");
	printf("#\n# li == log int, rqm1= ln2/(1/?-1), qi=?, rqp1= ln2/(1/?+1)\n");

	for (k=nsteps+1LL; k <= nmax*nsteps; k++)
	{
		x = k * delta;
		long double ox = 1.0L / x;
		long double qi = question_inverse(ox);
		qacc += lg2 * delta*qi;

		// term  = lg2 / ((1.0L / qi) - 1.0L);
		term  = lg2 * qi / (1.0L - qi);
		rqm1 += delta*term;

		term  = lg2 * qi / (1.0L + qi);
		rqp1 += delta*term;

		li += delta / logl(x);
		if (k%psteps == 0)
		{
			printf("%Ld	%Lf	%Lg	%Lg	%Lg	%Lg\n", k, x, li, rqm1, qacc, rqp1);
			fflush (stdout);
		}
	}

	/* Print the last one no matter what */
	printf("%Ld	%Lf	%Lg	%Lg	%Lg	%Lg\n", k, x, li, rqm1, qacc, rqp1);
}

