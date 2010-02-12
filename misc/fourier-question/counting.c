/* 
 * sym.c
 * Attempt to discover symmetries of integral of question mark
 * 
 * Linas Vepstas February 2010
 */

#include <stdio.h>
#include "question.h"

main()
{
	int nmax, k;
	int n = 10;

	nmax = 1 <<n;

	double acc = 0.0;
	double delta = 1.0 / ((double) nmax);
	for (k=0; k<=nmax; k++)
	{
		double x = k * delta;
		double term  = question_mark(k,nmax);
		acc += delta*term;
		printf("%d	%f %g	%g\n", k, x, term, acc);
	}
}

