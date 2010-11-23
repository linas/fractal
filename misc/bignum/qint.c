
/*
 * Integral of question mark measure.
 *
 * Linas March 2010
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mp-quest.h"

void qint(mpf_t center, double scale, int prec)
{
	int i, npts;
	double fmean, qfmean;

	mpf_t x, eps, step, qlo, qmid, qhi;
	mpf_init(x);
	mpf_init(step);
	mpf_init(eps);
	mpf_init(qlo);
	mpf_init(qmid);
	mpf_init(qhi);

	question_mark(qmid, center, prec);

	fmean = mpf_get_d(center);
	qfmean = mpf_get_d(qmid);
	printf("#\n# Mean=%g ?(Mean)=%g\n", fmean, qfmean); 

	npts = 1600;
	
	mpf_set_d(step, 1.41333);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
#if 1
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
#endif
	mpf_ui_div(step, 1, step);

	mpf_set_ui(eps, 1);
	mpf_div_ui(eps, eps, 100);
	// mpf_set_d(eps, 1.0e-200);
 
	for (i=0; i<npts; i++)
	{
		double feps, flo, fhi, fdiff, r;

		mpf_sub(x, center, eps);
		question_mark(qlo, x, prec);

		mpf_add(x, center, eps);
		question_mark(qhi, x, prec);

		mpf_sub(qlo, qmid, qlo);
		mpf_sub(qhi, qhi, qmid);

		mpf_div(qlo, qlo, eps);
		mpf_div(qhi, qhi, eps);

		feps = mpf_get_d(eps);
		flo = mpf_get_d(qlo);
		fhi = mpf_get_d(qhi);
		fdiff = fhi-flo;

#if 1
		r = pow(feps, scale);
		flo *= r;
		fhi *= r;
		fdiff *= r;
#endif
		printf("%d	%g	%g	%g	%g\n", i, feps, flo, fhi, fdiff);
		fflush(stdout);

		mpf_mul(eps, eps, step);
	}
}

void goldy(int prec)
{
	mpf_t golden;
	mpf_init(golden);

	// mpf_set_ui (x, 1);
	// mpf_div_ui (x, x, 4);

	// Golden ratio
	// golden = 0.5*(sqrt(5.0) - 1.0);
	mpf_sqrt_ui(golden, 5);
	mpf_sub_ui(golden, golden, 1);
	mpf_div_ui(golden, golden, 2);

	/* Wow!  eps^0.28 provides an excellent fit. */
	// pow(feps, 0.2798); even ebetter!
	// qint(golden, 0.2798, prec);

	// Silver mean
	mpf_sqrt_ui(golden, 2);
	mpf_sub_ui(golden, golden, 1);
	qint(golden, 0.213, prec);
}

int main (int argc, char * argv[])
{
	int prec, nbits;

	if (2 > argc)
	{
		fprintf(stderr, "Usage: %s <decimal-precision>\n", argv[0]);
		exit(1);
	}

	/* prec is decimal-places of precision */
	prec = 50;
	prec = atoi(argv[1]);

	/* Set the precision (number of binary bits) */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits);

	goldy(prec);
	
	return 0;
}
