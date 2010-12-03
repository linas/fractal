
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
#include "mp-consts.h"
#include "mp-zeta.h"

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

	npts = 5600;
	
	mpf_set_d(step, 1.41333);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
#if 0
	mpf_sqrt(step, step);
	mpf_sqrt(step, step);
#endif
	mpf_ui_div(step, 1, step);

	mpf_set_ui(eps, 1);
	mpf_div_ui(eps, eps, 100);
	// mpf_set_d(eps, 1.0e-100);
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

void goldy(double scale, int prec)
{
	mpf_t golden;
	mpf_init(golden);

#if GOLDEN
	// Golden ratio
	// golden = 0.5*(sqrt(5.0) - 1.0);
	mpf_sqrt_ui(golden, 5);
	mpf_sub_ui(golden, golden, 1);
	mpf_div_ui(golden, golden, 2);

	/* Wow!  eps^0.28 provides an excellent fit. */
	// pow(feps, 0.2798); even ebetter!
	qint(golden, 0.2798, prec);
#endif

#if SILVER
	// Silver mean 0.2141 is a good fit
	mpf_sqrt_ui(golden, 2);
	mpf_sub_ui(golden, golden, 1);
	qint(golden, 0.2141, prec);
#endif

#if 0
	mpf_set_ui(golden, 1);
	mpf_div_ui(golden, golden, 12);
	qint(golden, scale, prec);
#endif 

	fp_pi(golden, prec);
	mpf_sub_ui(golden, golden, 3);

#if 0
	fp_e(golden, prec);
	mpf_sub_ui(golden, golden, 2);

	fp_half_sqrt_three(golden);

	fp_pi_half(golden, prec);
	mpf_sub_ui(golden, golden, 1);

	fp_sqrt_two_pi(golden, prec);
	mpf_sub_ui(golden, golden, 2);

	fp_log_two_pi(golden, prec);
   mpf_sub_ui(golden, golden, 1);

	fp_two_over_pi(golden, prec);
	fp_log2(golden, prec);
	fp_euler_mascheroni(golden, prec);

	fp_zeta_even(golden, 2, prec);
	fp_zeta(golden, 3, prec);
	fp_zeta(golden, 5, prec);
   mpf_sub_ui(golden, golden, 1);
#endif

	qint(golden, scale, prec);
}

int main (int argc, char * argv[])
{
	int prec, nbits;
	double scale;

	if (3 > argc)
	{
		fprintf(stderr, "Usage: %s <decimal-precision> <scale>\n", argv[0]);
		exit(1);
	}

	/* prec is decimal-places of precision */
	prec = 50;
	prec = atoi(argv[1]);
	scale = atoi(argv[2]);

	/* Set the precision (number of binary bits) */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits);

	goldy(scale, prec);
	
	return 0;
}
