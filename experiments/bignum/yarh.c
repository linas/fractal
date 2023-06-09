/**
 * yarh.c
 * Port of generate/swap.C to bignum
 * actually, of yarh/yarh.c -- see also yarh/README
 *
 * FUNCTION:
 * Integral the permuation group of continued fractions
 * Expect to get Riemann zeta in the Gauss map case
 * and that is what we seem to get ... need high integration 
 * order though to get anything on the r=1/2 axis ... 
 *
 * Also includes integral of mobius ...
 *
 * Results written up in yarh.lyx
 *
 * Linas Feb 2005
 * Linas Vepstas December 2010
 * Linas Vepstas October 2015
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mp-complex.h"
#include "mp-misc.h"
#include "mp-trig.h"
#include "mp-zerofind.h"

/**
 * Swap the first and second digits of the continued fraction
 *
 * nprec == number of deciman places of precision.
 */
void swap_1_2 (mpf_t y, mpf_t x, int nprec)
{
	static int cur_prec = 0;
	static mpf_t epsi;
	if (cur_prec != nprec)
	{
		cur_prec = nprec;
		mpf_init(epsi);
		mpf_set_ui(epsi, 1);
		mpf_div_2exp(epsi, epsi, (int) (3.3 * ((double) nprec)));
	}

	/* a1 and a2 are the first two digitis of the 
	 * continued fraction */
	mpf_t ox, a1, a2;
	mpf_init (ox);
	mpf_init (a1);
	mpf_init (a2);

	mpf_ui_div (ox, 1, x);
	mpf_floor (a1, ox);
	mpf_sub(y, ox, a1);

	/* If 1/x-a_1 is zero, this means a2 is infinity.
	 * So exchanging a1 and a2 gives y=1/infty = 0 */
	if (0 < mpf_cmp (epsi, y))
		goto done;

	/* Sometimes, rounding errors above give 0.99999...
	 * Deal with these properly, as if above was a zero. */
	mpf_ui_sub (ox, 1, y);
	if (0 < mpf_cmp (epsi, ox))
	{
		mpf_set (y, ox);
		goto done;
	}

	/* Now get the second digit */
	mpf_ui_div (ox, 1, y);
	mpf_floor (a2, ox);
	mpf_sub(y, ox, a2);

	/* re-assemble the continued fraction */
	mpf_add(ox, y, a1);
	mpf_ui_div(y, 1, ox);
	mpf_add (ox, y, a2);
	mpf_ui_div(y, 1, ox);

done:
	mpf_clear (ox);
	mpf_clear (a1);
	mpf_clear (a2);
}

/**
 * The integrand, which is swap(x) * x^s 
 */
void integrand(cpx_t y, mpf_t x, cpx_t s, int nprec)
{
	mpf_t perm;
	mpf_init (perm);

	swap_1_2 (perm, x, nprec);

	cpx_mpf_pow (y, x, s, nprec);
	cpx_times_mpf (y, y, perm);

	mpf_clear (perm);
}

void test_parabola(cpx_t y, unsigned int nsteps, cpx_t s, int nprec)
{
	cpx_t cent;
	cpx_init (cent);
	cpx_set_d(cent, 0.4980812345, 18.313412345);
	cpx_sub(cent, cent, s);
	cpx_mul(y, cent, cent);
	cpx_clear (cent);
}

/**
 * Compute single integral of the integrand.
 * actually compute 
 * zeta = s/(s-1) - s \int_0^1 swap(x) x^{s-1} dx 
 *
 * Done using the simplest Newton sum possible.
 */
void integral(cpx_t y, unsigned int nsteps, cpx_t s, int nprec)
{
	int i;
// #define TESTING
#ifdef TESTING
test_parabola(y,nsteps,s,nprec);
return;
#endif

	mpf_t step, x;
	cpx_t term, ess, essm1;

	mpf_init (step);
	mpf_init (x);

	cpx_init (term);
	cpx_init (ess);
	cpx_init (essm1);

	cpx_set(ess, s);
	cpx_sub_ui (essm1, ess, 1, 0);

	/* Integration stepsize */
	mpf_set_ui (step, nsteps);
	mpf_ui_div (step, 1, step);

	/* initial value */
	mpf_div_ui (x, step, 2);
	mpf_neg (x, x);
	mpf_add_ui (x, x, 1);

	/* integration loop */
	cpx_set_ui (y, 0, 0);
	for (i=0; i<nsteps; i++)
	{
		integrand (term, x, essm1, nprec);
		cpx_add (y, y, term);

		mpf_sub (x, x, step);
	}

	/* Divide by the actual number of samples */
	cpx_div_ui (y, y, nsteps);

	/* integral times s */
	cpx_mul (y, y, ess);

	/* compute 1/(s-1) */
	cpx_recip (ess, essm1);
	cpx_sub (y, ess, y);

	/* s/(s-1) = 1/(s-1) + 1 so add 1 now */
	cpx_add_ui(y, y, 1, 0);

	mpf_clear (step);
	mpf_clear (x);
	cpx_clear (term);
	cpx_clear (ess);
	cpx_clear (essm1);
}

static unsigned int global_nsteps;
void integral_f(cpx_t y, cpx_t s, int nprec)
{
	integral(y, global_nsteps, s, nprec);
}
/* ========================================================= */
/*
	Some results: swap of 1 and 2:
*/

int main (int argc, char * argv[])
{
	unsigned int nsteps;
	int prec, nbits;

	if (argc < 3)
	{
		fprintf (stderr, "%s: <prec> <nsteps>\n", argv[0]);
		return -1;
	}

	prec = atoi(argv[1]);
	nsteps = atoi(argv[2]);

   /* Set the precision (number of binary bits) for calculations */
   nbits = 3.3*(prec + 8);
   mpf_set_default_prec (nbits);

#define ZERO_FINDER 
#ifdef ZERO_FINDER
	/* Set the precision to which we want the zero */
	int ndigits = 10;

	printf ("#\n# intermediate decimal precision = %d\n", prec);
	printf ("#\n# zero decimal precision = %d\n#\n", ndigits);
	fflush (stdout);

	cpx_t zero, s0, e1, e2;
	cpx_init(zero);
	cpx_init(s0);
	cpx_init(e1);
	cpx_init(e2);

	/* Initial guess */
	cpx_set_d (s0, 0.5, 15.0);

	/* Initial directions */
	cpx_set_d (e1, 0.05, 0.0);
	cpx_set_d (e2, 0.0, 0.1);

	double rr = sqrt(sqrt(2.0));
	while(1)
	{
		// printf ("#\n# num steps = %d\n", nsteps);
		global_nsteps = nsteps;
		int invalid = cpx_find_zero(zero, integral_f, s0, e1, e2, ndigits, prec);

		/* Print only valid results. */
		if (!invalid)
		{
			double re = cpx_get_re(zero);
			double im = cpx_get_im(zero);
			printf ("%d	%21.18g	%21.18g\n", nsteps, re, im);
			fflush (stdout);
		}

		nsteps ++;
#define LOG_STEPS
#ifdef LOG_STEPS
		nsteps = (int) (((double) nsteps * rr) + 3);
		nsteps += nsteps%2 + 1;  /* make it odd, always */
#endif
	}

	cpx_clear(zero);
#endif

// #define WALK_THE_LINE
#ifdef WALK_THE_LINE
	/* Walk up the imaginary axix at re=1/2, and print the results */
	printf ("#\n# decimal precision = %d\n", prec);
	printf ("#\n# num steps = %d\n#\n", nsteps);
	fflush (stdout);

	cpx_t y, s;
	cpx_init(y);
	cpx_init(s);

	double t = 0.0;
	int i = 1;
	while (t < 100)
	{
		cpx_set_d (s, 0.5, t);
		integral(y, nsteps, s, prec);

		double re = cpx_get_re(y);
		double im = cpx_get_im(y);

		printf("%d\t%g\t%g\t%g\n", i, t, re, im);
		fflush (stdout);
		t += 0.1;
		i++;
	}
#endif

	return 0;
}
