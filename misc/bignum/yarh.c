/**
 * yarh.c
 * Port of swap.C to bignum
 *
 * FUNCTION:
 * Integral the permuation group of continued fractions
 * Expect to get Riemann zeta in the Gauss map case
 * and that is what we seem to get ... need high integration 
 * order though to get anything on the r=1/2 axis ... 
 *
 * Results written up in yarh.lyx
 *
 * Linas Feb 2005
 * Linas Vepstas December 2010
 */

#include <stdio.h>
#include <stdlib.h>

#include "mp-complex.h"
#include "mp-trig.h"
#include "mp-misc.h"

/**
 * Swap the first and second digits of the continued fraction
 *
 * nprec == number of deciman places of precision.
 */
void swap_1_2 (mpf_t y, mpf_t x, int nprec)
{
	static int init = 0;
	static mpf_t zero;
	if (!init)
	{
		init = 1;
		mpf_init(zero);
		mpf_set_ui(zero, 0);
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
	if (mpf_eq (y, zero, 3.32*nprec))
		goto done;

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
#if TESTING
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
	cpx_div_ui (y, y, nsteps-1);

	/* integral times s */
	cpx_mul (y, y, ess);

	/* compute 1/(s-1) */
	cpx_recip (ess, essm1);
	cpx_sub (y, y, ess);

	/* s/(s-1) = 1/(s-1) + 1 so add 1 now */
	cpx_add_ui(y, y, 1, 0);

	mpf_clear (step);
	mpf_clear (x);
	cpx_clear (term);
	cpx_clear (ess);
	cpx_clear (essm1);
}

/* ---------------------------------------------- */

/*
 * Find bottom of parabola.  Given three points,
 * fit a parabola to that, then find the minimum.
 * 
 * return the bottom in "loc"
 */
void quad_min(mpf_t loc, mpf_t a, mpf_t b, mpf_t c,
              mpf_t fa, mpf_t fb, mpf_t fc)
{
	mpf_t ba, bc, fba, fbc, deno, numer;
	mpf_init (ba);
	mpf_init (bc);
	mpf_init (fba);
	mpf_init (fbc);
	mpf_init (deno);
	mpf_init (numer);

	/* differences */
	mpf_sub (ba, b, a);
	mpf_sub (bc, b, c);
	mpf_sub (fba, fb, fa);
	mpf_sub (fbc, fb, fc);

	/* cross product */
	mpf_mul(fbc, fbc, ba);
	mpf_mul(fba, fba, bc);

	/* denominator */
	mpf_sub(deno, fbc, fba);

	/* cross again */
	mpf_mul(fbc, fbc, ba);
	mpf_mul(fba, fba, bc);

	/* numerator */
	mpf_sub(numer, fbc, fba);
	mpf_div (numer, numer, deno);
	mpf_div_ui(numer, numer, 2);

	mpf_sub(loc, b, numer);

	mpf_clear (ba);
	mpf_clear (bc);
	mpf_clear (fba);
	mpf_clear (fbc);
	mpf_clear (deno);
	mpf_clear (numer);
}

void interp(mpf_t yoc, mpf_t xoc, mpf_t a, mpf_t b, mpf_t ya, mpf_t yb)
{
	mpf_t ba, yba;

	mpf_init(ba);
	mpf_init(yba);

	mpf_sub (ba, b, a);
	mpf_sub (yba, yb, ya);
	mpf_div (yba, yba, ba);
	
	mpf_sub (ba, b, xoc);
	mpf_mul (yba, yba, ba);
	mpf_add (yoc, yba, ya);

	mpf_clear(ba);
	mpf_clear(yba);
}

void find_zero(int nsteps, int prec)
{
	mpf_t fa, fb, fc;
	mpf_init (fa);
	mpf_init (fb);
	mpf_init (fc);

	cpx_t sa, sb, sc;
	cpx_t ya, yb, yc;
	cpx_init (sa);
	cpx_init (sb);
	cpx_init (sc);

	cpx_init (ya);
	cpx_init (yb);
	cpx_init (yc);

	/* Initial guess */
	cpx_set_d (sa, 0.49, 18.4);
	cpx_set_d (sb, 0.5, 18.3);
	cpx_set_d (sc, 0.51, 18.2);

	integral (ya, nsteps, sa, prec);
	integral (yb, nsteps, sb, prec);
	integral (yc, nsteps, sc, prec);

	cpx_abs(fa, ya);
	cpx_abs(fb, yb);
	cpx_abs(fc, yc);

#if 0
cpx_prt("ya = ", ya); printf("\n");
cpx_prt("yb = ", yb); printf("\n");
cpx_prt("yc = ", yc); printf("\n");
#endif

	mpf_t ioc, roc;
	mpf_init (ioc);
	mpf_init (roc);

	int i;

	/* This loop finds zero, adding one bit of accuracy
	 * per loop iteration. */
	for (i=0; i<80; i++)
	{
		/* First, do the imaginary */
		quad_min (ioc, sa[0].im, sb[0].im, sc[0].im,
              	fa, fb, fc);
	
		interp (roc, ioc, sa[0].im, sb[0].im, sa[0].re, sb[0].re);
		double zero_im = mpf_get_d(ioc);

		/* Throw away the farther point */
		if (0 < mpf_cmp(fb, fc))
		{
			// printf("duude b biggest\n");
			cpx_set(sb, sa);
			cpx_set(yb, ya);
			mpf_set(fb, fa);
		}
		else
		{
			// printf("duude c biggest\n");
			cpx_set(sc, sa);
			cpx_set(yc, ya);
			mpf_set(fc, fa);
		}

		cpx_set_mpf(sa, roc, ioc);
		integral (ya, nsteps, sa, prec);
		cpx_abs(fa, ya);

		/* Next the real */
		quad_min (roc, sa[0].re, sb[0].re, sc[0].re,
              	fa, fb, fc);

		interp (ioc, roc, sa[0].re, sb[0].re, sa[0].im, sb[0].im);
		double zero_re = mpf_get_d(roc);

		/* Throw away the farther point */
		if (0 < mpf_cmp(fb, fc))
		{
			// printf("duude b biggest\n");
			cpx_set(sb, sa);
			cpx_set(yb, ya);
			mpf_set(fb, fa);
		}
		else
		{
			// printf("duude c biggest\n");
			cpx_set(sc, sa);
			cpx_set(yc, ya);
			mpf_set(fc, fa);
		}

		cpx_set_mpf(sa, roc, ioc);
		integral (ya, nsteps, sa, prec);
		cpx_abs(fa, ya);

		cpx_abs(roc, ya);
		double mini = mpf_get_d(roc);
		printf("%d location= %15.12g +i %15.12g  min=%g\n", 
			i, zero_re, zero_im, mini);
	}

	cpx_clear (sa);
	cpx_clear (sb);
	cpx_clear (sc);

	cpx_clear (ya);
	cpx_clear (yb);
	cpx_clear (yc);

	mpf_clear (fa);
	mpf_clear (fb);
	mpf_clear (fc);

	mpf_clear (ioc);
	mpf_clear (roc);
}

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

   /* Set the precision (number of binary bits) */
   nbits = 3.3*prec + 10;
   mpf_set_default_prec (nbits);

	find_zero(nsteps, prec);

#if WALK_THE_LINE
	printf ("#\n# decimal precision = %d\n", prec);
	printf ("#\n# num steps = %d\n#\n", nsteps);

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
