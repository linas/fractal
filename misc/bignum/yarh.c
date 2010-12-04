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
#define TESTING
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

/* Interpolate x between points a dn b, returning y */
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

/* Swap a pair of points */
#define SWAP(s,y,f,t,w,g) { \
	cpx_set(se, s); \
	cpx_set(ye, y); \
	mpf_set(fe, f); \
	 \
	cpx_set(s, t); \
	cpx_set(y, w); \
	mpf_set(f, g); \
	 \
	cpx_set(t, se); \
	cpx_set(w, ye); \
	mpf_set(g, fe); \
}

/* Sort, so that the fa < fb < fc */
#define SORT \
	if (0 < mpf_cmp(fa, fb)) \
	{ \
		SWAP(sa, ya, fa, sb, yb, fb); \
	} \
	if (0 < mpf_cmp(fb, fc)) \
	{ \
		SWAP(sc, yc, fc, sb, yb, fb); \
	} \
	if (0 < mpf_cmp(fa, fb)) \
	{ \
		SWAP(sa, ya, fa, sb, yb, fb); \
	}

/* Crude zero-finder */
void find_zero_borken(int nsteps, int prec)
{
	mpf_t fa, fb, fc, fd, fe;
	mpf_init (fa);
	mpf_init (fb);
	mpf_init (fc);
	mpf_init (fd);
	mpf_init (fe);

	cpx_t sa, sb, sc, sd, se;
	cpx_t ya, yb, yc, yd, ye;
	cpx_init (sa);
	cpx_init (sb);
	cpx_init (sc);
	cpx_init (sd);
	cpx_init (se);

	cpx_init (ya);
	cpx_init (yb);
	cpx_init (yc);
	cpx_init (yd);
	cpx_init (ye);

	/* Initial guess */
	cpx_set_d (sa, 0.5, 18.3);
	cpx_set_d (sb, 0.495, 18.4);
	cpx_set_d (sc, 0.505, 18.2);

	integral (ya, nsteps, sa, prec);
	integral (yb, nsteps, sb, prec);
	integral (yc, nsteps, sc, prec);

	cpx_abs(fa, ya);
	cpx_abs(fb, yb);
	cpx_abs(fc, yc);

	/* Sort, so that the fa < fb < fc */
	SORT;

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
	for (i=0; i<20; i++)
	{
		/* First, do the imaginary */
		quad_min (ioc, sa[0].im, sb[0].im, sc[0].im,
              	fa, fb, fc);
		interp (roc, ioc, sa[0].im, sb[0].im, sa[0].re, sb[0].re);
		cpx_set_mpf(sd, roc, ioc);

		/* standard eval from here on */
		integral (yd, nsteps, sd, prec);
		cpx_abs(fd, yd);

		if (0 < mpf_cmp(fd, fc))
		{
			/* fd is a worse estimate than the previous three.
			 * Discard, try again.  */
			printf("bad guess! Walking!\n");
			cpx_sub(sd, sa, sc);
			cpx_div_ui(sd, sd, 2);
			cpx_add(sd, sd, sa);

			integral (yd, nsteps, sd, prec);
			cpx_abs(fd, yd);

			if (0 < mpf_cmp(fd, fc))
			{
				printf("bad again! Exit\n");
				exit(-1);
			}
			else
			{
				/* fd's closer than fc */
				SWAP (sd, yd, fd, sc, yc, fc);
				SORT;
			}
		}
		else
		{
			/* fd's closer than fc */
			SWAP (sd, yd, fd, sc, yc, fc);
			SORT;
		}

		/* Next the real */
		quad_min (roc, sa[0].re, sb[0].re, sc[0].re,
              	fa, fb, fc);
		interp (ioc, roc, sa[0].re, sb[0].re, sa[0].im, sb[0].im);
		cpx_set_mpf(sd, roc, ioc);

		/* standard eval from here on. */
		integral (yd, nsteps, sd, prec);
		cpx_abs(fd, yd);

		if (0 < mpf_cmp(fd, fc))
		{
			/* fd is a worse estimate than the previous three.
			 * Discard, try again.  */
			printf("bad guess! Walking!\n");
			cpx_sub(sd, sa, sc);
			cpx_div_ui(sd, sd, 2);
			cpx_add(sd, sd, sa);

			integral (yd, nsteps, sd, prec);
			cpx_abs(fd, yd);

			if (0 < mpf_cmp(fd, fc))
			{
				printf("bad again! Exit\n");
				exit(-1);
			}
			else
			{
				/* fd's closer than fc */
				SWAP (sd, yd, fd, sc, yc, fc);
				SORT;
			}
		}
		else
		{
			/* fd's closer than fc */
			SWAP (sd, yd, fd, sc, yc, fc);
			SORT;
		}

		double zero_re = cpx_get_re(sa);
		double zero_im = cpx_get_im(sa);
		double mini = mpf_get_d(fa);
		printf("%d location= %15.12g +i %15.12g  min=%g\n", 
			i, zero_re, zero_im, mini);
	}

	cpx_clear (sa);
	cpx_clear (sb);
	cpx_clear (sc);
	cpx_clear (sd);
	cpx_clear (se);

	cpx_clear (ya);
	cpx_clear (yb);
	cpx_clear (yc);
	cpx_clear (yd);
	cpx_clear (ye);

	mpf_clear (fa);
	mpf_clear (fb);
	mpf_clear (fc);
	mpf_clear (fd);
	mpf_clear (fe);

	mpf_clear (ioc);
	mpf_clear (roc);
}

/* =============================================== */
/* Powell's method, approximaely. */

void find_zero(int nsteps, int prec)
{
	mpf_t zero, epsi;
	mpf_init (zero);
	mpf_init (epsi);

	/* compute tolerance */
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.3*prec));

	cpx_t s0, s1, s2, sa, sb;
	cpx_init (s0);
	cpx_init (s1);
	cpx_init (s2);
	cpx_init (sa);
	cpx_init (sb);

	cpx_t na, nb;
	cpx_init (na);
	cpx_init (nb);

	cpx_t y0, y1, y2;
	cpx_init (y0);
	cpx_init (y1);
	cpx_init (y2);

	mpf_t f0, f1, f2;
	mpf_init (f0);
	mpf_init (f1);
	mpf_init (f2);

	mpf_t loc, lam0, lam1, lam2;
	mpf_init (loc);
	mpf_init (lam0);
	mpf_init (lam1);
	mpf_init (lam2);

	/* Initial guess */
	cpx_set_d (s0, 0.5, 18.3);

	/* Initial directions */
	cpx_set_d (na, 0.05, 0.0);
	cpx_set_d (nb, 0.0, 0.1);

	int i;
	for (i=0; i<20; i++)
	{
		/* Three colinear points to start with */
		cpx_add(s1, s0, na);
		cpx_add(s2, s1, na);
		
		integral (y0, nsteps, s0, prec);
		integral (y1, nsteps, sa, prec);
		integral (y2, nsteps, s2, prec);

		cpx_abs(f0, y0);
		cpx_abs(f1, y1);
		cpx_abs(f2, y2);

		mpf_set_ui (lam0, 0);
		mpf_set_ui (lam1, 1);
		mpf_set_ui (lam2, 2);

		/* loc provides new minimum, along direction a */
		quad_min (loc, lam0, lam1, lam2, f0, f1, f2);

		/* Move to that location */
		mpf_abs (zero, loc);
		if (0 < mpf_cmp(zero, epsi))
		{
			cpx_times_mpf (na, na, loc);
			cpx_add (sa, s0, na);
		}

		/* Repeat for direction b */
		cpx_add(s1, sa, nb);
		cpx_add(s2, s1, nb);
		
		integral (y0, nsteps, sa, prec);
		integral (y1, nsteps, s1, prec);
		integral (y2, nsteps, s2, prec);

		cpx_abs(f0, y0);
		cpx_abs(f1, y1);
		cpx_abs(f2, y2);

		mpf_set_ui (lam0, 0);
		mpf_set_ui (lam1, 1);
		mpf_set_ui (lam2, 2);

		/* loc provides new minimum, along direction a */
		quad_min (loc, lam0, lam1, lam2, f0, f1, f2);

		/* Move to that location */
		/* Move to that location */
		mpf_abs (zero, loc);
		if (0 < mpf_cmp(zero, epsi))
		{
			cpx_times_mpf (nb, nb, loc);
			cpx_add (sb, sa, nb);
		}

		/* Shuffle down */
		// cpx_set(na, nb);
		// cpx_sub(nb, sb, s0);
		cpx_set(s0, sb);

printf("\n%d\n", i);
cpx_prt("s0 = ", s0); printf("\n");
cpx_prt("na = ", na); printf("\n");
cpx_prt("nb = ", nb); printf("\n");
	}

	cpx_clear (s0);
	cpx_clear (s1);
	cpx_clear (s2);
	cpx_clear (sa);
	cpx_clear (sb);

	cpx_clear (na);
	cpx_clear (nb);

	cpx_clear (y0);
	cpx_clear (y1);
	cpx_clear (y2);

	mpf_clear (f0);
	mpf_clear (f1);
	mpf_clear (f2);

	mpf_clear (loc);
	mpf_clear (lam0);
	mpf_clear (lam1);
	mpf_clear (lam2);

	mpf_clear (zero);
	mpf_clear (epsi);
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
