/*
 * exact-bignum.c
 *
 * Implementation of exact-series. but in GMP
 *
 * Linas Vepstas December 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-binomial.h>
#include <mp-complex.h>

/* Integral of (x^(s-1))/(x+alpha) dx
 */
void eff (cpx_t sum, cpx_t s, mpf_t alpha, mpf_t x, int nprec)
{
	int k;
	mpf_t zero, epsi;
	mpf_init (zero);
	mpf_init (epsi);

	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.321*nprec));

	cpx_t term, ess;
	cpx_init (term);
	cpx_init (ess);

	mpf_t opxa, opxak, tmp;
	mpf_init (opxa);
	mpf_init (opxak);
	mpf_init (tmp);

	cpx_set(ess, s);
	cpx_set_ui(sum, 0, 0);

	mpf_div (opxa, x, alpha);
	mpf_add_ui (opxa, opxa, 1);
	mpf_set (opxak, opxa);

	for (k=1; k<155456123; k++)
	{
		cpx_binomial (term, ess, k);
		cpx_times_mpf (term, term, opxak);
		cpx_div_ui (term, term, k);
		if (k%2)
		{
			cpx_sub(sum, sum, term);
		}
		else
		{
			cpx_add (sum, sum, term);
		}

		/* Check to see if the last term is small enough */
		cpx_abs(zero, term);
		if (0 > mpf_cmp(zero, epsi)) break;

		mpf_mul (opxak, opxak, opxa);
	}

	/* Add in the logarithmic part */
	mpf_add(tmp, x, alpha);
	if (0< mpf_sgn(tmp))
	{
		fp_log(tmp, tmp, nprec);
		mpf_add(sum[0].re, sum[0].re, tmp);
	}
	else
	{
		mpf_neg(tmp, tmp);
		fp_log(tmp, tmp, nprec);
		mpf_sub(sum[0].re, sum[0].re, tmp);
	}

	/* scale by (-alpha)^s */
	mpf_neg(tmp, alpha);
	cpx_mpf_pow(term, tmp, ess, nprec);
	cpx_mul(sum, sum, term);

	cpx_clear(term);
	cpx_clear (ess);
	mpf_clear (opxa);
	mpf_clear (opxak);
	mpf_clear (epsi);
	mpf_clear (zero);
	mpf_clear (tmp);
}

/* Integral of (x^s) dx = (x^(s+1))/(s+1) */
void effo (cpx_t f, cpx_t s, mpf_t x, int prec)
{
	cpx_t essp1;
	cpx_init(essp1);
	cpx_add_ui(essp1, s, 1, 0);

	cpx_mpf_pow(f, x, essp1, prec);
	cpx_div (f, f, essp1);

	cpx_clear(essp1);
}

/* Integral of (x^(s-1))/(x+alpha) dx
 * however, converges differently
 */
void special_eff (cpx_t f, cpx_t s, mpf_t alpha, mpf_t x, int prec)
{
	exit(1);  // not yet implemented
}


void gral_s12 (cpx_t sum, cpx_t s, int ndigits, int a1max, int prec)
{
	unsigned int na1, na2;
	mpf_t a1, a2, a,b,c,d, rat, xlo, xhi;
	mpf_init(a1);
	mpf_init(a2);
	mpf_init(a);
	mpf_init(b);
	mpf_init(c);
	mpf_init(d);
	mpf_init(rat);
	mpf_init(xlo);
	mpf_init(xhi);

	cpx_t ess, essm1, thi, tlo, term;
	cpx_init(ess);
	cpx_init(essm1);
	cpx_init(thi);
	cpx_init(tlo);
	cpx_init(term);

	cpx_set (ess, s);
	cpx_sub_ui (essm1, s, 1, 0);
	cpx_set_ui (sum, 0, 0);

	for (na1=1; na1<a1max; na1++)
	{
		mpf_set_ui (a1, na1);
		for (na2=1; 1; na2++)
		{
			mpf_set_ui (a2, na2);
			mpf_sub(b, a1, a2);    // b = a1 - a2
			mpf_mul(c, a1, a2); 
			mpf_add_ui(c, c, 1);
			mpf_div(xlo, a2, c);   // xlo = a2 / (1 + a1 * a2)
			mpf_add_ui(xhi, a2, 1);
			mpf_add(rat, a1, c);
			mpf_div(xhi, xhi, rat); // xhi = (1 + a2) / (1 + a1 + a1 * a2)
			mpf_mul(c, c, b);
			mpf_neg(c, c);         // c = -(1 + a1 * a2) * b
			mpf_mul(d, a2, b);
			mpf_add_ui(d, d, 1);   // d = 1+a2*b
			mpf_mul(a, a1, b);
			mpf_ui_sub(a, 1, a);   // a = 1-a1*b
			mpf_div(rat, d, c);    // rat = d/c

			if ((1 == na1) && (2 == na2))
			{
			}
			else if (na1 == na2)
			{
				/* special case, where c=d=0, so general formula does not apply */
				effo(term, s, xhi, prec);
				cpx_add (sum, sum, term);

				effo(term, s, xlo, prec);
				cpx_sub (sum, sum, term);
			}
			else
			{
				eff(thi, ess, rat, xhi, prec);
				eff(tlo, ess, rat, xlo, prec);
				cpx_sub(thi, thi, tlo);
				cpx_mul(term, a, thi);
			}
		}
	}

	mpf_clear(a1);
	mpf_clear(a2);
	mpf_clear(a);
	mpf_clear(b);
	mpf_clear(c);
	mpf_clear(d);
	mpf_clear(rat);
	mpf_clear(xlo);
	mpf_clear(xhi);

	cpx_clear(ess);
	cpx_clear(essm1);
	cpx_clear(thi);
	cpx_clear(tlo);
	cpx_clear(term);
}

int main (int argc, char * argv[])
{
}
