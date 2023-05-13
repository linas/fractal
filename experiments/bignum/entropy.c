
/* 
 * entropy.c
 *
 * Compute the Gauss-Kuzmin entropy
 * sum_k p_k log(p_k)
 * where
 * p_k = log_2 (1-1/(k+1)^2)
 *
 * Linas Vepstas June 2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mp-consts.h"
#include "mp-misc.h"
#include "mp-trig.h"

int
main (int argc, char * argv[])
{
	int prec = 40;

	/* Set the precision (number of binary bits) */
	int nbits = 3.3*prec+30;
	mpf_set_default_prec (nbits);

	printf ("percision(base-10)=%d nbits=%d\n", prec, nbits);

	// log2 and log log 2
	mpf_t lg2, lglg2;
	mpf_init (lg2);  // log 2
	mpf_init (lglg2);  // log log 2
	fp_log2 (lg2, prec);
	fp_log(lglg2, lg2, prec);

	// terms in the series.
	mpf_t k, acc, prob, term, p_k;
	mpf_init (k);
	mpf_init (p_k);
	mpf_init (acc);
	mpf_init (term);
	mpf_init (prob);

	// Aitken delta-squared method for series acceleration.
	mpf_t s_n, s_nm1, s_nm2, delta, dprev, tmp;
	mpf_init (s_n);
	mpf_init (s_nm1);
	mpf_init (s_nm2);
	mpf_init (delta);
	mpf_init (dprev);
	mpf_init (tmp);

	// quadratic extrapolation
	mpf_t tee, ent, x_0, x_1, x_2, y_0, y_1, y_2, d01, d02, d12;
	mpf_init(tee);
	mpf_init(ent);
	mpf_init(x_0);
	mpf_init(x_1);
	mpf_init(x_2);
	mpf_init(y_0);
	mpf_init(y_1);
	mpf_init(y_2);
	mpf_init(d01);
	mpf_init(d02);
	mpf_init(d12);
	mpf_set_ui(x_1, 1);
	mpf_set_ui(x_2, 2);

	unsigned long int pcnt = 1;
	
	mpf_set_ui(k, 2);
	while(1)
	{
		// Compute p_k = - log (1-1/(k+1)^2)
		mpf_mul(p_k, k, k);
		mpf_ui_div(p_k, 1, p_k);
		fp_log_m1 (p_k, p_k, prec);

		// Compute sum_k p_k
		mpf_add (prob, prob, p_k);

		// Compute entropy
		fp_log(term, p_k, prec);
		mpf_mul(term, p_k, term);
		mpf_add (acc, acc, term);

		if (pcnt%640123 == 0)
		{
			fp_prt("k = ", k);
			printf ("\n");

			// total probability is tending to 1.0,
			// so print difference from 1.0
			mpf_div(term, prob, lg2);
			mpf_ui_sub (tee, 1, term);
			fp_prt("prob = ", tee);
			printf ("\n");

			// print the regular entropy
			mpf_div(term, acc, lg2);
			mpf_sub(ent, lglg2, term);
			fp_prt("H = ", ent);
			printf ("\n");

			// Aitken method
			// s = s_nm2 + (s_nm2-s_nm1)^2 / (s_n + s_nm2 - 2*s_nm1)
			// (sucks, for this particular problem)
			mpf_set (s_nm2, s_nm1);
			mpf_set (s_nm1, s_n);
			mpf_set (s_n, ent);
			mpf_set (dprev, delta);
			mpf_sub (delta, s_n, s_nm1);

			mpf_mul(term, dprev, dprev);
			mpf_sub(tmp, delta, dprev);
			mpf_div(term, term, tmp);
			fp_prt("delta = ", term);
			printf ("\n");
			mpf_sub(term, s_n, term);
			fp_prt("H aitken = ", term);
			printf ("\n");

			// print the information-theoretic entropy
			// which differs by a log2.
			mpf_div(term, term, lg2);
			fp_prt("H_2 (aitken) = ", term);
			printf ("\n");

			// quadratic exprapolation
			// quadratic extrapolation as a function of probability 
			// extrapolated to probability 1 should work very very well.
			// The graph appears to be very nearly linear. I guess its
			// xlogx in shape...
			mpf_set(x_2, x_1);
			mpf_set(x_1, x_0);
			mpf_set(x_0, tee);

			mpf_set(y_2, y_1);
			mpf_set(y_1, y_0);
			mpf_set(y_0, ent);

			mpf_sub(d01, x_0, x_1);
			mpf_sub(d02, x_0, x_2);
			mpf_sub(d12, x_1, x_2);

			mpf_mul(term, x_1, x_2);
			mpf_mul(term, term, y_0);
			mpf_div(term, term, d01);
			mpf_div(term, term, d02);

			mpf_mul(tmp, x_0, x_2);
			mpf_mul(tmp, tmp, y_1);
			mpf_div(tmp, tmp, d01);
			mpf_div(tmp, tmp, d12);
			mpf_sub(term, term, tmp);

			mpf_mul(tmp, x_0, x_1);
			mpf_mul(tmp, tmp, y_2);
			mpf_div(tmp, tmp, d02);
			mpf_div(tmp, tmp, d12);
			mpf_add(term, term, tmp);

			fp_prt("H quad= ", term);
			printf ("\n");

			// print the information-theoretic entropy
			// which differs by a log2.
			mpf_div(term, term, lg2);
			fp_prt("H_2 (quad) = ", term);
			printf ("\n");

			printf ("\n");
			fflush (stdout);
		}
		pcnt ++;

		mpf_add_ui(k, k, 1);
	}
}
