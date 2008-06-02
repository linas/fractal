
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


	unsigned short pcnt = 1;
	
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

		if (pcnt == 0)
		{
			fp_prt("k = ", k);
			printf ("\n");

			// total probability is tending to 1.0,
			// so print difference from 1.0
			mpf_div(term, prob, lg2);
			mpf_ui_sub (term, 1, term);
			fp_prt("prob = ", term);
			printf ("\n");

			// print the regular entropy
			mpf_div(term, acc, lg2);
			mpf_sub(term, lglg2, term);
			fp_prt("H = ", term);
			printf ("\n");

			// Aitken method
			// s = s_nm2 + (s_nm2-s_nm1)^2 / (s_n + s_nm2 - 2*s_nm1)
			mpf_set (s_nm2, s_nm1);
			mpf_set (s_nm1, s_n);
			mpf_set (s_n, term);
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

			printf ("\n");
			fflush (stdout);
		}
		pcnt ++;

		mpf_add_ui(k, k, 1);
	}
}
