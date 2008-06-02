
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
	int prec = 60;

	/* Set the precision (number of binary bits) */
	int nbits = 3.3*prec+100;
	mpf_set_default_prec (nbits);

	mpf_t lg2;
	mpf_init (lg2);  // log 2
	fp_log2 (lg2, prec);

	mpf_t k, acc, prob, term, p_k;
	mpf_init (prob);
	mpf_init (acc);
	mpf_init (p_k);
	mpf_init (term);
	mpf_init (k);


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

		if (pcnt != 0)
		{
			fp_prt("k = ", k);
			mpf_div(term, prob, lg2);
			mpf_ui_sub (term, 1, term);
			fp_prt("prob = ", term);
			printf ("\n");
		}
		pcnt ++;

		mpf_add_ui(k, k, 1);
	}
}
