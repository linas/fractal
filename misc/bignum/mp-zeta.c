/*
 * mp_zeta.c
 *
 * High-precison Riemann zeta function, using the 
 * Gnu Multiple-precision library.
 * 
 * Linas Vepstas July 2005
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>

void fp_prt (char * str, mpf_t val)
{
	printf (str);
	mpf_out_str (stdout, 10, 60, val);
	printf ("\n");
}

/* ============================================================================= */
/* i_poch_rising
 * rising pochhammer symbol, for integer values.
 *
 * Brute force, simple.
 */

void i_poch_rising (mpz_t poch, unsigned int k, unsigned int n)
{
	mpz_t acc;
	
	mpz_init (acc);

	mpz_set_ui (poch, 1);
	unsigned int i;
	for (i=0; i<n; i++)
	{
		mpz_mul_ui (acc, poch, i+k);
		mpz_set (poch, acc);
	}

	mpz_clear (acc);
}

/* i_factorial -- the factorial
 */
void i_factorial (mpz_t fact, unsigned int n)
{
	i_poch_rising (fact, 1, n);
}

/* ============================================================================= */
/* fp_binomial
 * Binomial coefficient
 */

void i_binomial (mpz_t bin, unsigned int n, unsigned int k)
{
	mpz_t top, bot;

	if (2*k < n) k = n-k;

	mpz_init (top);
	mpz_init (bot);
	i_poch_rising (top, k+1, n-k);
	i_factorial (bot, n-k); 

	mpz_divexact (bin, top, bot);
	
	mpz_clear (top);
	mpz_clear (bot);
}

/* ============================================================================= */
/* fp_euler
 * return Euler-Mascheroni const
 */
void fp_euler_mascheroni (mpf_t gam)
{
	char * g="0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495e0";
	
	mpf_set_str (gam, g, 10);
}

void fp_pi (mpf_t pi)
{
	char *p = "3.1415926535 8979323846 2643383279 5028841971 6939937510"
			   "5820974944 5923078164 0628620899 8628034825 3421170679"
				 "8214808651 3282306647 0938446095 5058223172 5359408128"
				  "4811174502 8410270193 8521105559 6446229489 5493038196"
				   "4428810975 6659334461 2847564823 3786783165 2712019091"
					 "4564856692 3460348610 4543266482 1339360726 0249141273"
					  "7245870066 0631558817 4881520920 9628292540 9171536436"
					   "7892590360 0113305305 4882046652 1384146951 9415116094"
						 "3305727036 5759591953 0921861173 8193261179 3105118548"
						  "0744623799 6274956735 1885752724 8912279381 8301194912"
						   "9833673362 4406566430 8602139494 6395224737 1907021798"
							 "6094370277 0539217176 2931767523 8467481846 7669405132"
							  "0005681271 4526356082 7785771342 7577896091 7363717872"
							   "1468440901 2249534301 4654958537 1050792279 6892589235"
								 "4201995611 2129021960 8640344181 5981362977 4771309960"
								  "5187072113 4999999837 2978049951 0597317328 1609631859"
								   "5024459455 3469083026 4252230825 3344685035 2619311881";
	mpf_set_str (pi, p, 10);
}

void fp_zeta2 (mpf_t zeta)
{
	mpf_t pi, pisq;
	mpf_init (pi);
	mpf_init (pisq);
	
	fp_pi (pi);
	mpf_mul (pisq, pi, pi);
	mpf_div_ui (zeta, pisq, 6);

	mpf_clear (pi);
	mpf_clear (pisq);
}

void fp_zeta3 (mpf_t zeta)
{
	// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
	// char * g="1.202056903159594285399738161511449990764986292";
	char * g="1.2020569031595942853997381615114499907649862923404988817922715553418382057863 130901864558736093352581461991577952607194184919959986732832137763968372079001 614539417829493600667191915755222424942439615639096641032911590957809655146512 799184051057152559880154371097811020398275325667876035223369849416618110570147 157786394997375237852779370309560257018531827900030765471075630488433208697115";

	mpf_set_str (zeta, g, 10);
}

void fp_zeta5 (mpf_t zeta)
{
	// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
	char * g="1.036927755143369926331365486457034168057080919501912811974192677 9038035897862814845600431065571333363796203414665566090428009617 7915597084183511072180087644866286337180353598363962365128888981 3352767752398275032022436845766444665958115993917977745039244643 9196666159664016205325205021519226713512567859748692860197447984 3200672681297530919900774656558601526573730037561532683149897971 9350398378581319922884886425335104251602510849904346402941172432 7576341508162332245618649927144272264614113007580868316916497918";

	mpf_set_str (zeta, g, 10);
}

void fp_zeta7 (mpf_t zeta)
{
	// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
	char * g="1.008349277381922826839797549849796759599863560565238706417283136 5716014783173557353460969689138513239689614536514910748872867774 1984033544031579830103398456212106946358524390658335396467699756 7696691427804314333947495215378902800259045551979353108370084210 7329399046107085641235605890622599776098694754076320000481632951 2586769250630734413632555601360305007373302413187037951026624779 3954650225467042015510405582224239250510868837727077426002177100 0195455778989836046745406121952650765461161356548679150080858554";

	mpf_set_str (zeta, g, 10);
}

void fp_zeta9 (mpf_t zeta)
{
	// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
	char * g="1.002008392826082214417852769232412060485605851394888756548596615 9097850533902583989503930691271695861574086047658470602614253739 7072243015306913249876425109092948687676545396979415407826022964 1544836250668629056707364521601531424421326337598815558052591454 0848901539527747456133451028740613274660692763390016294270864220 1123162209241265753326205462293215454665179945038662778223564776 1660330281492364570399301119383985017167926002064923069795850945 8457966548540026945118759481561430375776154443343398399851419383";

	mpf_set_str (zeta, g, 10);
}

void fp_zeta4 (mpf_t zeta)
{
	mpf_t pi, pisq;
	mpf_init (pi);
	mpf_init (pisq);
	
	fp_pi (pi);
	mpf_mul (pisq, pi, pi);
	mpf_mul (pi, pisq, pisq);
	mpf_div_ui (zeta, pi, 90);

	mpf_clear (pi);
	mpf_clear (pisq);
}

void fp_zeta6 (mpf_t zeta)
{
	mpf_t pi, pisq, ph;
	mpf_init (pi);
	mpf_init (pisq);
	mpf_init (ph);
	
	fp_pi (pi);
	mpf_mul (pisq, pi, pi);
	mpf_mul (pi, pisq, pisq);
	mpf_mul (ph, pi, pisq);
	mpf_div_ui (zeta, ph, 945);

	mpf_clear (pi);
	mpf_clear (pisq);
	mpf_clear (ph);
}

void fp_zeta8 (mpf_t zeta)
{
	mpf_t pi, pisq;
	mpf_init (pi);
	mpf_init (pisq);
	
	fp_pi (pi);
	mpf_mul (pisq, pi, pi);
	mpf_mul (pi, pisq, pisq);
	mpf_mul (pisq, pi, pi);
	mpf_div_ui (zeta, pisq, 9450);

	mpf_clear (pi);
	mpf_clear (pisq);
}
/* ============================================================================= */
/* fp_zeta
 * Floating-point-valued Riemann zeta for positive integer arguments 
 * return value placed in the arg "zeta".
 *
 * Simple-minded algo, carries out math to prec decimal digits
 */
void fp_zeta (mpf_t zeta, unsigned int s, int prec)
{
	unsigned long int us = s;

	if (s<2) return;
	switch (s)
	{
		case 2: fp_zeta2 (zeta); return;
		case 3: fp_zeta3 (zeta); return;
		case 4: fp_zeta4 (zeta); return;
		case 5: fp_zeta5 (zeta); return;
		case 6: fp_zeta6 (zeta); return;
		case 7: fp_zeta7 (zeta); return;
		case 8: fp_zeta8 (zeta); return;
		case 9: fp_zeta9 (zeta); return;
	}
	
	mpf_t acc;
	mpf_t term;
	mpf_t en;
	mpf_t inv;
	
	mpf_init (acc);
	mpf_init (term);
	mpf_init (en);
	mpf_init (inv);
	
	mpf_set_ui (zeta, 1);

	/* Compute number of terms to be carried out 
	 * However, this estimate is wrong; it stops 
	 * when next term is "smaller than" rather 
	 * than when its converged.
	 */
	double fprec = prec;
	fprec /= (double) s;
	double dig = pow (10.0, fprec);
	if (1.0e10 < dig)
	{
		fprintf (stderr, "Sorry bucko, can't do it\n");
		return;
	}
	int nmax = dig+1.0;
	// printf ("zeta will be computed with %d terms\n", nmax);
	
	int n;
	for (n=2; n<= nmax; n++)
	{
		mpf_set_ui (en, n);
		mpf_ui_div (inv, 1, en);  /* inv = 1/n */
		mpf_pow_ui (term, inv, us); /* term = 1/n^s */
		mpf_add (acc, zeta, term);
		mpf_set (zeta, acc);
	}

	mpf_clear (acc);
	mpf_clear (term);
	mpf_clear (en);
	mpf_clear (inv);
}

/* ============================================================================= */
/* compute a_sub_n
 */
void a_sub_n (mpf_t a_n, unsigned int n, unsigned int prec)
{
	int k;
	mpf_t fbin, term, zt, ok, one, acc, zeta;
	mpf_t gam;

	mpf_init (term);
	mpf_init (acc);
	mpf_init (zeta);
	mpf_init (zt);
	mpf_init (ok);
	mpf_init (one);
	mpf_init (fbin);
	mpf_init (gam);
	
	mpf_set_ui (one, 1);

	mpz_t ibin;
	mpz_init (ibin);
	mpf_set_ui (a_n, 0);

	for (k=1; k<= n; k++)
	{
		fp_zeta (zeta, k+1, prec);
		mpf_div_ui (zt, zeta, k+1);
		mpf_div_ui (ok, one, k);
		mpf_sub (term, ok, zt);

		i_binomial (ibin, n, k);
		mpf_set_z (fbin, ibin);

		mpf_mul (zeta, term, fbin);

		if (k%2) mpf_neg (term, zeta);
		else mpf_set (term, zeta);
		
		mpf_add (acc, a_n, term);
		mpf_set (a_n, acc);
	}

	/* add const terms */
	mpf_add_ui (term, a_n, 1);
	fp_euler_mascheroni (gam);
	mpf_sub (a_n, term, gam);

	/* subtract 1/2(n+1) */
	mpf_div_ui (ok, one, 2*(n+1));
	mpf_sub (term, a_n, ok);
	mpf_set (a_n, term);
	
	mpf_clear (term);
	mpf_clear (acc);
	mpf_clear (zeta);
	mpf_clear (zt);
	mpf_clear (ok);
	mpf_clear (one);
	mpf_clear (fbin);
	mpf_clear (gam);

	mpz_clear (ibin);
}

void a_bound_n (mpf_t b_n, unsigned int n)
{
	mpf_t en, sq_en;
	mpf_init (en);
	mpf_init (sq_en);

	mpf_set_ui (en, n+1);
	mpf_sqrt (sq_en, en);
	mpf_neg (en, sq_en);

	mpf_clear (en);
	mpf_clear (sq_en);
}

/* ============================================================================= */
main ()
{
	char str[4000];

#ifdef FACT_TEST
	mpz_t fact;
	mpz_init (fact);

	i_factorial (fact, 5);
	mpz_get_str (str, 10, fact);
	printf ("fact = %s\n", str);
#endif

#ifdef BINOMIAL_TEST
	int n, k;
	mpz_t bin;
	mpz_init (bin);

	for (n=1; n<7; n++)
	{
		for (k=0; k<=n; k++)
		{
			i_binomial (bin, n ,k);
			mpz_get_str (str, 10, bin);
			printf ("bin (%d %d) = %s\n", n, k, str);
		}
		printf ("---\n");
	}
#endif
	
	/* set the precision */
	mpf_set_default_prec (800);
	
#ifdef ZETA_STUFF
	mpf_t zeta;
	mpf_init (zeta);
	
	printf ("           000000000011111111112222222222333333333344444444445555555555666666666677777777778\n");
	printf ("           012345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
	fp_zeta (zeta, 2, 13);
	fp_prt ("13 digs= ", zeta);
	fp_zeta (zeta, 8, 30);
	fp_prt ("30 digs= ", zeta);
	fp_zeta (zeta, 8, 40);
	fp_prt ("40 digs= ", zeta);
	fp_zeta (zeta, 8, 50);
	fp_prt ("50 digs= ", zeta);
	fp_zeta (zeta, 8, 60);
	fp_prt ("60 digs= ", zeta);
	fp_zeta (zeta, 8, 70);
	fp_prt ("70 digs= ", zeta);
	fp_zeta (zeta, 8, 80);
	fp_prt ("0 digs= ", zeta);
#endif
	
	mpf_t a_n, b_n, prod;
	mpf_init (a_n);
	mpf_init (b_n);
	mpf_init (prod);

	fp_pi (a_n);
	fp_prt ("duude pi ", a_n);

	int prec = 35;
	int n;
	for (n=0; n<150; n++)
	{
		a_sub_n (a_n, n, prec+n);

		double dbn = 1.0/exp (-4.0*sqrt (n+1));
		mpf_set_d (b_n, dbn);
		mpf_mul(prod, a_n, b_n);
		
		printf ("a(%d) ",n);
		fp_prt ("= ", prod);
		fflush (stdout);
	}

}

