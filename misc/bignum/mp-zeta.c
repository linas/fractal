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
#include <stdlib.h>

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
/* floating point exponential */

void fp_exp (mpf_t ex, mpf_t z, unsigned int prec)
{
	mpf_t z_n, fact, term, acc, tmp;

	mpf_init (z_n);
	mpf_init (fact);
	mpf_init (term);
	mpf_init (acc);
	mpf_init (tmp);

	mpf_set_ui (ex, 1);
	mpf_set_ui (fact, 1);
	mpf_set (z_n, z);
	
	double mex = ((double) prec) * log (10.0) / log(2.0);
	unsigned int imax = mex +1.0;
	mpf_t maxterm, one;
	mpf_init (maxterm);
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_div_2exp (maxterm, one, imax);

	unsigned int n=1;
	while(1)
	{
		mpf_div (term, z_n, fact);
		mpf_add (acc, ex, term);
		mpf_set (ex,acc);
		
		/* don't go no father than this */
		mpf_abs (tmp, term);
		if (mpf_cmp (tmp, maxterm) < 0) break;
		
		n++;
		mpf_mul (tmp, z_n, z);
		mpf_set (z_n, tmp);
		mpf_mul_ui (tmp, fact, n);
		mpf_set (fact,tmp);
	}
	
	mpf_clear (z_n);
	mpf_clear (fact);
	mpf_clear (term);
	mpf_clear (acc);
	mpf_clear (tmp);

	mpf_clear (one);
	mpf_clear (maxterm);
}

/* ============================================================================= */
/* fp_euler
 * return Euler-Mascheroni const
 */
void fp_euler_mascheroni (mpf_t gam)
{
	// char * g="0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495e0";
	char * g="0.57721566490153286060651209008240243104215933593992359880576723488486772677766 467093694706329174674951463144724980708248096050401448654283622417399764492353 625350033374293733773767394279259525824709491600873520394816567085323315177661 152862119950150798479374508570574002992135478614669402960432542151905877553526 733139925401296742051375413954911168510280798423487758720503843109399736137255 306088933126760017247953783675927135157722610273492913940798430103417771778088 154957066107501016191663340152278935867965497252036212879226555953669628176388";
	
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

/* return e^pi */
void fp_e_pi (mpf_t e_pi)
{
	char *p= "23.140692632779269005729086367948547380266106242600211993445046409524342350690 452783516971997067549219675952704801087773144428044414693835844717445879609849 365327965863669242230268991013741764684401410395183868477243068059588162449844 491430966778413671631963414784038216511287637731470347353833162821294047891936 224820221006032065443362736557271823744989618858059591684872645479013397834026 595101499643792422968160799565381423536206957600770590460899883002254304871211 791300849327379580729427301931042601691939325853203428968661895283290521711157 185185506802254197204566370865568386830544799278170407497768540367556534957218 867882563994384718224585889428535247260568210271076018491534518468064887386774 439630514005169440540665265430968869063937315359837311042174433023967896690035";
	mpf_set_str (e_pi, p, 10);
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

void fp_zeta_even (mpf_t zeta, unsigned int n, unsigned int div)
{
	mpf_t pi, pip;
	mpf_init (pi);
	mpf_init (pip);
	
	fp_pi (pi);
	mpf_pow_ui (pip, pi, n);
	mpf_div_ui (zeta, pip, div);

	mpf_clear (pi);
	mpf_clear (pip);
}

void fp_zeta_even_str (mpf_t zeta, unsigned int n, char * snum, char * sdenom)
{
	mpf_t pi, pip, num, denom;
	mpf_init (pi);
	mpf_init (pip);
	mpf_init (num);
	mpf_init (denom);

	mpf_set_str (num,snum, 10);
	mpf_set_str (denom, sdenom, 10);
	
	fp_pi (pi);
	mpf_pow_ui (pip, pi, n);
	mpf_mul(pi, pip, num);
	mpf_div(zeta, pi, denom);

	mpf_clear (pi);
	mpf_clear (pip);
	mpf_clear (num);
	mpf_clear (denom);
}

/* return sum_n (n^k (e^{\pi k} \pm 1)^{-1}
 * The Simon Plouffe Ramanujan inspired thingy
 */
void fp_ess (mpf_t ess_plus, mpf_t ess_minus, unsigned int k, unsigned int prec)
{
	mpf_t e_pi, en, enp, epip, eppos, epneg, term, oterm, acc;

	mpf_init (e_pi);
	mpf_init (en);
	mpf_init (enp);
	mpf_init (epip);
	mpf_init (eppos);
	mpf_init (epneg);
	mpf_init (term);
	mpf_init (oterm);
	mpf_init (acc);

	fp_e_pi (e_pi);
	mpf_set_ui (ess_plus, 0);
	mpf_set_ui (ess_minus, 0);

	double mex = ((double) prec) * log (10.0) / log(2.0);
	unsigned int imax = mex +1.0;
	mpf_t maxterm, one;
	mpf_init (maxterm);
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_mul_2exp (maxterm, one, imax);
	
	int n;
	for (n=1; n<1000000000; n++)
	{
		mpf_set_ui (en, n);
		mpf_pow_ui (enp, en, k);
		mpf_pow_ui (epip, e_pi, 2*n);
		
		mpf_add_ui (eppos, epip, 1);
		mpf_sub_ui (epneg, epip, 1);

		mpf_mul (term, enp, eppos);
		mpf_ui_div (oterm, 1, term);
		mpf_add (acc, ess_plus, oterm);
		mpf_set (ess_plus, acc);

		mpf_mul (term, enp, epneg);
		mpf_ui_div (oterm, 1, term);
		mpf_add (acc, ess_minus, oterm);
		mpf_set (ess_minus, acc);

		/* don't go no father than this */
		if (mpf_cmp (term, maxterm) > 0) break;
	}

	mpf_clear (e_pi);
	mpf_clear (en);
	mpf_clear (enp);
	mpf_clear (epip);
	mpf_clear (eppos);
	mpf_clear (epneg);
	mpf_clear (term);
	mpf_clear (oterm);
	mpf_clear (acc);

	mpf_clear (one);
	mpf_clear (maxterm);
}

/* Implement Simon Plouffe odd-zeta sums */
void fp_zeta_odd (mpf_t zeta, unsigned int n, 
					 char *sdiv, char * spi, char * sminus, char * splus,  
					 unsigned int prec)
{
	mpf_t pi, pip, piterm, spos, sneg, spos_term, sneg_term, tmp;
	mpf_init (pi);
	mpf_init (pip);
	mpf_init (piterm);
	mpf_init (spos);
	mpf_init (sneg);
	mpf_init (spos_term);
	mpf_init (sneg_term);
	mpf_init (tmp);

	mpf_t div, c_pi, c_plus, c_minus;
	mpf_init (div);
	mpf_init (c_pi);
	mpf_init (c_plus);
	mpf_init (c_minus);
	
	mpf_set_str (div, sdiv, 10);
	mpf_set_str (c_pi, spi, 10);
	mpf_set_str (c_plus, splus, 10);
	mpf_set_str (c_minus, sminus, 10);
	
	fp_ess (spos, sneg, n, prec);
	mpf_mul (spos_term, spos, c_plus);
	mpf_mul (sneg_term, sneg, c_minus);
			  
	fp_pi (pi);
	mpf_pow_ui (pip, pi, n);
	mpf_mul (piterm, pip, c_pi);

	mpf_set (tmp, piterm);
	mpf_sub (zeta, tmp, spos_term);
	mpf_sub (tmp, zeta, sneg_term);
	mpf_div (zeta, tmp, div);
	
	mpf_clear (pi);
	mpf_clear (pip);
	mpf_clear (piterm);
	mpf_clear (spos);
	mpf_clear (sneg);
	mpf_clear (spos_term);
	mpf_clear (sneg_term);
	mpf_clear (tmp);
	
	mpf_clear (div);
	mpf_clear (c_pi);
	mpf_clear (c_plus);
	mpf_clear (c_minus);
}

/* ============================================================================= */
/* fp_zeta
 * Floating-point-valued Riemann zeta for positive integer arguments 
 * return value placed in the arg "zeta".
 *
 * Simple-minded algo, carries out math to prec decimal digits
 *
 * OLEIS A046988
 *  Numerator of zeta(2n)/Pi^(2n)
 * 0,1,1,1,1,1,691,2,3617,43867,174611,155366,236364091,
 * 1315862,6785560294,6892673020804,7709321041217,151628697551,
 * 26315271553053477373,308420411983322,261082718496449122051,
 * 3040195287836141605382,5060594468963822588186
 *
 * OLEIS A002432 
 *  Denominator of zeta(2n)/Pi^(2n)
 *  6,90,945,9450,93555,638512875,18243225,325641566250,
 *  38979295480125,1531329465290625,13447856940643125,
 *  201919571963756521875,11094481976030578125,
 *  564653660170076273671875,5660878804669082674070015625,
 *  62490220571022341207266406250,12130454581433748587292890625
 */
void fp_zeta (mpf_t zeta, unsigned int s, int prec)
{
	unsigned long int us = s;

#define ARRSZ 4000
	static int first_time = 1;
	static mpf_t zeta_cache[ARRSZ];
	static int zprec[ARRSZ];
	static unsigned int last_term[ARRSZ];

	/* Cache of best values so far */
	if(first_time)
	{
		first_time = 0;
		int i;
		for (i=0; i<ARRSZ; i++)
		{
			zprec[i] = 0;
		}
	}

	if (s<2) return;
	switch (s)
	{
		case 2: fp_zeta2 (zeta); return;
		case 3: fp_zeta3 (zeta); return;
		case 4: fp_zeta_even (zeta, 4, 90); return;
		case 5: fp_zeta5 (zeta); return;
		case 6: fp_zeta_even (zeta, 6, 945); return;
		case 7: fp_zeta7 (zeta); return;
		case 8: fp_zeta_even (zeta, 8, 9450); return;
		case 9: fp_zeta9 (zeta); return;
		case 10: fp_zeta_even (zeta, 10, 93555); return;
					
		case 12: fp_zeta_even_str (zeta, 12, "691", "638512875"); return;
		case 14: fp_zeta_even_str (zeta, 14, "2", "18243225"); return;
		case 16: fp_zeta_even_str (zeta, 16, "3617", "325641566250"); return;
		case 18: fp_zeta_even_str (zeta, 18, "43867", "38979295480125"); return;
		case 20: fp_zeta_even_str (zeta, 20, "174611", "1531329465290625"); return;
		case 22: fp_zeta_even_str (zeta, 22, "155366", "13447856940643125"); return;
		case 24: fp_zeta_even_str (zeta, 24, "236364091", "201919571963756521875"); return;
		case 26: fp_zeta_even_str (zeta, 26, "1315862", "11094481976030578125"); return;
		case 28: fp_zeta_even_str (zeta, 28, "6785560294", "564653660170076273671875"); return;
		case 30: fp_zeta_even_str (zeta, 30, "6892673020804", "5660878804669082674070015625"); return;
		case 32: fp_zeta_even_str (zeta, 32, "7709321041217", "62490220571022341207266406250"); return;
		case 34: fp_zeta_even_str (zeta, 34, "151628697551", "12130454581433748587292890625"); return;
		case 11: 
			fp_zeta_odd (zeta, 11, "425675250", "1453", "851350500", "0", prec); 
			return;
		case 13: 
			fp_zeta_odd (zeta, 13, "257432175", "89", "514926720", "62370", prec); 
			return;
		case 15: 
			fp_zeta_odd (zeta, 15, "390769879500", "13687", "781539759000", "0", prec); 
			return;
		case 17: 
			fp_zeta_odd (zeta, 17, "1904417007743250", "6758333", "3808863131673600", "29116187100", prec); 
			return;
		case 19: 
			fp_zeta_odd (zeta, 19, "21438612514068750", "7708537", "42877225028137500", "0", prec); 
			return;
		case 21: 
			fp_zeta_odd (zeta, 21, "1881063815762259253125", "68529640373", "3762129424572110592000", "1793047592085750", prec); 
			return;
	}

	/* If we are here, well have to compute values using brute force */
	/* But first, lets see if we can get lucky with the cache. */
	if ((s < ARRSZ) && (prec < zprec[s]))
	{
		mpf_set (zeta, zeta_cache[s]);
		return;
	}

	/* initialize the cache line, if needed */
	if ((s < ARRSZ) && (0 == zprec[s]))
	{
		mpf_init (zeta_cache[s]);
		mpf_set_ui (zeta_cache[s], 1);
		last_term[s] = 2;
		zprec[s] = prec;
	}
	
	/* If we are here, well have to compute values using brute force */
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
	int nmax = 1.1*dig+1.0;
	// printf ("zeta will be computed with %d terms\n", nmax);
	
	/* Start computations where we last left off. */
	int nstart = 2;
	if (s < ARRSZ)
	{
		mpf_set (zeta, zeta_cache[s]);
		nstart = last_term[s];
	}
	
	int n;
	for (n=nstart; n< nmax; n++)
	{
		mpf_set_ui (en, n);
		mpf_ui_div (inv, 1, en);  /* inv = 1/n */
		mpf_pow_ui (term, inv, us); /* term = 1/n^s */
		mpf_add (acc, zeta, term);
		mpf_set (zeta, acc);
	}

	/* cache the results */
	if (s < ARRSZ)
	{
		mpf_set (zeta_cache[s], zeta);
		last_term[s] = nmax;
	}
	
	mpf_clear (acc);
	mpf_clear (term);
	mpf_clear (en);
	mpf_clear (inv);
}

/* ============================================================================= */
/* rough count of number of digits in a number */

static inline unsigned int num_digits (mpz_t num, mpz_t tmpa, mpz_t tmpb)
{
	unsigned int n=0;
	
	mpz_set (tmpb, num);
	while (1)
	{
		mpz_fdiv_q_ui (tmpa, tmpb, 100);
		mpz_set (tmpb, tmpa);
		if (0 == mpz_sgn  (tmpa)) break;
		n += 2;
	}
	return n;
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
	
	mpz_t tmpa, tmpb;
	mpz_init (tmpa);
	mpz_init (tmpb);
	
	mpf_set_ui (one, 1);

	mpz_t ibin;
	mpz_init (ibin);
	mpf_set_ui (a_n, 0);

	int maxbump = 0;
	for (k=1; k<= n; k++)
	{
		/* commputer the binomial */
		i_binomial (ibin, n, k);
		mpf_set_z (fbin, ibin);

		/* The terms will have laternating signes, and
		 * will mostly cancel one-another. Thus, we need 
		 * to increase precision for those terms with the 
		 * largest binomial coefficients. This is will
		 * increase precision for the killer terms, 
		 * while keeping the others in bearable range,
		 * in terms to cpu time consumed.
		 */
		int ndigits = num_digits (ibin, tmpa,tmpb);
		if (maxbump < ndigits) maxbump = ndigits;

		/* compute 1/k - zeta (k+1)/(k+1) */
		fp_zeta (zeta, k+1, prec+ndigits);
		mpf_div_ui (zt, zeta, k+1);
		mpf_div_ui (ok, one, k);
		mpf_sub (term, ok, zt);

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
	mpz_clear (tmpa);
	mpz_clear (tmpb);

	// printf ("# max precision bump=%d\n", maxbump);
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

main (int argc, char * argv[])
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
	
	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s [ndigits] [nterms]\n", argv[0]); 
		exit (1);
	}

	/* the decimal precison (number of decimal places) */
	int prec = atoi (argv[1]);

	/* number of an's to compute */
	int nterms = atoi (argv[2]);

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	int bits = v + 30 + nterms/3;
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
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
	
// #define TEST_ZETA
#ifdef TEST_ZETA
	mpf_t zeta;
	mpf_init (zeta);
	// fp_zeta_odd (zeta, 3, 180, 7, 360, 0, 60); 
	// fp_prt ("duude zeta3= ", zeta);
	int i;
	int pr = 205;
	for (i=23; i<100; i+=2 ) {
		fp_zeta (zeta, i, pr);
		printf ("char * zeta_%d_%d = \"", i, pr);
		mpf_out_str (stdout, 10, pr, zeta);
		printf ("\";\n");
		fflush (stdout);
	}
#endif

#ifdef TEST_EXP
	mpf_t ex, one;
	mpf_init (ex);
	mpf_init (one);
	mpf_set_ui(one, 1);
	fp_exp (ex, one, 50);
	fp_prt ("e= ", ex);
	mpf_clear (ex);
	mpf_clear(one);
#endif
	
#ifdef TEST_DIGIT_COUNT
	mpz_t ival, tmpa, tmpb;
	mpz_init (ival);
	mpz_init (tmpa);
	mpz_init (tmpb);
	mpz_set_ui (ival, 3000000);
	int nd = num_digits (ival, tmpa, tmpb);
	printf ("found %d digits\n", nd);
#endif
	
#define A_SUB_N
#ifdef A_SUB_N
	mpf_t a_n, en, sq, term, b_n, prod;
	mpf_init (a_n);
	mpf_init (en);
	mpf_init (sq);
	mpf_init (term);
	mpf_init (b_n);
	mpf_init (prod);

	int n;
	printf ("#\n# zeta expansion terms \n#\n");
	printf ("# computed with variable precision of %d decimal places\n", prec);
	printf ("# computed with %d bits of default mpf \n", bits);
	for (n=0; n<nterms; n++)
	{
		a_sub_n (a_n, n, prec);

		/* compute the bound */
		mpf_set_ui (en, n+1);
		mpf_sqrt (sq, en);
		mpf_mul_ui (term, sq, 4);
		mpf_neg (en, term);
		fp_exp (b_n, en, prec);
		mpf_div (prod, a_n, b_n);
		
#ifdef FLT_BND
		double dbn = 1.0/exp (-4.0*sqrt (n+1));
		mpf_set_d (b_n, dbn);
		mpf_mul(prod, a_n, b_n);
#endif
		
		
		printf ("%d\t",n);
		fp_prt ("", prod);
		fflush (stdout);
	}
#endif

}

