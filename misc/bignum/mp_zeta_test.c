/*
 * mp_zeta_test.c
 *
 * Small test suite for 
 * High-precison Riemann zeta function, using the 
 * Gnu Multiple-precision library.
 *
 * Currently not automated, done purely by visual inspection
 *
 * Linas Vepstas July 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include "mp-consts.h"
#include "mp-complex.h"
#include "mp-misc.h"
#include "mp-trig.h"
#include "mp_zeta.h"

/* ==================================================================== */

void fp_pi_string (mpf_t pi)
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
void fp_e_pi_string (mpf_t e_pi)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		char *p= "23.140692632779269005729086367948547380266106242600211993445046409524342350690 452783516971997067549219675952704801087773144428044414693835844717445879609849 365327965863669242230268991013741764684401410395183868477243068059588162449844 491430966778413671631963414784038216511287637731470347353833162821294047891936 224820221006032065443362736557271823744989618858059591684872645479013397834026 595101499643792422968160799565381423536206957600770590460899883002254304871211 791300849327379580729427301931042601691939325853203428968661895283290521711157 185185506802254197204566370865568386830544799278170407497768540367556534957218 867882563994384718224585889428535247260568210271076018491534518468064887386774 439630514005169440540665265430968869063937315359837311042174433023967896690035";
		mpf_set_str (e, p, 10);
	}
	mpf_set (e_pi, e);
}

/* fp_euler
 * return Euler-Mascheroni const
 */
void fp_euler_mascheroni_string (mpf_t gam)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// char * g="0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495e0";
	char * g="0.57721566490153286060651209008240243104215933593992359880576723488486772677766 467093694706329174674951463144724980708248096050401448654283622417399764492353 625350033374293733773767394279259525824709491600873520394816567085323315177661 152862119950150798479374508570574002992135478614669402960432542151905877553526 733139925401296742051375413954911168510280798423487758720503843109399736137255 306088933126760017247953783675927135157722610273492913940798430103417771778088 154957066107501016191663340152278935867965497252036212879226555953669628176388";
	
		mpf_set_str (e, g, 10);
	}
	mpf_set (gam, e);
}


void fp_zeta2 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t z;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (z);

		mpf_t pi, pisq;
		mpf_init (pi);
		mpf_init (pisq);
		
		fp_pi (pi, 1000);  // XXX  bad hard-coded prcision value
		mpf_mul (pisq, pi, pi);
		mpf_div_ui (z, pisq, 6);
	
		mpf_clear (pi);
		mpf_clear (pisq);
	}
	mpf_set (zeta, z);
}

void fp_zeta3 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
		// char * g="1.202056903159594285399738161511449990764986292";
		char * g="1.2020569031595942853997381615114499907649862923404988817922715553418382057863 130901864558736093352581461991577952607194184919959986732832137763968372079001 614539417829493600667191915755222424942439615639096641032911590957809655146512 799184051057152559880154371097811020398275325667876035223369849416618110570147 157786394997375237852779370309560257018531827900030765471075630488433208697115";
		mpf_set_str (e, g, 10);
	}
	mpf_set (zeta, e);
}

void fp_zeta5 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
		char * g="1.036927755143369926331365486457034168057080919501912811974192677 9038035897862814845600431065571333363796203414665566090428009617 7915597084183511072180087644866286337180353598363962365128888981 3352767752398275032022436845766444665958115993917977745039244643 9196666159664016205325205021519226713512567859748692860197447984 3200672681297530919900774656558601526573730037561532683149897971 9350398378581319922884886425335104251602510849904346402941172432 7576341508162332245618649927144272264614113007580868316916497918";

		mpf_set_str (e, g, 10);
	}
	mpf_set (zeta, e);
}

void fp_zeta7 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
		char * g="1.008349277381922826839797549849796759599863560565238706417283136 5716014783173557353460969689138513239689614536514910748872867774 1984033544031579830103398456212106946358524390658335396467699756 7696691427804314333947495215378902800259045551979353108370084210 7329399046107085641235605890622599776098694754076320000481632951 2586769250630734413632555601360305007373302413187037951026624779 3954650225467042015510405582224239250510868837727077426002177100 0195455778989836046745406121952650765461161356548679150080858554";
		mpf_set_str (e, g, 10);
	}
	mpf_set (zeta, e);
}

void fp_zeta9 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
		char * g="1.002008392826082214417852769232412060485605851394888756548596615 9097850533902583989503930691271695861574086047658470602614253739 7072243015306913249876425109092948687676545396979415407826022964 1544836250668629056707364521601531424421326337598815558052591454 0848901539527747456133451028740613274660692763390016294270864220 1123162209241265753326205462293215454665179945038662778223564776 1660330281492364570399301119383985017167926002064923069795850945 8457966548540026945118759481561430375776154443343398399851419383";

		mpf_set_str (e, g, 10);
	}
	mpf_set (zeta, e);
}

/* ==================================================================== */

int main (int argc, char * argv[])
{
#ifdef FACT_TEST
	char str[4000];
	mpz_t fact;
	mpz_init (fact);

	i_factorial (fact, 5);
	mpz_get_str (str, 10, fact);
	printf ("fact = %s\n", str);
#endif

#ifdef I_BINOMIAL_TEST
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

// #define I_STIRLING_TEST
#ifdef I_STIRLING_TEST
	int n, k;
	mpz_t sitrly;
	mpz_init (sitrly);

	for (n=0; n<21; n++)
	{
		for (k=0; k<=n; k++)
		{
			i_stirling_first (sitrly, n ,k);
			mpz_get_str (str, 10, sitrly);
			printf ("sitrly (%d %d) = %s\n", n, k, str);
		}
		printf ("---\n");
	}
#endif

// #define F_BINOMIAL_TEST
#ifdef F_BINOMIAL_TEST
	int n, k;
	mpf_t bin;
	mpf_init (bin);

	for (n=1; n<7; n++)
	{
		for (k=0; k<=n; k++)
		{
			fp_binomial (bin, (double)n ,k);
			printf ("bin (%d %d) = ", n, k);
			mpf_out_str (stdout, 10, 60, bin);
			printf ("\n");
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
	/* XXX this should be stirling's approx for binomial */ 
	int bits = (int) (v + 300 + 3*nterms);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
// #define TEST_LOG
#ifdef TEST_LOG
	mpf_t lg;
	mpf_init (lg);
	mpf_set_ui (lg, 1);
	mpf_div_ui (lg, lg, 2);
	fp_log (lg, lg, prec);
	fp_prt ("log 2= ", lg);
	mpf_clear (lg);
#endif
	
// #define TEST_PI
#ifdef TEST_PI
	mpf_t pi, pis;
	mpf_init (pi);
	mpf_init (pis);
	fp_pi (pi, prec);
	fp_pi_string (pis);
	fp_prt ("pi= ", pi);
	mpf_sub (pi, pi, pis);
	fp_prt ("diff= ", pi);
	mpf_clear (pi);
	mpf_clear(pis);
#endif
	
// #define TEST_EULER
#ifdef TEST_EULER
	mpf_t gam, gams;
	mpf_init (gam);
	mpf_init (gams);
	fp_euler_mascheroni (gam, prec);
	fp_euler_mascheroni_string (gams);
	fp_prt ("gamma= ", gam);
	mpf_sub (gam, gam, gams);
	fp_prt ("diff= ", gam);
	mpf_clear (gam);
	mpf_clear(gams);
#endif
	
// #define TEST_ZETA_INT
#ifdef TEST_ZETA_INT
	mpf_t zeta, zetas;
	mpf_init (zeta);
	mpf_init (zetas);
	fp_zeta5 (zetas);
	fp_zeta (zeta, 5, prec);
	fp_prt ("zeta= ", zeta);
	mpf_sub (zeta, zeta, zetas);
	fp_prt ("diff= ", zeta);
	mpf_clear (zeta);
	mpf_clear(zetas);
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

// #define TEST_SINE
#ifdef TEST_SINE
	mpf_t ex, pi2, pi2n;
	mpf_init (ex);
	mpf_init (pi2);
	mpf_init (pi2n);
	fp_pi (pi2, prec);
	mpf_div_ui (pi2, pi2, 2);
	mpf_set (pi2n, pi2);
	
	int n;
	for (n=1; n<130; n++)
	{
		fp_sine (ex, pi2n, prec);
		if (1==n%4)
		{
			mpf_sub_ui (ex,ex,1);
		}
		if (3==n%4)
		{
			mpf_add_ui (ex,ex,1);
		}
		printf ("sin(%d pi/2)= ", n);
		fp_prt ("", ex);
		printf ("\n");

		fp_cosine (ex, pi2n, prec);
		if (0==n%4)
		{
			mpf_sub_ui (ex,ex,1);
		}
		if (2==n%4)
		{
			mpf_add_ui (ex,ex,1);
		}
		printf ("cos(%d pi/2)= ", n);
		fp_prt ("", ex);
		printf ("\n");
		printf ("\n");

		mpf_add (pi2n, pi2n, pi2);
	}

	mpf_clear (ex);
	mpf_clear(pi2);
	mpf_clear(pi2n);
#endif
	
// #define ZETA_STUFF
#ifdef ZETA_STUFF
	mpf_t zeta;
	mpf_init (zeta);
	
	printf ("           000000000011111111112222222222333333333344444444445555555555666666666677777777778\n");
	printf ("           012345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
	fp_zeta (zeta, 2, 13);
	fp_prt ("13 digs= ", zeta);
	fp_zeta (zeta, 7, 30);
	fp_prt ("30 digs= ", zeta);
	fp_zeta (zeta, 7, 40);
	fp_prt ("40 digs= ", zeta);
	fp_zeta (zeta, 7, 50);
	fp_prt ("50 digs= ", zeta);
	fp_zeta (zeta, 7, 60);
	fp_prt ("60 digs= ", zeta);
	fp_zeta (zeta, 7, 70);
	fp_prt ("70 digs= ", zeta);
	fp_zeta (zeta, 7, 80);
	fp_prt ("80 digs= ", zeta);
#endif

// #define TEST_ZETA
#ifdef TEST_ZETA
	mpf_t zeta;
	mpf_init (zeta);
	// fp_zeta_odd (zeta, 3, 180, 7, 360, 0, 60); 
	// fp_prt ("duude zeta3= ", zeta);
	int i;
	int pr = prec;
	for (i=3; i<nterms; i++ ) {
		fp_zeta (zeta, i, pr);
		printf ("char * zeta_%d_%d = \"", i, pr);
		mpf_out_str (stdout, 10, pr, zeta);
		printf ("\";\n");
		fflush (stdout);
	}
#endif

#define TEST_COMPLEX_ZETA
#ifdef TEST_COMPLEX_ZETA
	mpf_t nzeta;
	mpf_init (nzeta);
	
	cpx_t zeta, ess;
	cpx_init (&zeta);
	cpx_init (&ess);
	mpf_set_ui (ess.im, 0);

	int i;
	int pr = prec;
	for (i=3; i<nterms; i++ ) {
		fp_zeta (nzeta, i, pr);
		
		mpf_set_ui (ess.re, i);
		fp_borwein_zeta_c (zeta, ess, pr);

		printf ("zeta(%d)_%d = ", i, pr);
		cpx_prt ("", &zeta);
		printf ("  ");
		mpf_out_str (stdout, 10, pr, nzeta);
		printf ("\n");
		fflush (stdout);
	}
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

// #define TEST_BERNOULLI
#ifdef TEST_BERNOULLI
	mpq_t bern;
	mpq_init (bern);
	int n = 4;
	for (n=8; n<30; n++) 
	{
		q_bernoulli (bern, n);
		printf ("bernoulli (%d)= ", n);
		mpq_out_str (stdout, 10, bern);
		printf ("\n");
	}
#endif /* TEST_BERNOULLI */

// #define TEST_STIELTJES
#ifdef TEST_STIELTJES
	mpf_t stie;
	mpf_init (stie);
	int i;
	for (i=0; i<40; i++ ) {
		stieltjes_gamma (stie, i);
		printf ("gamma[%d] = ", i);
		mpf_out_str (stdout, 10, 60, stie);
		printf (";\n");
		fflush (stdout);
	}
#endif

// #define TEST_B_SUB_N
#ifdef TEST_B_SUB_N
	mpf_t rbs, ibs, bs;
	mpf_init (rbs);
	mpf_init (ibs);
	mpf_init (bs);
	int i;
	for (i=0; i<40; i++ ) {
		double re_s = i;
		double im_s = 0.0;
		b_sub_s_d (rbs, ibs, re_s, im_s, prec, nterms, -1.0);
		b_sub_n (bs, i, prec);
		mpf_sub (bs, bs, rbs);
		printf (" res=%g ", re_s);
		mpf_out_str (stdout, 10, 60, bs);
		printf (";\n");
		fflush (stdout);
	}
#endif

	return 0;
}

/* =============================== END OF FILE =========================== */

