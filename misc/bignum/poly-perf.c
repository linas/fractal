/**
 * poly-perf.c
 *
 * Make comparative performance measurements for
 * polylog algorithms.
 *
 * January 2007
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/times.h>
#include <stdint.h>
#include <unistd.h>

#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-misc.h"
#include "mp-polylog.h"
#include "mp-trig.h"

/* ================================================================= */
/* precomputed values polylog for z=0.4+i0.3 and s=0.5+i14.134725 */
char * restr = "0.3265169662918768964336182652591014517220411368781520995604920255170776185579031804786888110719485462956930901910289111142598766376347780606708296111402670455130172866635523311896707850483146400920734480050255756189972165820840107281595908445625030532922827568067067869995676785602863335714887677984385356874048511151935924863776046693702384454246620785423336418435148737377022260531762688161413950250870186762427749762005673690824036245794922995409478494192028328328193349448816517632990656124116990467195852259279198165305951004181863536453554902228954253983520315817913319460925409498945397928685423738062742335394770507700612508344790898765965000197342009778253906646338744682729847460825372181883106496668588712169117761042712179998465843117438339904344025907258959826006642206327148231728224861631729003018397678637406511026633690031123778256597578709742358070721518970362888330766318596905323944805533563770998495828236244578494713227178453074359786653502756903512235116634630084763150174688259086367750072147932067574224267177329213120269581600804434146705616744376714607942835549398672815924093550347603548501222499845507430941793838480910629427952371668741775506199610743538136589499558091579219161412339103099370786996586834058593955365707089534350037764647711741656648131162662371300876602e0";
char * imstr = "0.1192602610093100065401128667202121398194310190397869083822554301022962127663963476629969860856346774645310837636554153820836650352254973136373745479380716188578359957969156436051959717578495535326786996267206207968058305242211626494468971345045328425252794970842575055127083991554693827561407342426007184204626603453209593878448357514248832368613553809324048616983515096401107070047816712040680356822916688838162911562419765681249059889947272646613432129072090060594805553609533995671865746923135430496474419310866698823653980107352715160049770909462201812117919355859225869123736848825457105175520071523186063051675724019495962679089580812155008410924222866654109902444393435775620497178304606368060071069089910041060084225992496194998419582883110686436188289148226627434439403557696516251684172609878478372069415284200446520758194942843877901562776022247543339092760971238721806794328289680986892351895486016174475181579059906946095496154861332507882761692487809393309760470128252160424420299372082057390365644069429920312091399063297952529850607665328587938511544943896083706291775141843254472323107451978563165941640384895428293821158337192293459314822720260766163151776916664823523917657406939010521473306234033585759522851577721304932894524474882704589419966446750663524251965193748052550969084e0";
/* ================================================================= */

void do_perf(void)
{
	int i;
	int prec = 40;
	int rex;

	printf ("#\n# graph of periodic zeta as function of precision \n#\n");
	int hz = sysconf (_SC_CLK_TCK);
	printf ("# clock ticks=%d\n#\n", hz);
	fflush (stdout);

	for (prec=10; prec <123123; prec *= 1.41)
	{
		/* Set the precision (number of binary bits) */
		int nbits = 3.32*prec+prec+100;
		mpf_set_default_prec (nbits);

		printf ("%d\t%d\t", prec, nbits);

		cpx_t ess, zeta, zee, plog, expected;
		cpx_init (ess);
		cpx_init (zeta);
		cpx_init (zee);
		cpx_init (plog);
		cpx_init (expected);
	
		mpf_t que;
		mpf_init (que);

		mpf_set_str (expected[0].re, restr, 10);
		mpf_set_str (expected[0].im, imstr, 10);
	
		struct tms start, end;

#define MEASURE_POLYLOG_PERFORMANCE 1
#ifdef MEASURE_POLYLOG_PERFORMANCE
		cpx_set_d (ess, 0.5, 14.134725);
		cpx_set_d (zee, 0.4, 0.3);

#if 0
		/* First we warm the cache */
		times (&start);
		cpx_polylog (zeta, ess, zee, prec);
		times (&end);

		cpx_sub (zeta, zeta, expected);
		rex = get_prec (zeta, prec);
		printf ("%d\t", rex);

		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
		
		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_polylog (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		cpx_sub (zeta, zeta, expected);
		rex = get_prec (zeta, prec);
		printf ("%d\t", rex);
#endif

#if 0
		/* First we warm the cache */
		times (&start);
		cpx_polylog_euler (zeta, ess, zee, prec);
		times (&end);

		cpx_sub (zeta, zeta, expected);
		rex = get_prec (zeta, prec);
		printf ("%d\t", rex);

		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
		
		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_polylog_euler (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		cpx_sub (zeta, zeta, expected);
		rex = get_prec (zeta, prec);
		printf ("%d\t", rex);
#endif

#if 1
		/* First we warm the cache */
		times (&start);
		cpx_polylog_sum (zeta, ess, zee, prec);
		times (&end);

		cpx_sub (zeta, zeta, expected);
		rex = get_prec (zeta, prec);
		printf ("%d\t", rex);

		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		/* Then with a hot cache */
		times (&start);
		for (i=0; i<100; i++)
			cpx_polylog_sum (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		cpx_sub (zeta, zeta, expected);
		rex = get_prec (zeta, prec);
		printf ("%d\t", rex);
#endif

		printf ("\n");
		fflush (stdout);
#endif

// #define MEASURE_HURWITZ_PERFORMANCE
#ifdef MEASURE_HURWITZ_PERFORMANCE
		cpx_set_d (ess, 0.5, 14.134725);
		cpx_set_d (zee, 0.2, 0.0);
		mpf_set_d (que, 0.2);

#if 0
		/* First we warm the cache */
		times (&start);
		cpx_hurwitz_zeta (zeta, ess, que, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_hurwitz_zeta (zeta, ess, que, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
#endif

#if 1
		/* First we warm the cache */
		times (&start);
		cpx_hurwitz_taylor (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));

		/* Then with a hot cache */
		times (&start);
		for (i=0; i<1000; i++)
			cpx_hurwitz_taylor (zeta, ess, zee, prec);
		times (&end);
		printf ("%jd\t", (intmax_t) (end.tms_utime - start.tms_utime));
#endif
		printf ("\n");
#endif

		cpx_clear (ess);
		cpx_clear (zeta);
		cpx_clear (zee);
		cpx_clear (plog);
		cpx_clear (expected);
		mpf_clear (que);

		fflush (stdout);
	}
}

/* ================================================================= */


/* Compute to high precision, monitor the precision */
void hiprec(cpx_t zeta, int prec)
{
	/* Set the precision (number of binary bits) */
	int nbits = 3.322*prec+2*3.14*3.32*prec+140;
	nbits = 3.322*prec+ 7*prec + 140;
	// nbits = 3.322*prec+ 140;
	nbits = 3.322*prec+ prec+ 140;
	mpf_set_default_prec (nbits);
	cpx_set_prec (zeta, nbits);

	cpx_t cq, ess, prevzeta;
	cpx_init (ess);
	cpx_init (cq);
	cpx_init (prevzeta);
	cpx_set (prevzeta, zeta);

	cpx_set_d (ess, 0.5, 14.134725);
	cpx_set_d (cq, 0.4, 0.3);
	
	// cpx_polylog_euler (zeta, ess, cq, prec);
	cpx_polylog (zeta, ess, cq, prec);
	cpx_sub (ess, zeta, prevzeta);

	printf ("prec=%d ", prec);

	long rex = get_prec (ess, prec);
	printf ("prev=%ld ", rex);

#if 0
	double vre = mpf_get_d (ess[0].re);
	printf ("vre = %g ", vre);

	double vim = mpf_get_d (ess[0].im);
	printf ("vim = %g ", vim);
#endif

#if 1
	// gmp_printf ("re= %.220Ff ", zeta[0].re);
	printf ("re= ");
	mpf_out_str (stdout, 10, prec, zeta[0].re);
	
	printf (" im= ");
	mpf_out_str (stdout, 10, prec, zeta[0].im);
#endif
	printf ("\n");
	fflush (stdout);

	cpx_clear (ess);
	cpx_clear (cq);
	cpx_clear (prevzeta);
}

/* ================================================================= */

int
main (int argc, char * argv[])
{
#if 0
	if (argc != 2)
	{
		fprintf (stderr, "Usage: %s <sim>\n", argv[0]);
		exit (1);
	}
	double sim = atof (argv[1]);
#endif

#if 0
	cpx_t zeta;
	cpx_init (zeta);
	int prec;
	for (prec=10; prec<15001; prec *= 1.4)
	{
		hiprec(zeta, prec);
	}
#endif
	
	do_perf();
	
	return 0;
}
