
/*
 * dsubn.c
 *
 * proide analytic fit to d_n
 *
 * Linas November 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_gamma.h>
#include "harmonic.h"

void dsubn(double res, double ims)
{
	gsl_sf_result lnr, arg;
	gsl_sf_lngamma_complex_e (-res, -ims, &lnr, &arg);

	double r = exp(lnr.val);
	printf("Gamma(%g+i%g)= %15.10g * exp(i %15.10g)\n", -res, -ims, r, arg.val);

	double rez, imz, rezz, imzz;
	double h=1.0e-8;
	riemann_zeta (res+h, ims, &rez, &imz);
	riemann_zeta (res, ims, &rezz, &imzz);
	double zpre = (rez-rezz)/h;
	double zpim = (imz-imzz)/h;
	printf ("zeta-prime(%g+i%g)= %15.10g +i %15.10g\n", res, ims, zpre, zpim);

	double zpm = sqrt(zpre*zpre + zpim*zpim);
	double zarg = atan2 (zpim, zpre);

	printf ("zeta-prime(%g+i%g)= %15.10g *exp(i %15.10g)\n", res, ims, zpm, zarg);
	
	double trm = r/zpm;
	double trph = arg.val - zarg;

	printf ("mobius(%g+i%g)= %15.10g *exp(i %15.10g)\n", res, ims, 2.0*trm, trph);
	printf ("\n");
	
	/* ---------------------------------------- */
	double re_two_zeta = 1.83673535340283418798751778022;
	double im_two_zeta = -0.651197596522268672482042863805;

	double mag_two_zeta = sqrt(re_two_zeta*re_two_zeta + im_two_zeta*im_two_zeta);
	double arg_two_zeta = atan2(im_two_zeta, re_two_zeta);

	printf ("zeta(2rho) = (%g+i%g)= %15.10g *exp(i %15.10g)\n", 
	         re_two_zeta, im_two_zeta, mag_two_zeta, arg_two_zeta);
	
	double lma = trm * mag_two_zeta;
	double lph = trph + arg_two_zeta;

	printf ("liouville(%g+i%g)= %15.10g *exp(i %15.10g)\n", res, ims, 2.0*lma, lph);
	printf ("\n");
	
	/* ---------------------------------------- */
	double re_zeta_rm1 = -1.18447431294678829017538726492;
	double im_zeta_rm1 = -0.314293332466740197824844171097;

	double mag_zeta_rm1 = sqrt(re_zeta_rm1*re_zeta_rm1 + im_zeta_rm1*im_zeta_rm1);
	double arg_zeta_rm1 = atan2(im_zeta_rm1, re_zeta_rm1);

	printf ("zeta(rho-1) = (%g+i%g)= %15.10g *exp(i %15.10g)\n", 
	         re_zeta_rm1, im_zeta_rm1, mag_zeta_rm1, arg_zeta_rm1);
	
	double tma = trm * mag_zeta_rm1;
	double tph = trph + arg_zeta_rm1;

	printf ("totient(%g+i%g)= %15.10g *exp(i %15.10g)\n", res, ims, 2.0*tma, tph);
	printf ("\n");
}

void dsubn_liouville(void)
{

}
	
int
main ()
{
	double res=0.5;
	double ims = 14.13472514173469379;
	
	dsubn (res, ims);
	
	ims = 21.02203963877;

	// dsubn_mobius (res, ims);

	return 0;
}
