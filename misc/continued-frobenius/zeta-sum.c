
/* zeta-sum.c
 *
 * Goal: numerically validate a sum over zeta's
 * */

#include <math.h>
#include <stdio.h>
#include "zetafn.h"

main ()
{
	int i;
	long double sign = 1.0L;
	long double acc = 0.0;
	for (i=0; i<70; i++)
	{
		long double val;
		val = zetam1(i+2);
		val *= (val - 2.0L);
		long double eye = (long double ) (i+1);
		val = sign * eye *val;
		acc += val;
		sign = - sign;
		printf ("%d %30.25Lg\n", i, acc);
	}
	
	long double t;
	for (t=0.25; t>1.0e-5; t *= 0.5)
	{
	acc = 0.0;
	sign = 1.0L;
	for (i=0; i<136470; i++)
	{
		long double eye = (long double ) i;
		long double reg = - t*t*eye*eye;
		reg = expl (reg);
		long double val;
		val = sign * eye;
		acc += val *reg;
		sign = - sign;
	}
	printf ("%Lg %30.25Lg\n", t, acc);
	}

	int m = 0;
	int k;



#ifdef WHATEVER
	for (m=0; m<20; m++)
	{
	double t = 1.0;
	for (t=0.25; t>1.0e-5; t *= 0.5)
	{
		double fc_m = 1.0;         // factorial m
		for (k=2; k<m+1; k++)
		{
			fc_m *= (double) k;
		}
		double fc_km = fc_m;  // factorial k+m+1

		double fc_kp1 = 1.0;  // factorial k+1
		double s = 1.0;
		
		double fcr = fc_m;   // (k+m+1)! / (k+1)!
		double acc = 0.0;
		for (k=0; k<10460; k++)
		{
			double reg = exp (-t*t * k*k);
			fc_kp1 *= ((double) k+1);
			fc_km *= ((double) k+m+1);
			fcr /= ((double) k+1);
			fcr *= ((double) k+m+1);
			
			// printf ("k=%d fc_k=%f fc_km=%f\n", k, fc_k, fc_km);
			// acc += reg * s * fc_km * (1.0 + zetam1(k+m+2) ) / fc_kp1;
			if (70 > k) {
				acc += reg * s * fcr * (1.0 + zetam1(k+m+2));
			} else {
				acc += reg * s * fcr;
			}
			s = -s;
		}
		acc /= fc_m;
		printf ("m=%d tee=%f acc=%f\n", m, t, acc);
		fflush (stdout);
	}
	printf ("====================================\n");
}
#endif

#ifdef SOMETHING
	double t = 1.0;
	for (t=0.25; t>1.0e-9; t *= 0.5)
	{
		int k;
		double ss = 1.0;
		double lam = 0.0;
		double fk = 1.0;  // factorial k
		for (k=0; k<40; k++)
		{
			double reg = exp (-t*t*k*k);
			double term = (zetam1(k+2)+1.0);
			term *= term;
			term -= 1.0;
			term *= ss * ((double) (k+1));
			lam += reg * term;
			ss = -ss;
			// printf ("k=%d term=%f sum=%f\n", k, term, lam); 
		}
		lam += 0.25;
		lam /= zetam1(2) + 1.0;
		printf ("t=%f lambda = %18.12g\n", t, lam);
		// printf ("========================\n");
	}
#endif
	
#if 0
	int m;
	double fm = 1.0;  // factorial m
	for (m=0; m<4; m++) 
	{
		if (m!=0) fm *= m;
		for (t=1.0; t>1.0e-6; t *= 0.5)
		{
			int k;
			double sign = 1.0;
			
			double acc = 0.0;
			double fk = 1.0;
			double fkm = fm;
			for (k=0; k<40; k++)
			{
				if (k!=0) fk *= k;
				if (k!=0) fkm *= k+m;

				double reg = exp (-t*t*k*k);
				double par = (1.0+zetam1(k+m+2))*(1.0+zetam1(k));
				par *= sign / fk;
				par *= fkm * (k+m+1);
				par /= (m+1)*fm;
				par /= zeta[m+2];
				
				acc += reg * par;
				sign = -sign;
			}
			printf ("m=%d t=%f orth = %f\n", m, t, acc);
		}
	}
#endif
	

}
