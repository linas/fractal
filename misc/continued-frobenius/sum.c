
/*
 * sum.c
 *
 * validate some identities at the y=1-x expansion
 */
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_sf_zeta.h>

#include "zetafn.h"
#include "ache.h"


int
main (int argc, char *argv[])
{
	int m,n,p;

	n=0;
	if (argc==2) 
	{
		n = atoi(argv[1]);
	}

#if PROVEN_SUM
	for (m=8; m<20; m++)
	{
	for (k=8; k<20; k++)
	{
		int p;
		long double acc = 0.0L;
		acc = factorial (k+m+1) * (zetam1(k+m+2)+1.0L);
		acc = -acc / factorial (m);
		for (p=m; p<160; p++)
		{
			long double val = binomial (p,m);
			val *= factorial (p+k+1) / factorial (p);
			val *= zetam1 (p+k+2);
			acc += val;
			// printf ("m=%d k=%d p=%d acc=%Lg\n", m,k,p,acc);
		}
		printf ("m=%d k=%d acc=%Lg\n", m,k,acc);
	}}
#endif 
	
#ifdef MATRIX_ELT
	double trix[15][15];
	for (m=0; m<15; m++)
	{
	for (k=0; k<15; k++)
	{
		int p;
		long double acc = 0.0L;
		long double sign = 1.0L;
		for (p=0; p<=k; p++)
		{
			long double val = binomial (k,p);
			val *= binomial (m+k+1, m);
			val *= zetam1 (p+m+2);
			acc += sign * val;
			sign = -sign;
			// printf ("m=%d k=%d p=%d acc=%Lg\n", m,k,p,acc);
		}
		printf ("m=%d k=%d acc=%Lg\n", m,k,acc);
		trix[m][k] = acc;
	}
	printf ("---\n");
	}
	printf ("++++++++++++++++++++++++++++\n");

	for (m=0; m<15; m++)
	{
	for (k=m; k<15; k++)
	{
		printf ("m=%d k=%d rat=%19.17g\n", m,k, trix[k][m]/trix[m][k]);
	}
	printf ("---\n");
	}
#endif
	
#ifdef MORE_MATRIX_ELT
	k=n;
	long double len = 0.0;
	for (m=k; m<150; m++)
	{
		int p;
		long double acc = 0.0L;
		long double sign = 1.0L;
		for (p=0; p<=k; p++)
		{
			long double val = binomial (k,p);
			val *= binomial (m+k+1, m);
			val *= zetam1 (p+m+2);
			acc += sign * val;
			sign = -sign;
			// printf ("m=%d k=%d p=%d acc=%Lg\n", m,k,p,acc);
		}
		len += acc*acc;
		printf ("m=%d k=%d acc=%Lg len=%Lg\n", m,k,acc, len);
	}
	printf ("---\n");
#endif
	
#if ILL_CONDITIONED_BUT_WORKS
	int n;
	for (n=0; n<20; n++)
	{
	long double t;
	for (t=0.8; t>1e-5; t *= 0.8)
	{
		int p;
		long double acc = 0.0L;
		// acc = -(zetam1(n+2)+1.0L);
		acc = -1.0L;
		acc = 0.0L;
		long double sign = 1.0L;
		for (p=0; p<160; p++)
		{
			long double val = binomial (p+n+1,p);
			val *= zetam1 (p+n+2) +1.0L;
			long double reg = (long double) p;
			reg = expl (-t*t*reg*reg);
			acc += sign*val*reg;
			sign = -sign;
		}
		printf ("-- n=%d t=%Lg acc=%Lg\n", n,t, acc);
	}
	}
#endif

#if BINOMIAL_IDENT
	int n,p;
	for (m=0; m<15; m++)
	{
		for (p=0; p<15; p++)
		{
			double acc = 0.0;
			double sign = 1.0;
			for (n=0; n<=p; n++)
			{
				acc += binomial (p,n) * binomial (m+n+1,n+1) * sign;
				sign = -sign;
			}
			printf ("m=%d p=%d sum=%f gu=%f\n", m,p,acc, (double)binomial (m+1,p+1));
		}
		printf ("---\n");
	}
#endif
	
#ifdef GENERALIZED_BINMOIAL_THEORM
	long double t = 0.0;
	long double pw = 0.97314;
	for (t=0.0; t<1.0; t+=0.06L)
	{
		long double acc = 0.0L;
		long double tn = 1.0;
		for (n=0; n<130; n++)
		{
			acc += tn * cbinomial (pw, n);
			tn *= t;
		}
		acc -= powl (1.0+t, pw);
		
		printf ("its binmoial t=%Lg acc=%Lg\n",t, acc); 
	}
	
	for (t=0.1; t<21.0; t+=0.1L)
	{
		long double acc = 0.0L;
		for (n=0; n<230; n++)
		{
			acc += cbinomial (t, n);
		}
		acc -= powl (2.0L, t);
		
		printf ("its pow2 t=%Lg acc=%Lg\n",t, acc); 
	}

	for (k=0; k<15; k++)
	{
		for (m=0; m<15; m++)
		{
			long double complex ess = m+0.5 +3.14I;
			ess = m;
			ess += 4.5+13.14LI;
			long double complex acc = 0.0L;
			for (n=k; n<260; n++)
			{
				long double complex term = cbinomial (ess, n);
				term *= binomial (n,k);
				acc += term;
				printf ("duude %d term=%Lg\t+i%Lg\t  acc=%Lg\t+i%Lg\n", n, creall(term), cimagl(term), creall(acc), cimagl(acc));
			}
			// acc -= binomial (m,k) * pow (2.0, m-k);
			long double complex sk = ess - (long double) k;
			acc -= cbinomial (ess,k) *cpowl (2.0L, sk);
			printf ("k=%d s=%Lg  acc=%Lg \n", k,creall(ess),creall(acc));
		}
		printf ("------\n");
	}
#endif

#if CHECK_MATH
	for (k=1; k<15; k++)
	{
		long double acc = 0.0L;
		for (m=0; m<130; m++)
		{
			long double val = binomial (m+k+1, m);
			val *= zetam1 (m+k+2);
			val /= (long double) ((m+1)*(m+2));
			acc += val;
		}
		long double g = -zetam1(k+1);
		g += 1.0L/ ((long double) k);
		g /= (long double) (k+1);
		acc -= g;
		printf ("k=%d acc=%Lg\n", k,acc);
	}
#endif

#ifdef EULERS_GAMMA
	long double acc = 0.0L;
	for (m=0; m<170; m++)
	{
		acc += zetam1(m+2) / ((long double) (m+2));
	}
	acc = 1.0L-acc;
	printf ("acc=%30.25Lg\n", acc);
#endif

#ifdef A_SUB_N_TERMS

	long double acc_prev = 1.0L;
	for (n=0; n<50; n++)
	{
		long double acc = 0.0L;
		for (m=0; m<190; m++)
		{
			long double term;
			term = ache_mp (m,n);
			term /= (long double) (m+1)*(m+2);
			acc += term;
		}
		// acc *= 2.0L* (long double) (n+1);
		acc -= 0.5L/ (long double) (n+1);
		// acc = a_sub_n_poor_convergence (n);
		// acc /= expl (-4.0L*sqrtl(n+1)); 
		acc *= exp (((double)(n+1))/log((double)(n+1)));
		// printf ("n=%d term=%20.15Lg\n", n, acc);
		printf ("%d\t%20.15Lg\n", n, acc);
		acc_prev = acc;
	}
#endif

#ifdef ZETA_IN_TERMS_OF_H
	long double ess=3.0L;

	for (ess=-1.05L; ess<15.0; ess+=0.2L)
	{
		long double complex acc = 0.0L;
		for (m=0; m<90; m++)
		{
			long double sn = 1.0L;
			for (n=0; n<90; n++)
			{
				long double complex term;
				term = ache_mp (m,n);
				term /= (long double) (m+1)*(m+2);
				term *= cbinomial (ess-1.0L, n);
				acc += sn*term;
				//long double complex z = -acc;
				//z += 1.0L / (ess-1.0L);
				// z *= ess;
				// printf ("m=%d n=%d term= %Lg  z=%Lg\n",
				//		m,n, creall(term), creall(z));
				sn = -sn;
			}
			// printf ("--- \n");
		}
		acc = -acc;
		acc += 1.0L / (ess-1.0L);
		acc *= ess;
		double calib_zeta = gsl_sf_zeta(ess);
		double delta = acc - calib_zeta;
		printf ("ess=%Lg result=%Lg calib=%g delta = %g\n", 
			 ess, creall(acc), calib_zeta, delta);	
	}
#endif

#ifdef ZETA_IN_TERMS_OF_H_HALF_SOLVED
	long double ess=3.0L;

	for (ess=-1.05L; ess<15.0; ess+=0.2L)
	{
		long double complex acc = 0.0L;
		long double sn = 1.0L;
		for (n=0; n<130; n++)
		{
			long double complex term;
			term = a_sub_n (n);
			term *= cbinomial (ess-1.0L, n);
			acc += sn*term;

#if 0
			long double complex z = -acc;
			z += 1.0L / (ess-1.0L);
			z *= ess;
			z -= 0.5;
			printf ("n=%d term= %Lg  z=%Lg\n", n, creall(term), creall(z));
#endif
			sn = -sn;
		}
		acc = -acc;
		acc += 1.0L / (ess-1.0L);
		acc *= ess;
		acc -= 0.5;
		double calib_zeta = gsl_sf_zeta(ess);
		double delta = acc - calib_zeta;
		printf ("ess=%Lg result=%Lg calib=%g delta = %g\n", 
			 ess, creall(acc), calib_zeta, delta);	
	}
#endif

#ifdef ALTERNATING_BINOMIAL
	long double ks = 1.0L;
	for (k=5; k<10; k++)
	{
		long double ess = 0.0L;
		for (ess = 0.5; ess<5.0; ess+=0.2L)
		{
			long double sign = ks;
			long double acc = 0.0L;
			for (n=k; n<230; n++)
			{
				long double term;
				term = binomial (n,k);
				term *= cbinomial (ess, n);
				term *= sign;
				acc += term;
				// printf ("duude k=%d s=%Lg n=%d term=%Lg acc=%Lg\n", k,ess,n,term,acc);
				sign = -sign;
			} 
		printf ("duude k=%d s=%Lg acc=%Lg\n", k,ess,acc);
		}
		ks = -ks;
	}
#endif
	
#ifdef EXPLOSIVE
	long double ess = 0.0L;
	int nterms = n;
	for (ess = -2.0; ess<=0.0; ess+=0.1L)
	{
		long double sign = 1.0;
		long double acc = 0.0L;
		for (n=0; n<nterms; n++)
		{
			long double term;
			term = cbinomial (ess, n);
			term *= sign;
			acc += term;
			// printf ("duude s=%Lg n=%d term=%Lg acc=%Lg\n", ess,n,term,acc);
			sign = -sign;
		} 
		printf ("duude s=%Lg acc=%Lg \n", ess,acc);
	}
#endif

#ifdef PRINT_OUT_TEE_SUB_N
	for (n=0; n<40; n++)
	{
		long double sk = -1.0L;
		long double acc = 1.0L - M_GAMMA;
		for (k=1; k<=n; k++)
		{
			long double term = 1.0L/ (long double) k;
			term -= zetam1(k+1);
			term /= (long double) (k+1);
			term *= binomial (n,k);
			term *= sk;
			acc += term;
			// printf ("duuude k=%d term=%Lg, acc=%Lg\n", k, term, acc);
			sk = -sk;
		}
		// printf ("for n=%d term=%12.10Lg\n",n, ((acc*((long double)n)-0.5L)*n+0.5L)*n);
		long double en = 1.0L/((long double) n);
		long double sympa = 0.0L;
		// sympa = 1.0L -en *(1.0L-en*(1.0L-en*(1.0L-en*(1.0L-en))));
		// sympa = 1.0L -en *(1.0L-en*(1.0L-en*(1.0L-en*(1.0L-en*sympa))));
		// sympa *= 0.5*en;
		sympa = 0.5L/((long double) (n+1));
		
#ifdef COFU
		// the continued fraction approx sucks
		long double asymp = 0.0L;
		asymp = 1.0L +en;
		for (k=0;k<20; k++)
		{
			asymp = 1.0L + en/asymp;
		}
		asymp = 0.5*en / asymp;
#endif
		
		long double guess = 1.0L/factorial (n+1);
		printf ("for n=%d term=%12.10Lg sa=%Lg\n",n, acc-sympa, guess );
	}
#endif
	
#ifdef ZETA_IN_TERMS_OF_H_IMAG_AXIS
	long double complex ess=0.5L;

	for (p=0; p<40; p++, ess+=1.0LI)
	{
		long double complex acc = 0.0L;
		long double sn = 1.0L;
		for (n=0; n<122130; n++)
		{
			long double complex term;
			term = a_sub_n (n);
			term *= cbinomial (ess-1.0L, n);
			acc += sn*term;

#if 0
			long double complex z = -acc;
			z += 1.0L / (ess-1.0L);
			z *= ess;
			z -= 0.5;
			printf ("zeta term n=%d term= %Lg  z=%Lg\n", n, cabsl(term), cabsl(z));
#endif

			double worst = cabsl(term);
			if (1.0e-6 > worst) break;
			if (10.0 < worst) break;
			sn = -sn;
		}
		acc = -acc;
		acc += 1.0L / (ess-1.0L);
		acc *= ess;
		acc -= 0.5;
		printf ("%Lg\t%Lg\t%Lg\t%Lg\t%d\n", 
			cimagl(ess), creall(acc), cimagl(acc), cabsl(acc), n);
		fflush (stdout);
	}
#endif

#ifdef CHECK_BINOMIAL_SUM_FOR_1_OVER_S
	long double complex ess=6.0;

	ess = 0.5L;
	for (p=0; p<40; p++, ess+= 0.1LI)
	{
		long double complex acc = 0.0L;
		long double sign = 1.0L;
		for (n=0; n<31130; n++)
		{
			acc += sign * cbinomial (ess-1.0L,n) / (long double) (n+1);
			sign = -sign;
		}
		acc -= 1.0/ess;
		printf ("ess=%Lg+i%Lg acc=%Lg\n", creall(ess), cimagl(ess), acc);
	}
#endif

#ifdef VERIFY_MORE_BINOMIAL_SUMS
	for (n=1; n<100; n++)
	{
		double acc = 0.0;
		double sign = -1.0;
		for (k=1; k<=n; k++)
		{
			acc += sign * binomial (n,k) / (double) (k*(k+1));
			sign = -sign;
		}
		double yack = 0.0;
		for (k=1; k<=n; k++) 
		{
			yack -= 1.0/(double) (k+1);
		}
		printf ("its %d, acc=%g, yack=%g\n", n, acc, acc-yack);
	}
#endif

#ifdef GUESS_EIGENVECTOR
	int nterms;
	long double vlen = 0.0;
	long double vec[1350];
	double sign = 1.0;
	for (p=0; p<150; p++)
	{
		// vec[p] = 1.0L /(long double) ((p+1)*(p+1));
		// vec[p] = zetam1(p+2);
		// vec[p] = 1.0L/pow(2, p);
		// vec[p] = (p+1)/pow(3.0, p+1);
		// vec[p] = (p+1)*zetam1(p+2);
		// vec[p] = (p)*zetam1(p+2);
		// vec[p] = sqrt(p+1)*zetam1(p+2);
		// vec[p] = (p+1)*zetam1(p+2);
		// vec[p] = (p+1)*((1.0L/pow(2.0, p+2)) + (1.0L/pow(3.0, p+2)) + (1.0L/pow(4.0, p+2)) );
		// vec[p] = sqrt(2*p+3)*zetam1(p+2);
		// vec[p] = zetam1(p+2)*pow(p+3, 0.66);
		// vec[p] = zetam1(p+3)*pow(p+2, 0.5);
		// printf ("duude %d %Lg\n", p, vec[p]);
		// vec[p] = zetam1(p+2) +zetam1(p+3);
		// vec[p] = creall(cbinomial ((long double complex) (1.5+(double)p), p));
		// vec[p] /= pow (2.0, 2.5+p);
		// vec[p] *= 1.75;
		// long double complex ess;
		// ess = 0.5L;
		// ess += I * 14.1347L;
		// vec[p] = sign * creall (cbinomial(ess-1.0L, p));
		// vec[p] *= 1.0/(p*p*p*p);
		// printf ("v[%d] = %15.10Lg;\n", p, vec[p]);
		// printf ("%d\t%15.10Lg;\n", p, vec[p]);
		// vec[p] = 1.0L /powl(2, p+1);
		vec[p] = 3.0L /powl(2, p+1);
		if (0 == p) vec[0] -= 2.0;
		sign = -sign;

		vlen += vec[p] * vec[p];

		if (1.0e-16 > fabsl(vec[p]))
		{
			nterms = p+1;
			printf ("evaluate guess to %d terms\n", p+1);
			printf ("vec length squared = %Lg\n", vlen);
			break;
		}
	}
	// vec[0] = -0.75 + 1.75/ pow (2.0, 2.5);
	// printf ("v[0] = %15.10Lg;\n", vec[0]);
	
	for (m=0; m<20; m++)
	{
		long double acc = 0.0;
		long double term;
		for (p=0; p<nterms; p++)
		{
			term = ache_mp (m,p) * vec[p];
			acc += term;
		}
		acc /= vec[m];
		printf ("%d lam=%Lg last term=%Lg\n", m, acc, term);
	}
#endif

#ifdef ITERATE_ACHE_ON_CANDIDATE_VECTOR
	int i, nterms=0;
	long double vlen = 0.0;
	long double vec[1350];
for (n=0; n<20; n++) {
	long double vv = 1.0L/sqrtl( powl (4.0, n+1) -3.0);
	
	for (p=0; p<150; p++)
	{
		vec[p] = -3.0 * powl(2.0, -(p+1));
		if (n == p) vec[n] += pow(2.0,n+1);
		vec[p] *= vv;

		vlen += vec[p] * vec[p];

		if (1.0e-16 > fabsl(vec[p]))
		{
			nterms = p+1;
			// printf ("n=%d evaluate guess to %d terms\n", n, p+1);
			// printf ("vec length squared = %Lg\n", vlen);
			break;
		}
	}
	
	long double ratio, dot, cc;
	for (i=0; i<19; i++)
	{
		long double hold[1350];
		for (m=0; m<nterms; m++)
		{
			long double acc = 0.0;
			long double term;
			for (p=0; p<nterms; p++)
			{
				term = ache_mp (m,p) * vec[p];
				acc += term;
			}
			hold[m] = acc;
		}
		ratio = hold[0] / vec[0];
		vlen = 0.0;
		dot = 0.0;
		for (m=0; m<nterms; m++)
		{
			vec[m] = hold[m];
			vlen += vec[m]*vec[m];
			dot += vec[m]*powl (2.0, -(m+1));
		}
		dot *= sqrt(3.0);
		cc = vec[0];
		// printf ("const=%Lg 1/const=%Lg    ratio-1=%Lg   1/len=%Lg   dot=%Lg \n", cc, 1.0L/cc, ratio-1.0L, 1.0L/sqrtl(vlen), dot);
	}
	printf ("%d const=%Lg 1/const=%Lg    ratio-1=%Lg   1/len=%Lg   dot=%Lg \n", n, cc, 1.0L/cc, ratio-1.0L, 1.0L/sqrtl(vlen), dot);
}
	for (m=0; m<20; m++)
	{
		// printf ("Post iter %d = %Lg\n", m, vec[m]/ vec[0]);
	}
#endif

#if MORE_BINOMIAL_SUMS
	for (k=0; k<30; k++)
	{
		double acc = 0.0;
		for(p=0; p<130; p++)
		{
			double term = binomial (p+k,k) * zetam1 (p+k+2);
			acc += term;
		}
		// double guess = zetam1 (2); // for k=1
		// double guess = 1.0 -zetam1(2) + zetam1(3);  // for k=2
		// double guess = 1.0 -(zetam1(2)+1.0) + (zetam1(3)+1.0) - (zetam1(4)+1.0);
		// guess = -guess;
		double guess = 1.0;
		double sign = -1.0;
		for (p=1; p<=k; p++)
		{
			guess += sign * (zetam1(p+1)+1.0);
			sign = -sign;
		}
		if (k%2 == 1) guess = -guess;
		acc -= guess;
		printf ("its %d %g %g\n", k, acc, guess);
	}
#endif

#if 0
	k = 2;
	double acc = 0.0;
	double sign = 1.0;
	for (n=0; n<60; n+=2)
	{
		double term = binomial (n+k, k) * zetam1 (n+k+2);
		// term *= sign;
		acc += term;
		sign = -sign;
	}
	printf ("k=%d acc=%g\n", k, acc);
#endif
	
#ifdef SHOW_BINOMIALS
	long double sign = 1.0L;
	for (p=0; p<50; p++)
	{
		int i;
		printf ("%d", p);
		for (i=0; i<10; i++)
		{
			long double val;
			long double complex ess;
			ess = 0.5L;
			ess += I * (long double) (i+1);
			val = sign * creall (cbinomial(ess-1.0L, p));
			val /= factorial (i+1);
			printf ("\t%Lg", val);
		}
		printf ("\n");
		sign = -sign;
	}
#endif

#if CHECK_LEFT_VECTOR
	complex long double ess = 10.5+14.0*I;
	long double sign = 1.0;
	complex long double acc = 0.0;
	complex long double term;
	for (m=0; m<450; m++)
	{
		term = cbinomial (ess-1.0,m);
		term /= (m+1)*(m+2);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	acc -= 1.0L / (ess+1.0L);
	printf ("its %Lg +I %Lg last term=%Lg\n", creall(acc), cimagl(acc), term);
#endif

#if SIMPLE_BIN
	long double acc = 0.0;
	k = 6;
	for (n=k; n<130; n++)
	{
		acc += binomial (n,k) * powl (2.0L, -(n+1));
	}
	// acc *= powl (2.0L, (k+1));
	printf ("its k=%d acc=%Lg\n", k, acc);
#endif
	
#ifdef ALTERNATING_SIGN_ZETA_SUM
	for (m=0; m<90; m++)
	{
	long double acc = 0.0;
	long double sign = 1.0;
	for (k=0; k<460; k++)
	{
		acc += sign * binomial (m+k+1,k+1) * zetam1(k+m+2);
		sign = -sign;
	}
	acc *= powl (2.0L, (m+1));
	printf ("its m=%d acc=%Lg\n", m, acc);
	}
#endif

#ifdef FRACTIONAL_EMM
	long double emm;
	for (emm=0.0; emm<20.0; emm +=0.12345){
	long double acc = 0.0;
	long double sign = 1.0;
	for (k=0; k<160; k++)
	{
		long double zz = gsl_sf_zeta (emm+k+2) -1.0L;
		acc += sign * fbinomial (emm+k+1,k+1) * zz;
		// sign = -sign;
	}
	// acc *= powl (2.0L, (emm+1));
	printf ("its m=%Lg acc=%Lg\n", emm, acc);
	}
#endif

#if junk
	long double acc = 0.0;
	long double sign = 1.0;
	for (m=0; m<60; m++)
	{
		long double term = powl (0.5L, m+1);
		term /= (m+1)*(m+2);
		acc += sign * term;
		// sign = -sign;
	}
	printf ("its m=%d acc=%Lg\n", m, acc);
#endif

#if 0
	long double acc = 0.0;
	long double term;
	for (m=0; m<61; m++)
	{
		// long double term = powl (0.5L, m);
		term = 1.0L;
		term *= a_sub_n(m);
		acc += term;
	}
	// acc -= 2.0L + 3.0L *logl (2.0L);
	acc -= logl (sqrtl(2.0*M_PI)) - 1.0;
	printf ("its m=%d acc=%Lg last term=%Lg\n", m, acc, term);
#endif
	
#ifdef SHOW_FROBENIUS_MEASURE
	long double acc = 0.0;
	for (m=0; m<200; m++)
	{
		for (k=0; k<=m; k++)
		{
			long double term = ache_mp (k, m-k);
			acc += term*term;
		}
		printf ("%d meas = %Lg\n", m, acc);
	}
#endif

#ifdef H_H_TRANSPOSE
	for (m=0; m<10; m++)
	{
		for (n=0; n<10; n++)
		{
			long double acc = 0.0;
			long double term;
			for (k=0; k<70; k++)
			{
				term = ache_mp (m, k);
				term *= ache_mp (k,n);
				acc += term;
			}
			// acc = ache_mp(m,n);
			long double app = acc - pow (0.5, ((double)m+n));
			printf ("m[%d][%d] = %Lg   nus=%Lg  last term = %Lg\n", m, n, acc, app, term);
		}
		printf ("\n");
	}
#endif

#if 0
	// for (n=0; n<10; n++)
	{
		long double acc = 0.0L;
		long double sign = 1.0L;
		for (k=0; k<60; k++)
		{
			// acc += sign * zetam1(k+2) * (k+1)* (zetam1(k+2)+1.0);
			acc += sign * zetam1(k+2) *  pow (0.5, k);
			sign = -sign;
		}
		printf ("%Lg   \n",  acc);
	}
#endif

#if GUESS_EIGENVALUE
	long double acc = 0.0;
	long double term = 1.0L/3.0L;
	long double sign = 1.0L;
#if 0
	int i;
	for (i=0; i<10; i++)
	{
		acc += sign*term;
		printf ("its %d %Lg    %Lg\n", i, term, acc);
		term *= 1.0/9.0;
		sign = -sign;
	}
#endif

	acc = 0.0;
	term = 1.0/3.0;
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
	
	term = -1.0/27.0;
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
	
	term = 2.0/(9.0*27.0);
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
	
	term = -4.0/(81.0*27.0);
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
	
	term = 8.0/(9.0*81.0*27.0);
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
	
	term = -16.0/(81.0*81.0*27.0);
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
	
	term = 32.0/(9.0*81.0*81.0*27.0);
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
#endif // GUESS_EIGENVALUE
	
#if 0
	acc = 0.0;
	term = 1.0/3.0;
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
	
	term = -1.0/25.0;
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
	
	term = 2.0/(9.0*25.0);
	acc += term;
	printf ("its %Lg    %Lg\n", term, acc);
#endif

#ifdef DOUBLE_ZETA
	for (m=0; m<30; m++)
	{
		long double acc = 0.0L;
		for (k=0; k<50; k++)
		{
			acc += binomial (m+k+1, m)*zetam1(m+k+2)*(k+1)*zetam1(k+2);
		}
		printf ("its m=%d %20.15Lg   %20.15Lg  %20.15Lg \n", m, acc, acc*pow(2.0,(m+1)/2.0), ache_mp(m,0));
	}
#endif

#ifdef ANOTHER_BINOMIAL
	for (n=0; n<10; n++)
	{
		for (k=0; k<10; k++)
		{
			long double sign = 1.0L;
			long double acc = 0.0L;
			int j;
			for (j=0; j<=n; j++)
			{
				acc += sign * binomial (n,j)* binomial (j+k+1, k);
				sign = -sign;
			}
			sign = 1.0;
			if (n%2==1) sign = -1.0;
			long double guess = sign * binomial (k+1, k-n) ;
			printf ("its k=%d n=%d sum=%Lg  guess =%Lg \n", k,n,acc, guess);
		}
	}
#endif

#ifdef H_SQUARED
	for (n=0; n<10; n++)
	{
		for (m=0; m<10; m++)
		{
			long double acc= zetam1 (n+m+2);
			acc -= pow (0.5,n+m+2);
			acc *= binomial (n+m+1, n+1);
			
			long double h2 = 0.0L;
			for (k=0; k<50; k++)
			{
				h2 += ache_mp(m,k) * ache_mp(k,n);
			}

			long double sign = 1.0L;
			long double rem=0.0L;
			for (k=0; k<50; k++)
			{
				rem += sign * binomial (m+k+1,m)*zetam1(m+k+2) * ache_mp(k,n);
				sign = -sign;
			}

			printf ("its m=%d n=%d sum=%Lg  h2=%Lg    r=%Lg\n", m,n,acc, h2, acc+rem);
		}
		printf ("\n");
	}
#endif

#if EULER_GAMMA_SUB_N
	n=1;
	
	for (m=10; m<300000; m+=10000)
	{
		long double acc = 0.0L;
		long double term = 0.0L;
		for (k=1; k<= m; k++)
		{
			term = powl (logl(k),n) / k;
			acc += term;
		}
		long double bal = powl (logl(m), n+1) / (n+1);
		acc -= bal;

		printf ("its %d  acc=%Lg ball/acc=%Lg term=%Lg\n", m, acc, bal/acc, term);
	}
#endif

#ifdef Q_FUNC
	
	long double acc = 0.0L;
	long double term = 0.0L;
	long double z = 2.0;
	for (z=2.0; z<= 6.1; z+= 0.1)
	{
		// compute the q function
		acc = 0.0L;
		for (m=0; m<50; m++)
		{
			term = zbinomial (m+1, z-1.0) * zetam1(m+2);
			acc += term;
		}

		// compute finite differences
		long double af =0.0L;
		int nstart = floor(z);
		long double w = z-nstart;
		af = zetam1(nstart) + 1.0L;
		long double sn = -1.0L;
		for (m=1; m<70; m++)
		{
			long double sign = 1.0L;
			long double diff = 0.0L;
			for (k=0; k<=m; k++)
			{
				diff += sign * binomial (m,k) * zetam1(k+nstart);
				sign = -sign;
			}
			term = sn*diff * fbinomial (w,m);
			af += term;
			sn = -sn;
		}

		long double zeta = gsl_sf_zeta (z);
		long double diff = zeta - acc;
		printf ("z=%Lg, q=%Lg  az=%Lg zeta=%Lg diff=%Lg term=%Lg\n", z, acc, af, zeta, diff, term);
		// printf ("%Lg\t%Lg\t%Lg\t%Lg\n", z, acc, zeta, af);
	}
#endif

#ifdef VERIFY_BERNOULLI_MATRIX

	long double g_bern (int k, int m)
	{
		if (k==m) return powl (0.5L, k);
		if (k>m) return 0.0L;
		return binomial (m,k) * powl(0.5,m+1);
	}
	
	long double eig(int m, int p)
	{
		if (m>p) return 0.0;
		return binomial (p,m) * bernoulli (p-m);
	}

	for (p=0; p<20; p++)
	{
	for (k=0; k<=p; k++)
	{
		long double acc = 0.0;
		for (m=0; m<50; m++)
		{
			acc += g_bern (k,m)* eig(m,p);
		}
		if (fabs(eig (k,p)) >1.0e-16) {
			acc /= eig (k,p);
			printf ("p=%d k=%d eig = %Lg\n", p, k, acc);
		}
	}
	printf ("\n");
	}
#endif

#ifdef POCHHAMMER_CHECK
	long double s = 1.0;
	for (s=0.2; s<4.0; s+= 0.37) {
		long double z;
		for (z=0.0; z<1.0; z+= 0.1) {
		long double acc = 0.0;
		long double term;
		long double zp = 1.0;
		for (k=0; k<50; k++)
		{
			term = fbinomial (s-1, k) *zp;
			acc += term;
			zp *= z;
		}
		printf ("duude s=%Lg, z=%Lg pow=%Lg, sum=%Lg\n", s, z, powl(1.0+z, s-1), acc);
		}
	}
#endif

#ifdef MAJOR_BERNOULII_BOGOSITY
	long double h_bern (int m, int p)
	{
		if (p<m) return 0.0;
		if (p==m) return pow (0.5, p);
		double sign = 1.0L;
		if ((m+p)%2==1) sign = -1.0L;
		return pow(0.5, p+1)* fbinomial (p,m);
	}

	long double s = 1.0;
#if 0
	for (s=1.2; s<4.0; s+= 0.17) 
	{
		long double bb = 0.0L;
		long double sign = 1.0L;
		long double term;
		for (m=0; m<20; m++)
		{
			bb = 0.0;
			for (n=m+1; n<650; n++)
			{
				term= binomial (n,m) * fbinomial (s-1.0,n);
				sign = 1.0;
				if (n%2==1) sign = -1.0;
				bb += sign * term;
				printf ("s=%Lg n=%d term=%Lg summo=%Lg\n", s, n, sign*term, bb);
			}
			printf ("duude s=%Lg m=%d expect=%Lg bb=%Lg\n", s, m, fbinomial (s-1,m), bb);
		}
	}
#endif

	for (n=0; n<20; n++)
	{
		long double acc = 0.0L;
		for (m=0; m<=n; m++)
		{
			acc += binomial (n,m)/((m+1)*(m+2));
		}
		long double pred = pow(2,n+2) - 1 - (n+2);
		pred /= (n+1)*(n+2);
		printf ("its %d pred=%Lg acc=%Lg\n", n, pred, acc);
	}

	for (s=0.2; s<4.0; s+= 0.17) 
	{
		long double cb = 0.0L;
		long double sign = 1.0L;
		long double term;
		for (m=0; m<50; m++)
		{
			sign = 1.0L;
			for (n=0; n<50; n++)
			{
				term = fbinomial (s-1, n);
				term *= sign;
				term *= h_bern (m,n);
				term /= (m+1)*(m+2);
				cb += term;
				sign = -sign;
			}
		}
		
		long double qc = 0.0L;
		for (m=0; m<50; m++)
		{
			sign = 1.0L;
			if (m%2==1) sign = -1.0;

			term = fbinomial (s-1, m);
			term *= sign;
			term *= pow (0.5, m);
			term /= (m+1)*(m+2);
			qc += term;
			
			sign = -sign;
			long double part = 0.0;
			for (n=m+1; n<50; n++)
			{
				term = fbinomial (s-1, n);
				term *= sign;
				// term *= h_bern (m,n);
				term *= pow(0.5, n+1)* fbinomial (n,m);
				part += term;
				sign = -sign;
			}
			part /= (m+1)*(m+2);
			qc += part;
		}
		
		qc = 0.0;
		for (m=0; m<50; m++)
		{
			sign = 1.0L;
			if (m%2==1) sign = -1.0;

			term = fbinomial (s-1, m);
			term *= sign;
			term *= pow (0.5, m);
			term /= (m+1)*(m+2);
			qc += term;
			// printf ("duude m=%d term=%Lg sum=%Lg\n", m, term, qc);
		}
		printf (" ------------------------------------------------------------first=%Lg\n", qc);
			
		for (m=0; m<50; m++)
		{
			long double part = 0.0;
			for (n=m; n<50; n++)
			{
				term = fbinomial (s-1, n);
				sign = 1.0L;
				if (n%2==1) sign = -1.0;

				term *= sign;
				// term *= h_bern (m,n);
				term *= pow(0.5, n)* fbinomial (n,m);
				part += term;
			}
			part /= (m+1)*(m+2);
			qc += part;
		}
		qc *= 0.5;
		
		long double acc = 0.0;
		sign = 1.0L;
		for (m=0; m<50; m++)
		{
			term = fbinomial (s-1, m);
			term *= sign;
			term *= pow (0.5, m);
			term /=  (m+1)*(m+2);
			term *= (pow(2.0,m+2) - (m+2));
			acc += term;
			// printf ("anti-duude m=%d term=%Lg sum=%Lg\n", m, term, acc);
			sign = -sign;
		}
		acc *= 0.5;
		long double ch = (s-1)/(s+1) + pow(0.5,s);
		ch /= s;

		
		printf ("its %Lg char=%Lg pre-sum=%Lg %Lg  sum=%Lg \n", s, ch, cb, qc, acc);
	}
#endif

#ifdef MORE_ZETA
#define LEN 650
	long double lam[LEN];

	for (m=0; m<LEN; m++)
	{
		lam[m] = pow(0.5,m);
	}

	long double complex s = 0.0;
	for (p=0; p<50; p++)
	{
		s += 0.2I;
		long double complex acc = 0.0;
		long double complex term;
		long double sign = 1.0L;
		for (m=0; m<LEN; m++)
		{
			term = cbinomial (s-1, m);
			term *= sign;
			term *= lam[m];
			term /=  (m+1)*(m+2);
			term *= (pow(2.0,m+2) - (m+2));
			acc += term;
			// printf ("anti-duude m=%d term=%Lg sum=%Lg\n", m, term, acc);
			sign = -sign;
		}
		acc *= 0.5;

		acc -= 2.0/(s+1.0);
		acc *= -s;
		
		long double complex ch = (s-1)/(s+1) + cpow(0.5,s);
		ch /= s;
		
		ch -= 2.0/(s+1.0);
		ch *= -s;
		
		
		printf ("s=%Lg +I %Lg    acc=%Lg +I %Lg   alt=%Lg +I %Lg \n", creall(s), cimagl(s), creall(acc), cimagl(acc),
							 creall(ch), cimagl(ch));
	}
#endif

#if ONE_WAY
	int j;
	// for (k=0; k<10; k++)
	// {
	for (j=0; j<10; j++)
	{
		long double duh=0.0;
		// for (n=k; n<40; n++)
		for (n=0; n<10; n++)
		{
			long double acc = 0.0;
			for (m=0; m<=j; m++)
			{
				acc += ache_mp(n,m) * stirling_first (j,m);
			}
			long double stir = stirling_second (n,k);
			long double term = acc * stir;
			duh += term;
			// printf ("duude j=%d k=%d n=%d acc=%Lg stir=%Lg term=%Lg duhh=%Lg\n", j, k, n, acc, stir, term, duh);
			printf ("duude j=%d n=%d acc=%Lg\n", j, n, acc);
		}
		printf ("---\n");
	// }	
	}
#endif

#if REV
	int j;
	for (j=0; j<8; j++)
	{
	for (k=0; k<8; k++)
	{
		long double duh=0.0;

		for (m=0; m<=j; m++)
		{
			long double acc = 0.0;
			for (n=k; n<40; n++)
			{
				acc += ache_mp(n,m) * stirling_second (n,k);
			}
			long double stir = stirling_first (j,m);
			long double term = acc * stir;
			duh += term;
			printf ("duude j=%d k=%d m=%d acc=%Lg stir=%Lg term=%Lg duhh=%Lg\n", j,k, m, acc, stir, term, duh);
		}
		printf ("---\n");
	}	
	}
#endif 

#ifdef STIR 
	for (n=0; n<30; n++)
	{
		int j;
		for (j=0; j<=n; j++)
		{
			long double ha = 1.0L;
			long double acc = 0.0;
			for (k=0; k<=n-j; k++)
			{
				long double term = stirling_first (n,j+k);
				term *= binomial (j+k,j);
				acc += ha*term;
				// printf ("its n=%d j=%d k=%d term=%Lg acc=%Lg\n", n,j,k,term,acc);
				ha *= -0.5;
			}
			acc /= factorial (n);
			printf ("---------- its n=%d j=%d acc=%Lg\n", n,j, acc);
		}
	}
#endif

#ifdef ZXOOOO
long double zoot (int s)
{
	int n;
	long double sign = 1.0L;
	long double acc=0.0L;
	for (n=0; n<s; n++)
	{
		acc += sign *binomial (s-1,n)*a_sub_n(n);
		sign = -sign;
	}
	acc *= -s;
	acc -= 0.5;
	acc += s / ((long double)(s-1));
	return acc;
}

	for (k=2; k<30; k++)
	{
		long double ch = zoot(k) - 1.0L - zetam1(k);
		printf ("its %d %Lg\n", k, ch);
	}
#endif

#if 0
	int i,j,k;
	i = 1;
	k = 2;
	for (j=0; j<10; j++)
	{
		double dij = ache_mp(i,j) / ache_mp(j,i);
		double djk = ache_mp(j,k) / ache_mp(k,j);
		printf ("Its %d %g \n", j, dij*djk);
	}
#endif

   double a=zetam1(2);
	double b=zetam1(2)-zetam1(3);
	double c=2.0*zetam1(3);
	double d=c-3.0*zetam1(4);
	printf ("a=%g  b=%g  c=%g  d=%g\n", a,b,c,d);

	double apd = zetam1(2)+2.0*zetam1(3)-3.0*zetam1(4);
	double amd = zetam1(2)-2.0*zetam1(3)+3.0*zetam1(4);
	double bc = 2.0*(zetam1(2) - zetam1(3))*zetam1(3);

	double lp=0.5*apd + 0.5*sqrt(amd*amd+4.0*bc);
	double lm=0.5*apd - 0.5*sqrt(amd*amd+4.0*bc);
  printf ("duude a+d=%g\n",apd);	
  printf ("duude a-d=%g\n",amd);	
  printf ("duude bc=%g\n",bc);	
  printf ("lp=%g    lm=%g\n", lp, lm);

  double cp = -(zetam1(2)-lp)/(zetam1(2)-zetam1(3));
  double cm = -(zetam1(2)-lm)/(zetam1(2)-zetam1(3));

  printf ("cp=%g   cm=%g\n", cp, cm);
  
  double csp = -0.5*(zetam1(2)-lp)/zetam1(3);
  double csm = -0.5*(zetam1(2)-lm)/zetam1(3);

  printf ("csp=%g   csm=%g\n", csp, csm);
  
  double np = 1.0/(1.0-0.5*(a-lp)*(1.0/b+1.0/c)+(a-lp)*(a-lp)/(3.0*b*c));
  double nm = 1.0/(1.0-0.5*(a-lm)*(1.0/b+1.0/c)+(a-lm)*(a-lm)/(3.0*b*c));

  printf ("np=%g   nm=%g\n", np, nm);

  printf ("nlp=%g   nlm=%g\n", np*lp, nm*lm);

  double A = 0.5*(np*lp+nm*lm);
  double B = (np*lp*csp+nm*lm*csm)/6.0;
  double C = (np*lp*cp+nm*lm*cm)/2.0;
  double D = (np*lp*cp*csp+nm*lm*cm*csm)/6.0;

  printf ("A=%g  B=%g  C=%g  D=%g\n", A,B,C,D);

  double P = 1.0 - A - B;
  double Q = 1.0 - C - D;
  double R = A+B+C+D;

  printf ("P=%g  Q=%g  R=%g\n", P,Q, R);

  double re= -0.5* Q/P;
  printf ("re=%g\n", re);
  double disc = Q*Q - 4.0*P*R;

  printf ("disc=%g\n", disc);
  double im=0.5*sqrt(-disc)/P;
  printf ("im=%g\n", im);


#if 0
	int i,j;
	for (i=0; i<10; i++)
	{
		for (j=0; j<7; j++) 
		{
			printf ("%8.6Lf\t", ache_mp(i,j));
		}
		printf ("\n");
	}
#endif
	

	return 0;
}
