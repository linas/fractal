/*
 * FILE: bag_ener.c
 * FUNCTION: compute energy levels of firee quarks in chiral bag
 * HISTORY: Created by Linas Vepstas <linas@linas.org> 1983
 *   July 2003 ported from FORTRAN to C
 */

#include <math.h>
#include <stdio.h>
#include "bessel.h"
#include "find_zero.h"

/* ============================================================= */
/* Implicit equations defining the quark energy levels 
 * energy levels occur at the zero's of these eqn's.
 * So we need to search for zero's.
 */

static double efn1_ptyk;
static double efn1_co;
static double efn1_si;
static int efn1_kay;

static void efninit (double theta, int ispect, int kay, double emax, int kpty)
{
	efn1_ptyk = kpty;
	efn1_kay = kay;
	efn1_co = cos (theta);
	efn1_si = sin (theta);
}

/* ============================================================= */

static double efn (double omega)
{
	double bkm1, bk, bkp1;
	double val, bks;

	if (0 < efn1_kay)
	{
		quickbessel (efn1_kay, omega, &bkm1, &bk, &bkp1);
		bks = bk*bk;
		val = efn1_co * (bkp1 * bkm1 - bks) 
			+ efn1_ptyk * bk * (bkp1 - bkm1)
			+ efn1_si * bks / omega;
	}
	else
	{
		double sn, cs, b0, b1;
		sn = sin (omega);
		cs = cos (omega);
		b0 = sn/omega;
		b1 = (b0 - cs) / omega;
		val = efn1_co * b1 - (efn1_ptyk - efn1_si) * b0;
	}
	return val;
}

/* ============================================================= */

#if 0
static double yefn (double omega)
{
	double bkm1, bk, bkp1;
	double ykm1, yk, ykp1;
	double val, yks;

	if (0 < efn1_kay)
	{
		quickbessneu (efn1_kay, omega, &bkm1, &bk, &bkp1, &ykm1, &yk, &ykp1);

		yks = yk*yk;
		val = efn1_co * (ykp1 * ykm1-yks) 
		     + efn1_ptyk * yk * (ykp1-ykm1)
		     + efn1_si * yks / omega;
	}
	else
	{
		double sn, cs, y0, y1;
		sn = sin (omega);
		cs = cos (omega);

		y0 = -cs/omega;
		y1 = (y0 - sn)/omega;
		val = efn1_co * y1 - (efn1_ptyk - efn1_si) * y0;
	}
	return val;
}
#endif

/* ============================================================= */
/*
        SUBROUTINE ENERGY (theta, ispect, kay, emax, kpty)

C       THIS SUBROUTINE RETURNS zeros of an ENERGY EIGENVALUE equation,
c       with parameters theta, ispect, kay, emax

C       ******* INPUT *******

C       THETA- CHIRAL ANGLE
C       KPTY-   PARAMETER INDICATING PARITY
C       KAY-    =L+S+T    K MUST BE INTEGER.
C       ISPECT- INTEGER, +1 OR -1, IS THE SIGN OF THE SPECTRUM
C               I.E. WE LOOK FOR POS. OR NEG. ENERGY SEA SOLUTIONS

C       ******* OUTPUT ******

C       ENER-   NENL zeros were found below EMAX. these are returned in common

C               WHEN WORKING WITH K = L+S+T, WE TAKE
C               1-- EQUATION FOR K=L            KPTY = PTYK = +1
C               2-- EQUATION FOR K=L+/-1        KPTY = PTYK = -1

C               WHEN WORKING WITH J = L+S, WE TAKE
C               1-- EQUATION FOR J = L+1/2      KPTY = PTYK = +1
C               2-- EQUATION FOR J = L-1/2      KPTY = PTYK = -1
*/

int
quark_energy (double *ener, double theta, int ispect, int kay, 
               double emax,  int kpty)
{
/*
C       MAXIT-  MAXIMUM NUMBER OF ITERATIONS THE CONVERGENCE ALGORITHM
C       IS TO PERFORM WHEN SEARCHING FOR THE EIGENVALUES.  IN PRACTICE,
C       IT SEEMS TO TAKE LESS THAN TEN ITERATIONS-- SUGGEST MAXIT =15
C       if there is failure, an error message will be printed
C       takes about 40 for K=1000, moer for larger K's
*/
	int maxit = 100;
/*
C       NSIG-   NUMBER OF SIGNIFICANT DIGITS FOR THE SEARCH ROUTINE,
C       I.E. THAT THE ENERGIES ARE TO BE CALCULATED TO.
c       NSIG = 15   -- this will work on the Vax and the PC but ...
*/
	int nsig = 14;

/*
C       LOOK FOR ZEROS of efn IN THE INTERVAL A .LE. ZERO .LE. B

C       WE WILL NOT SEARCH FOR ZEROS NEAR ZERO, EXCEPT FOR K=0.
C       WHY? FOR LARGE KAY, THERE ARE DEFINETLY NO
C       ZEROS NEAR ZERO.  HOWEVER, FOR LARGE KAY, EFN (OMEGA) IS
C       *VERY* *SMALL* NEAR OMEGA = ZERO, AND, IF ONE IS NOT CAREFUL,
C       THE ALGORITHM WILL RETURN SPURIOUS ZEROS. THUS, WE START THE
C       SEARCH AWAY FROM ZERO.

C       THERE ARE FURTHER PITFALLS:  AS ENERGY INCREASES, ONE TENDS
C       TO FIND PAIRS OF EIGENVALUES THAT ARE VERY CLOSE TO EACH OTHER.
C       IF ONE IS NOT CAREFUL, THEN ONE WILL MISS BOTH. AS A CURE,
C       ONE MAY ATTEMPT TO USE A SMALL INCREMENT IN THE SEARCH FOR
C       ZEROS, OR, ALTERNATLY, DECREASE THE INTERVAL AS THE SEARCH
C       GOES ON. I'VE CHOOSEN A DECREASING STEP SIZE.  I HOPE IT
C       DECREASES SUFICIENTLY RAPIDLY TO CATCH ALL THE ZEROS.

C       GOT IT: FOR LARGE ENERGIES, ENERGY LEVELS OCCUR IN APPROXIMATELY
C       DEGENERATE PAIRS.  THE SPLITTING BETWEEN PAIRS OCCURS AT INTERVALS
C       OF PI + O(1/ENERGY). THE INTRA-PAIR SPLITTING GOES LIKE
C       DELTA ENERGY = 1/( (PI/2) * N ) + O(1/N**2)
C       WHERE N DENOTES THE NTH ENERGY LEVEL.  THEREFORE, WE NEED THE STEP
C       SIZE TO BE .LT. DELTA ENERGY IF WE ARE TO CATCH BOTH MEMEBERS OF
C       THE PAIR.
*/

	double a,b, bb,sygnus, step;
	double funa, funb;
	int numfound, justfound;

	/* initialize the equation whos parameters are to be found */
	efninit (theta, ispect, kay, emax, kpty);

	a = 0.0;
	sygnus = ispect;
	b = fmax ((double)(kay-2), 3.2e-6);
	step = 0.1;
	funb = efn(sygnus * b);
	numfound = 0;
	justfound = -99;

	/* CHECK TO SEE IF WE WANT TO FIND MORE ENERGIES */
	while (a < emax)
	{
/*
C	       DEFINE THE NEXT INTERVAL
c	       if we have not recently found a zero, step along with
c	       step size step. but if we have found a zero recently,
c	       take a big step and save some cpu time.
*/
		if (0 > justfound)
		{
			a = b;
			b = b + step;
			funa = funb;
			funb = efn (sygnus * b);
		}
		else
		{
			a = fabs (ener [numfound-1]) + step;
			if (20 > numfound)
			{
				b = a + 2.2;
			}
			else
			{
				b = a + 2.95;
			}
			funa = efn (sygnus * a);
			funb = efn (sygnus * b);
			justfound = -99;
		}

		if (0.0 > funa*funb)
		{
			/* FOUND A ZERO IN INTERVAL.  CALL CONVERGENCE ALGORITHM. */
			bb = FindZero (efn, 0.0, nsig, sygnus*a, sygnus*b, maxit);

			/* STORE AND PRINT THE ENERGIES FOUND */
			ener [numfound] = bb;
			numfound++;

			/*  ***** DECREASE step SIZE PER COMMENTS ABOVE */
			justfound = 99;
			if (1 < numfound)
			{
				int nnn = numfound-1;
				double dil = ener[nnn] - ener [nnn-1];
				dil = 0.4 * fabs (dil);
				step = fmin (step, dil);
			}
		}
	}
	return numfound;
}

/* ============================================================= */
