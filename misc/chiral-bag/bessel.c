/*******************************************************************
 * FILE: bessel.c
 * FUNCTION: Compute spherical bessel functions
 * HISTORY: Created by Linas Vepstas <linas@linas.org> in 1983
 * c	1983 VAX/VMS 11/780 at SUNYSBNP
 * C	1987- Prime os 19.1  Sacaly, France.  The Goddamn Compiler
 * C		doesn't like for-- end for structures.  Revert to 
 * c		using line numbers.
 * C	1987- IBM PC -- MSDOS (Saclay, France; Krakow, Poland)
 * c	1988 vax/vmx 11/780 (visitng StonyBrook)
 * C	1989- IBM RT/ AIX 2.2.1   -- F77 doesn't like underscores
 * July 2003 -- port from FORTRAN to C language
******************************************************************/

#include <math.h>
#include <stdio.h>
#include "bessel.h"

/* ================================================================== */
/*
        SUBROUTINE BesselRecurrenceStartRegion (Z, NORDER, NSTART)
C       GIVEN AN ARGUMENT X AND ORDER NORDER, THE PARENT ROUTINE USES A
C       BACKWARDS RECURRENCE RELATION TO COMPUTE BESSEL FUNCTIONS.
C       IN ORDER TO OBTAIN AND ACCURATE ANSWER, THE RECURRENCE MUST BE
C       STARTED AT A HIGH ENOUGH ORDER SUCH THAT SUCESSIVE TERMS BECOME
C       BETTER AND BETTER APPROXIMATIONS TO THE TRUE RECURRENCE RELATION
C       J(N+1)(Z) + J(N-1)(Z) = ((2N+1)/Z) * JN(Z)
C       I.E. ONE NEEDS TO START IN THE "ASYMPTOTIC REGION" OF SMALL X
C       AND LARGE ORDER, WHERE THE LARGE ORDER BESSEL FUN. HAS NOT YET
C       ENCOUNTERED ITS FIRST ZERO. OUR VERY CRUDE ESTIMATE FOR THIS
C       REGION IS NSTART = MAX (2*X, N+10)
C       (ABROMOWIZ + STEGUN, PAGE 452FF)

C       THIS ROUTINE ACCEPTS THE VALUES X AND NORDER, AND RETURNS THE
C       DESIRED STARTING POINT NSTART.

C       NOTES:
C       1) THE BESSEL ROUTINE WILL NOT WORK RELIABLY FOR ARGUMENTS
C               LESS THAN 2.0E-6 OR SO
C       2) THE IRREGULAR BESSEL FUNCTIONS WILL OVERFLOW (1.0E-30 OR SO)
C               FOR EXCESSIVELY SMALL ARGUMENTS (LIKEWISE, THE REGULAR
C               FUNCTIONS WILL UNDERFLOW FOR THE SAME VALUES).  THEREFORE,
C               A LOWER BOUND ON THE ARGUMENT HAS BEEN FIXED, DETERMINED
C               IN PART BY THE SMALL ARGUMENT EXPANSIONS FOR THE BESSEL
C               FUNCTIONS (REF: ABROMOWITZ + STEGUN). THE ARGUMENT RETURNED
C               BY THIS ROUTINE HAS BEEN BOUNDED BELOW BY THIS "SMALL"
C               BOUNDING ARGUMENT (WHICH MAY BE QUITE LARGE, FOR LARGE ORDERS)
*/

#define fmax(x,y) (((x)>(y))?(x):(y))

static int 
BesselRecurrenceStartRegion (double z, int norder)
{
	int nstart, iaz;
	double az, an, zmin, peritr;

	az = fabs(z);
	an = (double)(norder + 1);
	zmin = (-60.0 + an * (log(an)-1.0)) / an;
	zmin = 2.0 * exp (zmin + 12.5 / (an*an)) + 1.0e-6;
	az = fmax (zmin, az);
	z = copysign (az, z);

	peritr = 2.0 * an / az;
	if (3.0 < peritr)
	{
		double extra = 20.0 / log (peritr);
		nstart = norder + 8 + (int) (extra);
	}
	else 
	{
		iaz = (int) az;
		if (33 <= iaz) { nstart = 2 * iaz; }
		else if (13 <= iaz) { nstart = 3 * iaz; }
		else if (7 <= iaz) { nstart = 4 * iaz; }
		else { nstart = 6 * iaz +20; }
	}
	return nstart;
}

/* ================================================================== */
/*
C       THIS SUBROUTINE CALCULATES THE SPHERICAL BESSEL FUNCTION FOR
C       ALL ORDERS .LE. N FOR (POSITIVE) ARGUMENT Z
C       THE RESULTS ARE RETURNED IN THE ARRAY BESS (0:N)
*/
void 
bessel (int n, double x, double *bess)
{
#define ARRSZ 3650
	double tmp [ARRSZ];
	int i, j, m, nstart, mon[ARRSZ];
	double s, c, r, z, oneoz, scale;

/*
C       WE USE BACKWARD RECURENCE TO FIND THE BESSEL FUNCTIONS
C       J(N+1)(Z) + J(N-1)(Z) = ((2N+1)/Z) * JN(Z)
C       (ABROMOWIZ + STEGUN, PAGE 452FF)
C       THE RECURRENCE RELATION MUST BE STARTED AT A HIGH ORDER.
C       SEE COMMENT CARDS IN THE SUBROUTINE "BesselRecurrenceStartRegion"
*/
	z = x;
	nstart = BesselRecurrenceStartRegion (z,n);
	if (ARRSZ <= nstart)
	{
		printf ("BESSEL TEMPORARY ARRAY OUT OF BOUNDS\n");
		return;
	} 
	tmp[nstart] = 0.0;
	tmp[nstart-1] = 1.0E-300;
	mon[nstart-1] = 0;
	oneoz = 1.0/z;
	for (i=nstart-1; i>0; i--)
	{
		tmp [i-1] = oneoz * ((double) (2*i+1)) * tmp[i] - tmp [i+1];
		mon [i-1] = mon [i];
		if (fabs (tmp [i-1]) > 1.0E+250)
		{
			tmp [i] = tmp [i] * 1.0E-500;
			tmp [i-1] = tmp [i-1] * 1.0E-500;
			mon [i] = mon [i+1] + 500;
			mon [i-1] = mon [i];
		}
	}

	/* Compute r = sin (z) * oneoz / tmp[0]; from first principles */
	s = z * tmp[0];
	c = tmp[0] - z * tmp[1];

	i = 0;
	while (c > 1.0e100) 
	{
		i++;
		c *= 1.0e-100;
		s *= 1.0e-100;
	}
	r = 1.0 / sqrt (c*c + s*s);
	while (0 < i)
	{
		r *= 1.0e-100;
		i--;
	}
	/* compute sign */
	r = copysign (r, s);
	s = z / M_PI;
	i = (int) floor(s);
	if (i%2 != 0) { r = -r; }

/*
C       RENORMALIZE THE VALUES.  SINCE EXPONENTIATION TAKES ABOUT A FACTOR
C       OF 20 MORE CPU TIME THAN IF STATEMENTS, OR MULTIPLICATION,
C       I'VE RIGGED A LITTLE SCHEME HERE WHEN ALL I REALLY WANT TO DO IS
C       BESS (I) = C * TMP(I) * (10.0D0 ** (MON(I)-M))
*/
   m = mon [0];
	scale = r;
	bess [0] = scale * tmp[0];
	for (i=1; i<=n; i++)
	{
		if (mon[i] != mon[i-1]) 
		{ 
			double p = 1.0;
			/* scale = r * pow (10.0, (double) (mon[i]-m)); */
			for (j=0; j<mon[i]-m; j+=500) p *=1.0e500;
			scale = r * p;
		}
		bess[i] = scale * tmp[i];
	}
}

/* ============================================================== */
/*
        SUBROUTINE QUICKBESSEL (N, X, BESSM1, BESS, BESSP1)

C       THIS SUBROUTINE CALCULATES THE
C       REGULAR SPHERICAL BESSEL FUNCTIONS BESS AND THE
C       IRREGULAR SPHERICAL BESSEL FUNCTIONS NEU
C       FOR THE ORDERS N-1, N AND N+1, FOR (POSITIVE) ARGUMENT Z.
* Bugs: not as fast as it could be ...
*/

void 
quickbessel (int n, double x, 
             double *bessm1, double *bess, double *bessp1)
{
	double val [ARRSZ];
	bessel (n,x,val);

	*bessm1 = val[n-1];
	*bess = val[n];
	*bessp1 = val[n+1];
}

#if 0
/***
C       ******************************************************************

        SUBROUTINE BESSNEU (N, X, BREGULAR, BIRREG)

C       THIS SUBROUTINE CALCULATES THE SPHERICAL BESSEL FUNCTIONS,
C       BOTH THE ONES REGULAR AT THE ORIGEN, AND THE IRREGULARS (NEUMANN)
C       FOR ALL ORDERS .LE. N FOR (POSITIVE) ARGUMENT Z
C       THE RESULTS ARE RETURNED IN THE ARRAYS BREGULAR (0:N), BIRREG (0:N)
c	Bugs fixed:
c	May 4 1988
c	renomalization statement was  IF (MON (I) .NE. MON (I-1))
c	changed to IF (MON (I) .NE. MON (I+1))

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DOUBLE PRECISION BREGULAR (0:N), BIRREG (0:N)
        DOUBLE PRECISION TMP (-200:400)
        INTEGER * 4 MON (-200:400)

C       WE USE BACKWARD RECURENCE TO FIND THE BESSEL FUNCTIONS
C       J(N+1)(Z) + J(N-1)(Z) = ((2N+1)/Z) * JN(Z)
C       (ABROMOWIZ + STEGUN, PAGE 452FF)
C       THE RECURRENCE RELATION MUST BE STARTED AT A HIGH ORDER.
C       SEE COMMENT CARDS IN THE SUBROUTINE 'BesselRecurrenceStartRegion'

        Z = X
        CALL BesselRecurrenceStartRegion (Z, N, NSTART)
        IF (NSTART .GT. 400) STOP 'BESSEL TEMPORARY ARRAY OUT OF BOUNDS'
        TMP (NSTART) = 0.0
        TMP (NSTART-1) = 1.0D-30
        MON (NSTART-1) = 0
        ONEOZ = 1.0D0 / Z
C       WE'VE GOT TO GO TO NEGATIVE ORDERS TO GET THE IRREGULAR SOLUTIONS
        DO 1001 I = NSTART-1, -(N+2), -1
                TMP (I-1) = ONEOZ * DBLE (2*I+1) * TMP(I) - TMP (I+1)
                MON (I-1) = MON (I)
C       PREVENT OVERFLOW BY KEEPING TRACK OF EXPONENT MANUALLY
                IF (DABS (TMP (I-1)) .GT. 1.0D+20) THEN
                        TMP (I) = TMP (I) * 1.0D-30
                        TMP (I-1) = TMP (I-1) * 1.0D-30
                        MON (I) = MON (I+1) + 30
                        MON (I-1) = MON (I)
                END IF
1001    CONTINUE
        C = DSIN (Z) * ONEOZ
        C = C / TMP(0)
        M = MON (0)
C       RENORMALIZE THE VALUES.  SINCE EXPONENTIATION TAKES ABOUT A FACTOR
C       OF 20 MORE CPU TIME THAN IF STATEMENTS, OR MULTIPLICATION,
C       I'VE RIGGED A LITTLE SCHEME HERE WHEN ALL I REALLY WANT TO DO IS
C       BESS (I) = C * TMP(I) * (10.0D0 ** (MON(I)-M))
        SCALE = C
        BREGULAR (0) = SCALE * TMP(0)
        DO 1002 I = 1, N
                IF (MON (I) .NE. MON (I-1))
     +           SCALE = C * (10.0D0 ** (MON(I)-M))
                BREGULAR (I) = SCALE * TMP(I)
1002    CONTINUE

C       RESCALE BESSEL FUNCTTIONS IRREGULAR AT THE ORIGEN
        SCALE = C
        DO 1003 I = -1, -(N+1), -1
c	bug May 4, 1988
                IF (MON (I) .NE. MON (I+1))
     +           SCALE = C * (10.0D0 ** (MON(I)-M))
                BIRREG (-(I+1)) = SCALE * TMP(I)
1003    CONTINUE
C       AND THIER SIGN
        DO 1004 I = 0,N,2
                BIRREG (I) = - BIRREG (I)
1004    CONTINUE
        RETURN
        END

C       **************************************************************

        SUBROUTINE QUICKBESSEL (N, X, BESSM1, BESS, BESSP1)

C       THIS SUBROUTINE CALCULATES THE
C       REGULAR SPHERICAL BESSEL FUNCTIONS BESS AND THE
C       IRREGULAR SPHERICAL BESSEL FUNCTIONS NEU
C       FOR THE ORDERS N-1, N AND N+1, FOR (POSITIVE) ARGUMENT Z.

C       NOTE:****
C               FOR N=0, THE BESSEL FUNS FOR N=-1 ARE CONSIDERED TO BE
C               UNDEFINED, AND THE VALUE OF 0.0 IS RETURNED FOR THE FUNCTIONS.
C               THIS IS A PROGRAMING CONSIDERATION **ONLY**, AS BESSEL FUNS.
C               FOR NEGATIVE ARGUMENTS REALLY ARE DEFINED AND ARE
C               RELATED TO THE IRREGULAR SOLUTIONS

C       THE RESULTS ARE RETURNED AS SHOWN ABOVE

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DOUBLE PRECISION TMP (-2:2000)
        INTEGER * 4 MON (-2:2000)

C       WE USE BACKWARD RECURENCE TO FIND THE BESSEL FUNCTIONS
C       J(N+1)(Z) + J(N-1)(Z) = ((2N+1)/Z) * JN(Z)
C       (ABROMOWIZ + STEGUN, PAGE 452FF)
C       THE RECURRENCE RELATION MUST BE STARTED AT A HIGH ORDER.
C       SEE COMMENT CARDS IN THE SUBROUTINE 'BesselRecurrenceStartRegion'

        Z = X
        CALL BesselRecurrenceStartRegion (Z, N, NSTART)
        IF (NSTART .GT. 2000) STOP 'BESSEL WORK ARRAY OUT OF BOUNDS'
        TMP (NSTART) = 0.0
        TMP (NSTART-1) = 1.0D-28
        MON (NSTART-1) = 0
        ONEOZ = 1.0D0 / Z
        DO 1001 I = NSTART-1, 1, -1
                TMP (I-1) = ONEOZ * DBLE (2*I+1) * TMP(I) - TMP (I+1)
                MON (I-1) = MON (I)
                IF (DABS (TMP (I-1)) .GT. 1.0D+15) THEN
                        TMP (I) = TMP (I) * 1.0D-30
                        TMP (I-1) = TMP (I-1) * 1.0D-30
                        MON (I) = MON (I+1) + 30
                        MON (I-1) = MON (I)
                END IF
1001    CONTINUE
        C = DSIN (Z) * ONEOZ
        C = C/ TMP(0)
        M = MON (0)
C       N=0 SPECIAL CASE CONSTANTS
        MON (-1) = M
        MON (-2) = M
        TMP (-1) = 0.0D0
C       RENORMALIZE THE VALUES.  SINCE EXPONENTIATION TAKES ABOUT A FACTOR
C       OF 20 MORE CPU TIME THAN IF STATEMENTS, OR MULTIPLICATION,
C       I'VE RIGGED A LITTLE SCHEME HERE WHEN ALL I REALLY WANT TO DO IS
C       BESS (I) = C * TMP(I) * (10.0D0 ** (MON(I)-M))
        SCALE = C * (10.0D0 ** (MON (N-2) - M))
        DO 1002 I = N-1, N+1
                IF (MON (I) .NE. MON (I-1))
     +           SCALE = C * (10.0D0 ** (MON(I)-M))
                TMP(I) = SCALE * TMP(I)
1002    CONTINUE
        BESSM1 = TMP (N-1)
        BESS = TMP (N)
        BESSP1 = TMP (N+1)
        RETURN
        END

C       **************************************************************


        SUBROUTINE QUICKBESSNEU
     +          (N,X, BESSM1, BESS, BESSP1, NEUM1, NEU, NEUP1)

C       THIS SUBROUTINE CALCULATES THE
C       REGULAR SPHERICAL BESSEL FUNCTIONS BESS AND THE
C       IRREGULAR SPHERICAL BESSEL FUNCTIONS NEU
C       FOR THE ORDERS N-1, N AND N+1, FOR (POSITIVE) ARGUMENT Z.

C       NOTE:****
C               FOR N=0, THE BESSEL FUNS FOR N=-1 ARE CONSIDERED TO BE
C               UNDEFINED, AND THE VALUE OF 0.0 IS RETURNED FOR THE FUNCTIONS.
C               THIS IS A PROGRAMING CONSIDERATION **ONLY**, AS BESSEL FUNS.
C               FOR NEGATIVE ARGUMENTS REALLY ARE DEFINED AND ARE
C               RELATED TO THE IRREGULAR SOLUTIONS

C       THE RESULTS ARE RETURNED AS SHOWN ABOVE

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DOUBLE PRECISION TMP (-200:400)
        DOUBLE PRECISION NEUM1, NEU, NEUP1
        INTEGER * 4 MON (-200:400)

C       WE USE BACKWARD RECURENCE TO FIND THE BESSEL FUNCTIONS
C       J(N+1)(Z) + J(N-1)(Z) = ((2N+1)/Z) * JN(Z)
C       (ABROMOWIZ + STEGUN, PAGE 452FF)
C       THE RECURRENCE RELATION MUST BE STARTED AT A HIGH ORDER.
C       SEE COMMENT CARDS IN THE SUBROUTINE 'BesselRecurrenceStartRegion'

        Z = X
        CALL BesselRecurrenceStartRegion (Z, N, NSTART)
        IF (NSTART .GT. 400) STOP 'BESSEL TEMPORARY ARRAY OUT OF BOUNDS'

        TMP (NSTART) = 0.0
        TMP (NSTART-1) = 1.0D-30
        MON (NSTART-1) = 0
        ONEOZ = 1.0D0 / Z
C       WE'VE GOT TO GO TO NEGATIVE ORDERS TO GET THE IRREGULAR SOLUTIONS
        DO 1001 I = NSTART-1, -(N+2), -1
                TMP (I-1) = ONEOZ * DBLE (2*I+1) * TMP(I) - TMP (I+1)
                MON (I-1) = MON (I)
C       PREVENT OVERFLOW BY KEEPING TRACK OF EXPONENT MANUALLY
                IF (DABS (TMP (I-1)) .GT. 1.0D+20) THEN
                        TMP (I) = TMP (I) * 1.0D-30
                        TMP (I-1) = TMP (I-1) * 1.0D-30
                        MON (I) = MON (I+1) + 30
                        MON (I-1) = MON (I)
                END IF
1001    CONTINUE
        C = DSIN (Z) * ONEOZ
        C = C / TMP(0)
        M = MON (0)
C       RENORMALIZE THE VALUES.  SINCE EXPONENTIATION TAKES ABOUT A FACTOR
C       OF 20 MORE CPU TIME THAN IF STATEMENTS, OR MULTIPLICATION,
C       I'VE RIGGED A LITTLE SCHEME HERE WHEN ALL I REALLY WANT TO DO IS
C       BESS (I) = C * TMP(I) * (10.0D0 ** (MON(I)-M))

C       RENORMALIZE THE POSITIVE ORDERS
C       FOR N=0, THEN BESS (-1) IS SET TO ZERO
C       BUT WE WANT TO KEEP THIS VALUE FOR THE NEGATIVE ORDERS, SO STORE IT
        TEMP = TMP (-1)
        TMP (-1) = 0.0D0
        SCALE = C * (10.0D0 ** (MON (N-2) - M))
        IF (MON (N-1) .NE. MON (N-2))
     +   SCALE = C * (10.0D0 ** (MON(N-1)-M))
        BESSM1 = SCALE * TMP (N-1)
        IF (MON (N) .NE. MON (N-1))
     +   SCALE = C * (10.0D0 ** (MON(N)-M))
        BESS = SCALE * TMP (N)
        IF (MON (N+1) .NE. MON (N))
     +   SCALE = C * (10.0D0 ** (MON(N+1)-M))
        BESSP1 = SCALE * TMP (N+1)
        TMP (-1) = TEMP
C       RENORMALIZE THE NEGATIVE ORDERS
C       THESE HAVE AN EXTRA FACTOR OF (-)**N IN THEM:: THIS IS SING
C       SING IS NEGATIVE FOR N EVEN
        SING = DBLE (2*MOD(N,2)-1)
C       FOR N=0, THEN NEU (-1) IS SET TO ZERO
        TMP (0) = 0.0D0
        SCALE = C * (10.0D0 ** (MON (-N+1) - M))
        IF (MON (-N) .NE. MON (-N+1))
     +   SCALE = C * (10.0D0 ** (MON(-N)-M))
        NEUM1 = - SING * SCALE * TMP(-N)
        IF (MON (-N-1) .NE. MON (-N))
     +   SCALE = C * (10.0D0 ** (MON(-N-1)-M))
        NEU = SING * SCALE * TMP(-N-1)
        IF (MON (-N-2) .NE. MON (-N-1))
     +   SCALE = C * (10.0D0 ** (MON(-N-2)-M))
        NEUP1 = - SING * SCALE * TMP(-N-2)
        RETURN
        END

C       ************************************************************
**/
#endif
