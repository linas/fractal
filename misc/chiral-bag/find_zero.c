/*
 * FILE: find_zero.c
 * FUNCTION: find zero of a function
 * HISTORY:
C   IMSL ROUTINE NAME   - ZBRENT
c	Vax/VMS 11/780
c      1984  heavily modified by L. Vepstas because the imsl routine was
c       just an absolute load of crap
c	
c	The blinking Pr1me Fortran won't shit take do-end do loops!!!!
c	1986 ported to IBM PC Krakow, Poland
c
c	now, its 1989 and ported to IBM RT AIX v2.2.1
c	The RT F77 compiler doesn't like underscores in variable names

 * July 2003 -- port from FORTRAN to C
 */

/*
C
C   LATEST REVISION     - 18 May 1986
C
C   PURPOSE	     - ZERO OF A FUNCTION WHICH CHANGES SIGN IN A
C			   GIVEN INTERVAL (BRENT ALGORITHM)
C
C   USAGE	       - FindZero (F,EPS,NSIG,A,B,MAXFN)
C
C   ARGUMENTS    F      - AN EXTERNAL FUNCTION SUBPROGRAM F(X)
C			   PROVIDED BY THE USER WHICH COMPUTES F FOR
C			   ANY X IN THE INTERVAL (A,B). (INPUT)
C			   F MUST APPEAR IN AN EXTERNAL STATEMENT IN
C			   THE CALLING PROGRAM
C		EPS    - FIRST CONVERGENCE CRITERION (INPUT).  A ROOT,
C			   B, IS ACCEPTED IF ABS(F(B)) IS LESS THAN OR
C			   EQUAL TO EPS.  EPS MAY BE SET TO ZERO.
C		NSIG   - SECOND CONVERGENCE CRITERION (INPUT).  A ROOT,
C			   B, IS ACCEPTED IF THE CURRENT APPROXIMATION
C			   AGREES WITH THE TRUE SOLUTION TO NSIG
C			   SIGNIFICANT DIGITS.
C		A,B    - ON INPUT, THE USER MUST SUPPLY TWO POINTS, A
C			   AND B, SUCH THAT F(A) AND F(B) ARE OPPOSITE
C			   IN SIGN.
C		MAXFN  - ON INPUT, MAXFN SHOULD CONTAIN AN UPPER BOUND
C			   ON THE NUMBER OF FUNCTION EVALUATIONS
C			   REQUIRED FOR CONVERGENCE.  ON OUTPUT, MAXFN
C			   WILL CONTAIN THE ACTUAL NUMBER OF FUNCTION
C			   EVALUATIONS USED.
C
C   REMARKS  1.  ON EXIT FROM ZBRENT, WHEN IER=0, A AND B SATISFY THE
C		FOLLOWING,
C		F(A)*F(B) .LE.0,
C		ABS(F(B)) .LE. ABS(F(A)), AND
C		EITHER ABS(F(B)) .LE. EPS OR
C		ABS(A-B) .LE. MAX(ABS(B),0.1)*10.0**(-NSIG).
C		THE PRESENCE OF 0.1 IN THIS ERROR CRITERION CAUSES
C		LEADING ZEROES TO THE RIGHT OF THE DECIMAL POINT TO BE
C		COUNTED AS SIGNIFICANT DIGITS. SCALING MAY BE REQUIRED
C		IN ORDER TO ACCURATELY DETERMINE A ZERO OF SMALL
C		MAGNITUDE.
C	    2.  ZBRENT IS GUARANTEED TO REACH CONVERGENCE WITHIN
C		K = (ALOG((B-A)/D)+1.0)**2 FUNCTION EVALUATIONS WHERE
C		  D=MIN(OVER X IN (A,B) OF
C		    MAX(ABS(X),0.1)*10.0**(-NSIG)).
C		THIS IS AN UPPER BOUND ON THE NUMBER OF EVALUATIONS.
C		RARELY DOES THE ACTUAL NUMBER OF EVALUATIONS USED BY
C		ZBRENT EXCEED SQRT(K). D CAN BE COMPUTED AS FOLLOWS,
C		  P = AMIN1(ABS(A),ABS(B))
C		  P = AMAX1(0.1,P)
C		  IF ((A-0.1)*(B-0.1).LT.0.0) P = 0.1
C		  D = P*10.0**(-NSIG)
C
*/

#define fmax(x,y) ((x)>(y))?(x):(y)

double 
FindZero ((double (*func)(double)), bound, int nsig, 
          double a, 
          double b, maxit)
{
	int case1, case2, case3;
	int NumOfIter;
	double sigfigs;
	double x,a,b,fb;
	double delta, bisector, tolernace;
	double p,q,r,s,tmp;

	sigfigs = pow (10.0, (-nsig));
	x = a;
	fa = func (x);
	x = b;
	fb = func (x);
	NumOfIter = 2;

	if (0.0 < fa*fb)
	{
      printf ("FINDZERO: no zero bracketed by the specified interval\n");
		return b;
	}

	/* Reorganize the endpoints of the interval containing the zero
	 * so that the zero is closer to endpoint b
	 */
	c = a;
	fc = fa;
	if (fabs(fc) < fabs(fb)) 
	{
		fa = fb;
		fb = fc;
		fc = fa;
		a = b;
		b = c;
		c = a;
	}
	delta = c-b;
	bisector = 0.5 * delta;

	/*   Now, we start looping and continue to do so until the
	 *   convergence criteria have been satisfied
	 */
	tolerance = sigfigs * fmax (fabs (b), 0.1);
	while (1)
	{
		if ((fabs (fb) < bound) ||
		    (fabs (c-b) < tolerance)) break;

/*
c	       if we have only two points,
c	       the best we can do is linear interpolation
c	       otherwise, try an inverse quadratic interpolation
*/
		s = fb / fa;
		if ((a == c) || (2 == NumOfIter)) 
		{
			p = (b-a) * s;
			q = 1.0 - s;
		}
		else
		{
			q = fa/fc;
			r = fb/fc;
			tmp = r - 1.0;
			p = s * ((c-b) * q * (q-r) - (b-a) * tmp);
			q = - (q - 1.0) * tmp * (s - 1.0);
		}

/*
c	       delta =p/q will be the distance to the next estimated zero.
c	       we will use this value, if 1), 2) and 3) hold
c	       1) if the new guess lies between b and c
c		       (this to make sure we don't hop out)
c	       2) if the new guess is closer to endpoint b than endpoint c
c		       (case2 is to prevent oscillatory convergence).
c	       3) if the new delta is less than half of the old delta
c		       (case3 to insure convergence).
c	       if any one not true, then we bisect

c	       We also have a fourth possibility:
c	       4) if the new delta is better than the required tolerance
c	       then, chances are, the true zero is within tolerance of b,
c	       so we will step along with steps of size tolerance.

c	       note that the testing is done so as to avoid division by zero
*/

		tmp = bisector * q;
		x = p * tmp;
		if ((0.0 < x) && (x <= tmp*tmp) && (fabs (p) < fabs (0.5 * delta * q)))
		{
			delta = p/q;
		}
		else 
		{
			delta = bisector;
		}
		x = 0.5 * tolerance;
		if (fabs(delta) <= x) delta = copysign (x, bisector);

		a = b;
		fa = fb;
		b = b + delta;
		x = b;
		fb = func (x);
		NumOfIter ++;

/*
c	       arrange the endpoints so that the zero is bracketed
c	       by the interval (b,c) and is closer to endpoint b
*/
		if ((fb*fc >= 0.0) && (0.0 != fc))
		{
			x = c;
			c = b;
			b = a;
			a = x;
			x = fc;
			fc = fb;
			fb = fa;
			fa = x;
		}
		if (fabs (fc) < fabs (fb))
		{
			x = b;
			b = c;
			c = x;
			x = fb;
			fb = fc;
			fc = x;
		}
		tolerance = sigfigs * fmax (fabs(b), 0.1);
		bisector = 0.5 * (c-b);

		if (NumOfIter > maxit)
		{
			printf ("exceeded allowed num of evals on zero finding\n");
			break;
		}
	}

	/* we have found the zero. bracket it and return */
	a = c;
	maxit = NumOfIter;
	return b;
}
