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
C
C   LATEST REVISION     - 18 May 1986
C
C   PURPOSE             - ZERO OF A FUNCTION WHICH CHANGES SIGN IN A
C                           GIVEN INTERVAL (BRENT ALGORITHM)
C
C   USAGE               - CALL FindZero (F,EPS,NSIG,A,B,MAXFN)
C
C   ARGUMENTS    F      - AN EXTERNAL FUNCTION SUBPROGRAM F(X)
C                           PROVIDED BY THE USER WHICH COMPUTES F FOR
C                           ANY X IN THE INTERVAL (A,B). (INPUT)
C                           F MUST APPEAR IN AN EXTERNAL STATEMENT IN
C                           THE CALLING PROGRAM
C                EPS    - FIRST CONVERGENCE CRITERION (INPUT).  A ROOT,
C                           B, IS ACCEPTED IF ABS(F(B)) IS LESS THAN OR
C                           EQUAL TO EPS.  EPS MAY BE SET TO ZERO.
C                NSIG   - SECOND CONVERGENCE CRITERION (INPUT).  A ROOT,
C                           B, IS ACCEPTED IF THE CURRENT APPROXIMATION
C                           AGREES WITH THE TRUE SOLUTION TO NSIG
C                           SIGNIFICANT DIGITS.
C                A,B    - ON INPUT, THE USER MUST SUPPLY TWO POINTS, A
C                           AND B, SUCH THAT F(A) AND F(B) ARE OPPOSITE
C                           IN SIGN.
C                           ON OUTPUT, BOTH A AND B ARE ALTERED.  B
C                           WILL CONTAIN THE BEST APPROXIMATION TO THE
C                           ROOT OF F. SEE REMARK 1.
C                MAXFN  - ON INPUT, MAXFN SHOULD CONTAIN AN UPPER BOUND
C                           ON THE NUMBER OF FUNCTION EVALUATIONS
C                           REQUIRED FOR CONVERGENCE.  ON OUTPUT, MAXFN
C                           WILL CONTAIN THE ACTUAL NUMBER OF FUNCTION
C                           EVALUATIONS USED.
C
C   REMARKS  1.  ON EXIT FROM ZBRENT, WHEN IER=0, A AND B SATISFY THE
C                FOLLOWING,
C                F(A)*F(B) .LE.0,
C                ABS(F(B)) .LE. ABS(F(A)), AND
C                EITHER ABS(F(B)) .LE. EPS OR
C                ABS(A-B) .LE. MAX(ABS(B),0.1)*10.0**(-NSIG).
C                THE PRESENCE OF 0.1 IN THIS ERROR CRITERION CAUSES
C                LEADING ZEROES TO THE RIGHT OF THE DECIMAL POINT TO BE
C                COUNTED AS SIGNIFICANT DIGITS. SCALING MAY BE REQUIRED
C                IN ORDER TO ACCURATELY DETERMINE A ZERO OF SMALL
C                MAGNITUDE.
C            2.  ZBRENT IS GUARANTEED TO REACH CONVERGENCE WITHIN
C                K = (ALOG((B-A)/D)+1.0)**2 FUNCTION EVALUATIONS WHERE
C                  D=MIN(OVER X IN (A,B) OF
C                    MAX(ABS(X),0.1)*10.0**(-NSIG)).
C                THIS IS AN UPPER BOUND ON THE NUMBER OF EVALUATIONS.
C                RARELY DOES THE ACTUAL NUMBER OF EVALUATIONS USED BY
C                ZBRENT EXCEED SQRT(K). D CAN BE COMPUTED AS FOLLOWS,
C                  P = AMIN1(ABS(A),ABS(B))
C                  P = AMAX1(0.1,P)
C                  IF ((A-0.1)*(B-0.1).LT.0.0) P = 0.1
C                  D = P*10.0**(-NSIG)
C
        SUBROUTINE FindZero (func, bound, NSIG, A,B, maxit)
        implicit double precision (a-h, o-z)
        logical *1 case1, case2, case3, happy
        double precision func
        external func
        sigfigs = 10.0d0 ** (-NSIG)
        x = a
        fa = func (x)
        x = b
        fb = func (x)
        NumOfIter = 2
c
        IF (FA*FB .GT. 0.0d0)
     +  stop 'FINDZERO: no zero bracketed by the specified interval'
c
c       reorganize the endpoints of the interval containing the zero
c       so that the zero is closer to endpoint b
        C = A
        FC = FA
        IF (DABS(FC) .LT. DABS(FB)) THEN
                FA = FB
                FB = FC
                FC = FA
                A = B
                B = C
                C = A
        END IF
        delta = c-b
        bisector = 0.5d0 * delta
c
c       now, we start looping and continue to do so until the
c       convergence criteria have been satified
        tolerance = sigfigs * dmax1 (dabs (b), 0.1d0)
        DO 1001 IIIIII=1,1000000
                happy = dabs (fb) .lt. bound
                happy = happy .or. (dabs (c-b) .le. tolerance)
                if (happy) goto 1002
c
c               if we have only two points,
c               the best we can do is linear interpolation
c               otherwise, try an inverse quadratic interpolation
                s = fb / fa
                IF ((A .EQ. C) .OR. (NumOfIter .EQ. 2)) THEN
                        p = (b-a) * s
                        Q = 1.0d0 - S
                ELSE
                        Q = FA/FC
                        R = FB/FC
                        TMP = R - 1.0d0
                        P = S * ((C-B) * Q * (Q-R) - (B-A) * TMP)
                        Q = - (Q - 1.0d0) * TMP * (S - 1.0d0)
                END IF
c
c               delta =p/q will be the distance to the next estimated zero.
c               we will use this value, if 1), 2) and 3) hold
c               1) if the new guess lies between b and c
c                       (this to make sure we don't hop out)
c               2) if the new guess is closer to endpoint b than endpoint c
c                       (case2 is to prevent oscillatory convergence).
c               3) if the new delta is less than half of the old delta
c                       (case3 to insure convergence).
c               if any one not true, then we bisect

c               We also have a fourth possibility:
c               4) if the new delta is better than the required tolerance
c               then, chances are, the true zero is within tolerance of b,
c               so we will step along with steps of size tolerance.

c               note that the testing is done so as to avoid division by zero

                tmp = bisector * q
                x = p * tmp
                case1 = 0.0d0 .lt. x
                case2 = x .le. tmp*tmp
                case3 = dabs (p) .lt. dabs (0.5d0 * delta * q)
                delta = bisector
                if (case1 .and. case2 .and. case3) delta = p/q
                x= 0.5d0 * tolerance
                if (dabs(delta) .le. x) delta = dsign (x, bisector)

                a = b
                fa = fb
                b = b + delta
                x = b
                fb = func (x)
                NumOfIter = NumOfIter + 1

c               arrange the endpoints so that the zero is bracketed
c               by the interval (b,c) and is closer to endpoint b
                if ((dsign(1.0d0,fc) * dsign(1.0d0,fb) .gt. 0.0d0)
     +                  .and. (fc .ne. 0.0d0)) then
                        x = c
                        c = b
                        b = a
                        a = x
                        x = fc
                        fc = fb
                        fb = fa
                        fa = x
                end if
                if (dabs (fc) .lt. dabs (fb)) then
                        x = b
                        b = c
                        c = x
                        x = fb
                        fb = fc
                        fc = x
                end if
                tolerance = sigfigs * DMAX1 (DABS(B), 0.1D0)
                bisector = 0.5d0 * (c-b)
c
                if (NumOfIter .ge. maxit)
     +          stop 'exceeded allowed num of evals on zero finding'
1001    CONTINUE
1002    CONTINUE
c
c       we have found the zero. bracket it and return
        A = C
        maxit = NumOfIter
        return
        end
