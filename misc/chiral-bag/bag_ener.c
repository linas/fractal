        SUBROUTINE ENERGY (theta, ispect, kay, emax, kpty)

C       THIS SUBROUTINE RETURNS zeros of an ENERGY EIGENVALUE equation,
c       with parameters theta, ispect, kay, emax

C       ******* INPUT *******

C       THETA- CHIRAL ANGLE
C       IPTY-   PARAMETER INDICATING PARITY
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


        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EXTERNAL EFN, YEFN
        COMMON /ENER/ numFound, ENER (2)

C       MAXIT-  MAXIMUM NUMBER OF ITERATIONS THE CONVERGENCE ALGORITHM
C       IS TO PERFORM WHEN SEARCHING FOR THE EIGENVALUES.  IN PRACTICE,
C       IT SEEMS TO TAKE LESS THAN TEN ITERATIONS-- SUGGEST MAXIT =15
C       if there is failure, an error message will be printed
        MAXIT = 40
C       NSIG-   NUMBER OF SIGNIFICANT DIGITS FOR THE SEARCH ROUTINE,
C       I.E. THAT THE ENERGIES ARE TO BE CALCULATED TO.
c       NSIG = 15   /* this will work on the Vax and the PC but ...*/
        NSIG = 14

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


C       initialize the equation whos parameters are to be found
        call efnInit (theta, ispect, kay, emax, kpty)

        A = 0.0D0
        sygnus = DBLE (ISPECT)
        B = dmax1 (DBLE (KAY-2), 3.2d-4)
        STEP = 0.1D0
        FUNB = DSIGN (1.0D0, EFN(sygnus * b))
        numfound = 0
        justfound = -99
c       CHECK TO SEE IF WE WANT TO FIND MORE ENERGIES
        DO 1001 IIIIII= 1,10000000
                IF (A .GT. EMAX) GOTO 1002
C               DEFINE THE NEXT INTERVAL
c               if we have not recently found a zero, step along with
c               step size step. but if we have found a zero recently,
c               take a big step and save some cpu time.
                if (justfound .lt. 0) then
                        a = b
                        b = b + step
                        funa = funb
                        funb = dsign (1.0d0, efn (sygnus * b))
                else
                        a = dabs (ener (numFound)) + step
                        if (numfound .le. 20) then
                                b = a + 2.2d0
                        else
                                b = a + 2.95d0
                        end if
                        funa = dsign (1.0d0, efn (sygnus * a))
                        funb = dsign (1.0d0, efn (sygnus * b))
                        justfound = -99
                end if

                IF (FUNA*FUNB .LE. 0.0D0) THEN
C               FOUND A ZERO IN INTERVAL.  CALL CONVERGENCE ALGORITHM.
C               ** NOTE ** AA, BB, AND M ARE MODIFIED BY find_zero.
C               THESE VARIABLES MUST BE RESET BEFORE EACH CALL.
                        AA = sygnus * A
                        BB = sygnus * B
                        M = MAXIT
                        CALL findzero (EFN, 0.0D0,NSIG,AA,BB,M)
c                       STORE AND PRINT THE ENERGIES FOUND
                        numfound = numfound  + 1
                        ENER (numfound) = BB
c	print *,' its ', bb
C                       ***** DECREASE STEP SIZE PER COMMENTS ABOVE
                        justfound = 99
                        if (numfound .gt. 1) then
                                nnn = numfound
                                dil = ener(nnn) - ener (nnn-1)
                                dil = dabs (dil) * 0.4d0
                                step = dmin1 (step, dil)
                        end if
                END IF
1001    CONTINUE
1002    CONTINUE

        RETURN
        END

c       *********************************************************

        subroutine efninit (theta, ispect, kay, emax, kpty)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMMON /EFN1/ PTYK, CO, SI, K

        PTYK = DBLE (KPTY)
        K = KAY
        CO = DCOS (THETA)
        SI = DSIN (THETA)

        return
        end

c       *********************************************************

        DOUBLE PRECISION FUNCTION EFN (OMEGA)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMMON /EFN1/ PTYK, CO, SI, KAY
        IF (KAY.GT.0) THEN
                CALL QUICKbessel (kay, omega, bkm1, bk, bkp1)
                BKS = BK*BK
                EFN = CO * (BKP1 * BKM1 - BKS)
     +               + PTYK * BK * (BKP1 - BKM1)
     +               + SI * BKS / OMEGA
        ELSE
                SN = DSIN (OMEGA)
                CS = DCOS (OMEGA)
                B0 = SN/OMEGA
                B1 = (B0 - CS) / OMEGA
                EFN = CO*B1-(PTYK-SI)*B0
        END IF
        RETURN
        END

c       *********************************************************

        DOUBLE PRECISION FUNCTION YEFN (OMEGA)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMMON /EFN1/ PTYK, CO, SI, KAY
        IF (KAY.GT.0) THEN
                CALL QUICKbessneu
     +          (kay, omega, bkm1, bk, bkp1, ykm1, yk, ykp1)
                yks = yk*yk
                YEFN = CO*(YKP1*YKM1-YKS)+PTYK*YK*(YKP1-YKM1)+SI*YKS/OMEGA
        ELSE
                SN = DSIN (OMEGA)
                CS = DCOS (OMEGA)
                B0 = -CS/OMEGA
                B1 = (B0 - SN)/OMEGA
                YEFN = CO*B1-(PTYK-SI)*B0
        END IF
        RETURN
        END
