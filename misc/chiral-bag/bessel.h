/*******************************************************************
 * FILE: bessel.h
 * FUNCTION: Compute spherical bessel functions
 * HISTORY: Created by Linas Vepstas <linas@linas.org> in 1983
 * c	1983 VAX/VMS 11/780 at SUNYSBNP
 * C	1987- Prime os 19.1  Sacaly, France.  The Goddamn Compiler
 * C		doesn't like for-- end for structures.  Revert to 
 * c		using line numbers.
 * C	1987- IBM PC -- MSDOS (Saclay, France; Krakow, Poland)
 * c	1988 vax/vmx 11/780 (visitng StonyBrook)
 * C	1989- IBM RT/ AIX 2.2.1   -- F77 doesn't like underscores
******************************************************************/

#ifndef _BESSEL_H_
#define _BESSEL_H_

/*
C       THIS SUBROUTINE CALCULATES THE SPHERICAL BESSEL FUNCTION FOR
C       ALL ORDERS .LE. N FOR (POSITIVE) ARGUMENT Z
C       THE RESULTS ARE RETURNED IN THE ARRAY BESS (0:N)
*/
void bessel (int n, double x, double *bess);

#endif
