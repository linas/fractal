/*
 * FILE: bag_ener.h
 * FUNCTION: compute energy levels of firee quarks in chiral bag
 * HISTORY: Created by Linas Vepstas <linas@linas.org> 1983
 *   July 2003 ported from FORTRAN to C
 */

int
quark_energy (double *ener, double theta, int ispect, int k, 
               double emax,  int kpty);
