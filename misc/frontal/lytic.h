/*
 *
 * Numerical exploration of the analytic continuation
 * given in frontal.lyx
 *
 * September 2015
 */

#include <complex.h>
#include "../sum/sum.h"

int bitcount(int k, int bits[LEN]);
int xbitcount(int k, double x);
double complex count_extend(int k, int c_k, int c_1, double complex u);
double complex count_extend(int k, int bits[LEN], double complex u);
double complex count_extend(int k, double x, double complex u);
double complex sum_extend(double x, double complex u);

