/*
 * Compute matrix elements of fourier of minkowski question mark 
 * by means of brute-force integration. This seems to be a viable
 * quick-n-dirty way of getting these.
 *
 * Linas June 2008
 */

void set_npts(int);
void make_g_plain_elt(int m, int n, long double *pre, long double *pim);
void make_full_elt(int m, int n, long double *pre, long double *pim);
