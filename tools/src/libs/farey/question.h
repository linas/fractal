
/*
 * FUNCTION:
 * question.h
 *
 * Minkowski question mark function
 * (wrapper around the C++ code)
 *
 * HISTORY:
 * Linas Vepstas September 2006
 */

#ifdef   __cplusplus
extern "C" {
#endif

double question_mark (int num, int denom);
double fquestion_mark (double);

/**
 * dyadic_to_stern_brocot -- convert real number to stern-brocot
 *
 * This function maps the unit interval, expressed as a dyadic
 * tree, to the full stern brocot tree.  The input is 0<x<1.
 * Perform binary expansion of x. Use L,R for binary digits 0,1.
 * Concatenate the L and R to form a matrix. Pad on the right with
 * zeros. The resulting fractional linear transform will map the
 * entire complex plane to b/d. This function returns b/d.
 *
 * Basically, this function maps the dyadic tree to the
 * Stern-Brocot tree.
 */
long double dyadic_to_stern_brocot (long double x);

/**
 * question_inverse - return the inverse of the question mark function
 *
 * This implements a rapid algorithm to compute the inverse of
 * the question mark function.
 */
long double question_inverse (long double x);

/**
 * Computes the Stern-Brocot tree.
 * It takes as input a dyadic fraction N / 2^level
 * assuming 0 <= N <= 2^level.
 * It returns the corresponding p/q value in p and q.
 *
 * In practical terms, given an input x = N / 2^level
 * this returns the inverse question-mark, given by
 * (as per usual) ?(p/q) = N / 2^level.
 */
void stern_brocot_tree(unsigned long N, int level,
                       unsigned long *p, unsigned long *q);
#ifdef   __cplusplus
};
#endif
