
/*
 * matrix.h
 *
 * Matrix utils
 *
 * Linas Jan 2004
 */

#ifndef __MATRIX_H__
#define __MATRIX_H__

void
multmatrix (int dim, matrix *prod, matrix *ml, matrix *mr);

void
identmatrix (int dim, matrix *e, long double val);

void
copymatrix (int dim, matrix *to, matrix *from);

void
addmatrix (int dim, matrix *to, matrix *afrom, matrix *bfrom);

void
lammatrix (int dim, matrix *to, long double lambda);

void
scalematrix (int dim, matrix *to, long double scale);

#endif /* __MATRIX_H__ */
