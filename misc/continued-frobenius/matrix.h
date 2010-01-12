
/*
 * matrix.h
 *
 * Matrix utils
 *
 * Linas Jan 2004
 */

#ifndef __MATRIX_H__
#define __MATRIX_H__

typedef long double matrix[MS][MS];
typedef long double vector[MS];

void
multmatrix (matrix *prod, matrix *ml, matrix *mr);

void
identmatrix (matrix *e, long double val);

void
copymatrix (matrix *to, matrix *from);

void
addmatrix (matrix *to, matrix *afrom, matrix *bfrom);

void
lammatrix (matrix *to, long double lambda);

void
scalematrix (matrix *to, long double scale);

#endif /* __MATRIX_H__ */
