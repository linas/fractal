
/*
 * Eigenvalue finding
 * Call lapack from C code.
 *
 * Linas Vepstas September 2004
 */

#include <complex.h>

extern void dgeev_ (char *jobvl,  // compute left eigenvecs
                    char *jobvr,  // compute right eigenvecs
                    int *n,       // order of the matrix
                    double *a,    // the matrix
                    int *lda,     // leading dimension ??
                    double *wr,   // dim n real part of eigenvalues
                    double *wi,   // dim n imag part of eigenvalues
                    double *vl,   // left eigenvectors 
                    int *ldvl,    // dim of left eignevecs
                    double *vr,   // right eignevecs
                    int * ldvr,   // dim of right eigenvec array
                    double *work, // workspace
                    int * lwork,  // size of workspace, >= 4*N
                    int * info);  // returned success code


extern void zgeev_ (char *jobvl,  // compute left eigenvecs
                    char *jobvr,  // compute right eigenvecs
                    int *n,       // order of the matrix
                    complex double *a,    // the matrix
                    int *lda,     // leading dimension ??
                    complex double *w,   // dim n eigenvalues
                    complex double *vl,   // left eigenvectors 
                    int *ldvl,    // dim of left eignevecs
                    complex double *vr,   // right eignevecs
                    int * ldvr,   // dim of right eigenvec array
                    complex double *work, // workspace
                    int * lwork,  // size of workspace, >= 4*N
                    double * rwork,  // dimension 2*N
                    int * info);  // returned success code

