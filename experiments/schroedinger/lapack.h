
/*
 * Eigenvalue finding
 * Call lapack from C code.
 *
 * Linas Vepstas September 2004
 */

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

extern void dstegr_(char *jobz,   // flag, compute eiegenvectors or not
                    char *range,  // flag, range of eigenvecs to be found
                    int *n,       // order of the matrix
                    double *d,    // diagonal elts, of length n
                    double *e,    // subdiagonal elts, of length n
                    double *vl,   // range low
                    double *vu,   // range high
                    int    *il,   // index range low
                    int    *iu,   // index range high
                    double *abstol,// absolute tolerance
                    int    *m,    // return number of eigenvalues found
                    double *w,    // length m array containing eigenvals
                    double *z,    // returned eigenvecs
                    int * ldz,    // leading dimension of array z
                    int *isuppz,  // eigenvec support
                    double *work, // work
                    int *lwork,   // dimension of work array
                    int *iwork,   // work
                    int *liwork,  // domension of work array
                    int *info);   // return code


