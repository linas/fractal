/**
 * stirling.h
 *
 * Stirling numbers of the first and second kind.
 */
#ifndef STIRLING_H__
#define STIRLING_H__
#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>

// Stirling function of the first kind (signed)
long double complex cstirling_first (long double complex z, int m);

// Stirling function of the second kind
long double complex cstirling_second (long double complex z, int m);

#ifdef __cplusplus
};
#endif

#endif /* STIRLING_H__ */

