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

// Stirling function of the first kind (signed)
long int stirling_first (int n, int k);

// Stirling function of the second kind
long int stirling_second (int n, int k);

#ifdef __cplusplus
};
#endif

#endif /* STIRLING_H__ */

