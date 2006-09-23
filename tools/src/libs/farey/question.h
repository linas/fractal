
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
double question_inverse (double);

#ifdef   __cplusplus
};
#endif
