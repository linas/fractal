
/*
 * question.c
 *
 * Minkowski question mark function
 * (wrapper around the C++ code)
 * 
 * Linas Vepstas September 2006
 */

#include "Farey.h"

static ContinuedFraction work;

double question_mark (int num, int denom)
{
	work.SetRatio (num,denom);
	return work.ToFarey ();
}

double fquestion_mark (double x)
{
	work.SetReal (x);
	return work.ToFarey ();
}
