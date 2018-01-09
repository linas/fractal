/*
 * julie.C
 *
 * Julia Set visualization
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/* Verbatim copu of julia.c */

#define NBITS 10
#define NPTS (1<<NBITS)

double map[NPTS+1];
void up(double v, double Kay, int loc, int lvl);
void down(double v, double Kay, int loc, int lvl)
{
	lvl --;
	v *= 2.0 * Kay;
   if (Kay < v) v = Kay;
	int cent = loc - (1<<lvl);
   map[cent] = v;

	if (lvl < 1) return;
	down(v, Kay, cent, lvl);
	up(v, Kay, cent, lvl);
}

void up(double v, double Kay, int loc, int lvl)
{
	lvl --;
	v *= 2.0 * Kay;
	v -= Kay;
   if (v< 0.0) v = 0.0;
	int cent = loc + (1<<lvl);
   map[cent] = v;

	if (lvl < 1) return;
	down(v, Kay, cent, lvl);
	up(v, Kay, cent, lvl);
}

void mkmap(double y, double Kay)
{
	for (int i=0; i<=NPTS; i++)
	{
		map[i] = 99999;
	}
	int lvl = NBITS-1;
	int cent = 1<<lvl;
	map[cent] = y;
	down (y, Kay, cent, lvl);
	up (y, Kay, cent, lvl);
}

static void julia_diagram (float *array,
                                 int array_size,
                                 double x_center,
                                 double x_width,
                                 double Kay,
                                 int itermax,
                                 double why)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	mkmap(Kay, Kay);
	for (int j=1; j<array_size; j++)
	{
		array[j] = map[j];
	}
}

DECL_MAKE_BIFUR(julia_diagram)
