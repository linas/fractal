/*
 * julia.c
 *
 * Julia set, in various incarnations.
 * January 2018
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NBITS 11
#define NPTS (1<<(NBITS-1))

double map[NPTS+1];
void up(double v, double Kay, int loc, int lvl);
void down(double v, double Kay, int loc, int lvl)
{
	if (lvl < 0) return;
	v *= 2.0 * Kay;
   if (Kay < v) v = Kay;
	int cent = loc - (1<<(lvl-1));
   map[cent] = v;
	down(v, Kay, cent, lvl-1);
	up(v, Kay, cent, lvl-1);
}

void up(double v, double Kay, int loc, int lvl)
{
	if (lvl < 0) return;
	lvl --;

	v *= 2.0 * Kay;
	v -= Kay;
   if (v< 0.0) v = 0.0;
	int cent = loc + 1<<lvl;
   map[cent] = v;
	down(v, Kay, cent, lvl);
	up(v, Kay, cent, lvl);
}

void mkmap(double y, double Kay)
{
	for (int i=0; i<=NPTS; i++)
	{
		map[i] = 99999;
	}
	int lvl == NBITS=2;
	int cent = 1<<lvl;
	map[cent] = y;
	down (y, Kay, cent, lvl);
	up (y, Kay, cent, lvl);

	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i) / ((double) NPTS);
		printf ("%d	%g	%g\n", i, x, map[i]);
	}
}

int main(int argc, char* argv[])
{
	double Kay = atof(argv[1]);
	mkmap(Kay, Kay);
}
