/*
 * Scatterplots
 *
 * June 2019
 */

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "brat.h"

void MakeHisto (
   char     *name,     /* value of argv[0] */
   float    *glob,
   int      sizex,
   int      sizey,
   double   x_center,
   double   y_center,
   double   width,
   double   height,
   int      itermax,
   double   renorm)
{
	int cnt = 0;
	float freq;
	float fmi;

	while (true)
	{
		int rc = scanf("%f\t%f", &freq, &fmi);
		// printf("duuude got %d %d EOF=%d %g %g\n", cnt, rc, EOF, freq, fmi);
		if (0 == rc)
		{
			char bug[50];
			fgets(bug, 50, stdin); // evil
			continue;
		}
		if (EOF == rc)
		{
			break;
		}
		cnt ++;

freq = - log(freq) / log(2.0);

		int ni = sizex * (0.5 + (freq - x_center) / width);
		if (ni < 0) ni = 0;
		if (sizex <= ni) ni = sizex - 1;

		int nj = sizey * (0.5 + (fmi - y_center) / height);
		if (nj < 0) nj = 0;
		if (sizey <= nj) nj = sizey - 1;

		glob [nj*sizex + ni] += 1.0;
	}
}
