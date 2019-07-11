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
   double   re_center,
   double   im_center,
   double   width,
   double   height,
   int      itermax,
   double   renorm)
{
	int cnt = 0;
	double freq;
	double fmi;
	while (true)
	{
		cnt ++;
		int rc = fscanf(stdin, "%g	%g\n", &freq, &fmi);
		if (EOF == rc)
		{
printf("duude fail cnt=%d rc=%d\n", cnt, rc);
		}
if (10 < cnt)break;
	}
}
