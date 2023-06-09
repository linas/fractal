
/*
 * orbit.C
 *
 * Orbits of vectors under the dyadic representation 
 * of the modular group.
 *
 * Linas Vepstas November 2004
 */
#include <stdio.h>
#include <malloc.h>

#include "Farey.h"
#include "question.h"

double gee[2][2];
double are[2][2];

int ic=0;
int nsz=0;
double *xv;
double *yv;

void 
init (void)
{
	gee[0][0] = 1.0;
	gee[0][1] = 0.0;
	gee[1][0] = 0.0;
	gee[1][1] = 0.5;

	are[0][0] = 1.0;
	are[0][1] = 0.0;
	are[1][0] = 1.0;
	are[1][1] = -1.0;
}

void
recur (double mat[2][2], int lvl, int width)
{
	double tmp[2][2];
	int i;
	lvl --;
	if (0 == lvl) return;

	for (i=0; i<width; i++)
	{
		tmp[0][0] = gee[0][0] * mat[0][0] + gee[0][1] * mat[1][0];
		tmp[0][1] = gee[0][0] * mat[0][1] + gee[0][1] * mat[1][1];
		tmp[1][0] = gee[1][0] * mat[0][0] + gee[1][1] * mat[1][0];
		tmp[1][1] = gee[1][0] * mat[0][1] + gee[1][1] * mat[1][1];

		mat[0][0] = tmp[0][0];
		mat[0][1] = tmp[0][1];
		mat[1][0] = tmp[1][0];
		mat[1][1] = tmp[1][1];

		// printf ("%8.6g	%8.6g	w=%d	d=%d\n", mat[1][0], mat[1][1], width, lvl);
		if (ic >= nsz) return;
		xv[ic] = mat[1][0];
		yv[ic] = mat[1][1];
		ic ++;

		tmp[0][0] = are[0][0] * mat[0][0] + are[0][1] * mat[1][0];
		tmp[0][1] = are[0][0] * mat[0][1] + are[0][1] * mat[1][1];
		tmp[1][0] = are[1][0] * mat[0][0] + are[1][1] * mat[1][0];
		tmp[1][1] = are[1][0] * mat[0][1] + are[1][1] * mat[1][1];

		// printf ("%8.6g	%8.6g\n", tmp[1][0], tmp[1][1]);
		if (ic >= nsz) return;
		xv[ic] = tmp[1][0];
		yv[ic] = tmp[1][1];
		ic ++;

		recur (tmp, lvl, (width-1));
	}
}

void 
sortme (void)
{
	int i,j;
	for (i=0; i<ic; i++) 
	{
		for (j=i; j<ic; j++)
		{
			if (xv[i] > xv[j])
			{
				double tmp;
				tmp = xv[i];
				xv[i] = xv[j];
				xv[j] = tmp;

				tmp = yv[i];
				yv[i] = yv[j];
				yv[j] = tmp;
			}
		}
	}
}

int main (int argc, char * argv[])
{

	double mat[2][2];
	
	mat[0][0] = 1.0;
	mat[0][1] = 0.0;
	mat[1][0] = 0.0;
	mat[1][1] = 1.0;

	init ();
#define NMAX 123123
	nsz = NMAX;
	ic = 0;
	xv = (double *) malloc (nsz*sizeof (double));
	yv = (double *) malloc (nsz*sizeof (double));

	recur (mat, 3, 3);

	// sortme ();

	int i;
	for (i=0; i<ic; i++)
	{
#if 0
		double xf = question_inverse (xv[i]);
		double yf = question_inverse (yv[i]);
		printf ("%8.6g	%8.6g	%8.6g	%8.6g\n", xv[i], 0.0, xf, 0.0);
		printf ("%8.6g	%8.6g	%8.6g	%8.6g\n", xv[i], yv[i], xf, yf);
		printf ("%8.6g	%8.6g	%8.6g	%8.6g\n", xv[i], 0.0, xf, 0.0);
#endif
		printf ("%8.6g	%8.6g\n", xv[i], yv[i]);
	}
}
