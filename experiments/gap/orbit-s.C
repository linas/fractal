
/*
 * orbit-s.C
 *
 * Orbits of vectors under the dyadic representation 
 * of the modular group.
 *
 * Linas Vepstas November 2004
 */
#include <stdio.h>
#include <malloc.h>
#include <math.h>

double gee[3][3];
double ige[3][3];
double are[3][3];

int ic=0;
int nsz=0;
double *xv;
double *yv;

void 
plain_dyadic_rep_init (void)
{
	gee[0][0] = 0.5;
	gee[0][1] = 0.0;
	gee[0][2] = 0.0;
	
	gee[1][0] = 0.0;
	gee[1][1] = 0.5;
	gee[1][2] = 0.0;
	
	gee[2][0] = 0.0;
	gee[2][1] = 0.0;
	gee[2][2] = 1.0;

	are[0][0] = -1.0;
	are[0][1] = 0.0;
	are[0][2] = 1.0;
	
	are[1][0] = 0.0;
	are[1][1] = 1.0;
	are[1][2] = 0.0;
	
	are[2][0] = 0.0;
	are[2][1] = 0.0;
	are[2][2] = 1.0;
}

// generates a butterfly moth thing at the origin
void 
init_butterfly (void)
{
	gee[0][0] = 0.5;
	gee[0][1] = 0.0;
	gee[0][2] = 0.0;
	
	gee[1][0] = 0.0;
	gee[1][1] = 0.5*sqrt(3);
	gee[1][2] = 0.0;
	
	gee[2][0] = 0.0;
	gee[2][1] = 0.0;
	gee[2][2] = 1.0;

	// reflection about line 30 degrees up 
	// for vector v, its reflection is
	// v_ref = v - (v.n) n 
	// where n is normal to the reflection plane
	// for triangle, n = (-1/2, -sqrt(3)/2)
	// matrix is outer product -2n x n
	//
	are[0][0] = 0.5;
	are[0][1] = 0.5*sqrt(3);
	are[0][2] = 0.0;
	
	are[1][0] = 0.5*sqrt(3);
	are[1][1] = -0.5;
	are[1][2] = 0.0;
	
	are[2][0] = 0.0;
	are[2][1] = 0.0;
	are[2][2] = 1.0;
}

// geneates dyadic/farey accumulating level sets
// wow !! density of accumulation may be farey-distribution-like !!
// bin this !!! Fantastic!
// Wow, the dist is almost the question mark, but not quite ... 
// weirdly distorted!
void 
init_farey_levelsets (void)
{
	gee[0][0] = 0.5;
	gee[0][1] = 0.0;
	gee[0][2] = 0.0;
	
	gee[1][0] = 0.0;
	gee[1][1] = 0.5*sqrt(3);
	gee[1][2] = 0.0;
	
	gee[2][0] = 0.0;
	gee[2][1] = 0.0;
	gee[2][2] = 1.0;

	// reflecion about x=0.5 line
	are[0][0] = -1.0;
	are[0][1] = 0.0;
	are[0][2] = 1.0;
	
	are[1][0] = 0.0;
	are[1][1] = 1.0;
	are[1][2] = 0.0;
	
	are[2][0] = 0.0;
	are[2][1] = 0.0;
	are[2][2] = 1.0;
}

void 
init (void)
{
	gee[0][0] = 0.5;
	gee[0][1] = 0.0;
	gee[0][2] = 0.0;
	
	gee[1][0] = 0.0;
	gee[1][1] = 0.5*sqrt(3);
	gee[1][2] = 0.0;
	
	gee[2][0] = 0.0;
	gee[2][1] = 0.0;
	gee[2][2] = 1.0;

	ige[0][0] = 2.0;
	ige[0][1] = 0.0;
	ige[0][2] = 0.0;
	
	ige[1][0] = 0.0;
	ige[1][1] = 2.0/sqrt(3);
	ige[1][2] = 0.0;
	
	ige[2][0] = 0.0;
	ige[2][1] = 0.0;
	ige[2][2] = 1.0;

	are[0][0] = -1.0;
	are[0][1] = 0.0;
	are[0][2] = 1.0;
	
	are[1][0] = 0.0;
	are[1][1] = 1.0;
	are[1][2] = 0.0;
	
	are[2][0] = 0.0;
	are[2][1] = 0.0;
	are[2][2] = 1.0;

#if 0
	// for vector v, its reflection is
	// v_ref = v - (v.n) n 
	// where n is normal to the reflection plane
	// for triangle, n = (-1/2, -sqrt(3)/2)
	// matrix is outer product -2n x n
	//
	are[0][0] = 0.5;
	are[0][1] = 0.5*sqrt(3);
	are[0][2] = 0.0;
	
	are[1][0] = 0.5*sqrt(3);
	are[1][1] = -0.5;
	are[1][2] = 0.0;
	
	are[2][0] = 0.0;
	are[2][1] = 0.0;
	are[2][2] = 1.0;
#endif
}

void
recur (double mat[3][3], int lvl, int width)
{
	double tmp[3][3];
	double sav[3][3];
	int i;
	if (0 == lvl) return;
	lvl --;

	sav[0][0] = mat[0][0];
	sav[0][1] = mat[0][1];
	sav[0][2] = mat[0][2];
	sav[1][0] = mat[1][0];
	sav[1][1] = mat[1][1];
	sav[1][2] = mat[1][2];
	sav[2][0] = mat[2][0];
	sav[2][1] = mat[2][1];
	sav[2][2] = mat[2][2];

	for (i=0; i<width; i++)
	{
		tmp[0][0] = gee[0][0] * mat[0][0] + gee[0][1] * mat[1][0] + gee[0][2] * mat[2][0];
		tmp[0][1] = gee[0][0] * mat[0][1] + gee[0][1] * mat[1][1] + gee[0][2] * mat[2][1];
		tmp[0][2] = gee[0][0] * mat[0][2] + gee[0][1] * mat[1][2] + gee[0][2] * mat[2][2];
		
		tmp[1][0] = gee[1][0] * mat[0][0] + gee[1][1] * mat[1][0] + gee[1][2] * mat[2][0];
		tmp[1][1] = gee[1][0] * mat[0][1] + gee[1][1] * mat[1][1] + gee[1][2] * mat[2][1];
		tmp[1][2] = gee[1][0] * mat[0][2] + gee[1][1] * mat[1][2] + gee[1][2] * mat[2][2];
		
		tmp[2][0] = gee[2][0] * mat[0][0] + gee[2][1] * mat[1][0] + gee[2][2] * mat[2][0];
		tmp[2][1] = gee[2][0] * mat[0][1] + gee[2][1] * mat[1][1] + gee[2][2] * mat[2][1];
		tmp[2][2] = gee[2][0] * mat[0][2] + gee[2][1] * mat[1][2] + gee[2][2] * mat[2][2];
		
		mat[0][0] = tmp[0][0];
		mat[0][1] = tmp[0][1];
		mat[0][2] = tmp[0][2];
		mat[1][0] = tmp[1][0];
		mat[1][1] = tmp[1][1];
		mat[1][2] = tmp[1][2];
		mat[2][0] = tmp[2][0];
		mat[2][1] = tmp[2][1];
		mat[2][2] = tmp[2][2];

		// printf ("%8.6g	%8.6g	w=%d	d=%d\n", mat[1][0], mat[1][1], width, lvl);
		if (ic >= nsz) return;
		xv[ic] = mat[0][2];
		yv[ic] = mat[1][2];
		ic ++;

		tmp[0][0] = are[0][0] * mat[0][0] + are[0][1] * mat[1][0] + are[0][2] * mat[2][0];
		tmp[0][1] = are[0][0] * mat[0][1] + are[0][1] * mat[1][1] + are[0][2] * mat[2][1];
		tmp[0][2] = are[0][0] * mat[0][2] + are[0][1] * mat[1][2] + are[0][2] * mat[2][2];
		
		tmp[1][0] = are[1][0] * mat[0][0] + are[1][1] * mat[1][0] + are[1][2] * mat[2][0];
		tmp[1][1] = are[1][0] * mat[0][1] + are[1][1] * mat[1][1] + are[1][2] * mat[2][1];
		tmp[1][2] = are[1][0] * mat[0][2] + are[1][1] * mat[1][2] + are[1][2] * mat[2][2];
		
		tmp[2][0] = are[2][0] * mat[0][0] + are[2][1] * mat[1][0] + are[2][2] * mat[2][0];
		tmp[2][1] = are[2][0] * mat[0][1] + are[2][1] * mat[1][1] + are[2][2] * mat[2][1];
		tmp[2][2] = are[2][0] * mat[0][2] + are[2][1] * mat[1][2] + are[2][2] * mat[2][2];
		
		// printf ("%8.6g	%8.6g\n", tmp[1][0], tmp[1][1]);
		if (ic >= nsz) return;
		xv[ic] = tmp[0][2];
		yv[ic] = tmp[1][2];
		ic ++;

		recur (tmp, lvl, (width-1));
	}

#if INVERSE_TOO
	mat[0][0] = sav[0][0];
	mat[0][1] = sav[0][1];
	mat[0][2] = sav[0][2];
	mat[1][0] = sav[1][0];
	mat[1][1] = sav[1][1];
	mat[1][2] = sav[1][2];
	mat[2][0] = sav[2][0];
	mat[2][1] = sav[2][1];
	mat[2][2] = sav[2][2];

	for (i=0; i<width; i++)
	{
		tmp[0][0] = ige[0][0] * mat[0][0] + ige[0][1] * mat[1][0] + ige[0][2] * mat[2][0];
		tmp[0][1] = ige[0][0] * mat[0][1] + ige[0][1] * mat[1][1] + ige[0][2] * mat[2][1];
		tmp[0][2] = ige[0][0] * mat[0][2] + ige[0][1] * mat[1][2] + ige[0][2] * mat[2][2];
		
		tmp[1][0] = ige[1][0] * mat[0][0] + ige[1][1] * mat[1][0] + ige[1][2] * mat[2][0];
		tmp[1][1] = ige[1][0] * mat[0][1] + ige[1][1] * mat[1][1] + ige[1][2] * mat[2][1];
		tmp[1][2] = ige[1][0] * mat[0][2] + ige[1][1] * mat[1][2] + ige[1][2] * mat[2][2];
		
		tmp[2][0] = ige[2][0] * mat[0][0] + ige[2][1] * mat[1][0] + ige[2][2] * mat[2][0];
		tmp[2][1] = ige[2][0] * mat[0][1] + ige[2][1] * mat[1][1] + ige[2][2] * mat[2][1];
		tmp[2][2] = ige[2][0] * mat[0][2] + ige[2][1] * mat[1][2] + ige[2][2] * mat[2][2];
		
		mat[0][0] = tmp[0][0];
		mat[0][1] = tmp[0][1];
		mat[0][2] = tmp[0][2];
		mat[1][0] = tmp[1][0];
		mat[1][1] = tmp[1][1];
		mat[1][2] = tmp[1][2];
		mat[2][0] = tmp[2][0];
		mat[2][1] = tmp[2][1];
		mat[2][2] = tmp[2][2];

		// printf ("%8.6g	%8.6g	w=%d	d=%d\n", mat[1][0], mat[1][1], width, lvl);
		if (ic >= nsz) return;
		xv[ic] = mat[0][2];
		yv[ic] = mat[1][2];
		ic ++;

		tmp[0][0] = are[0][0] * mat[0][0] + are[0][1] * mat[1][0] + are[0][2] * mat[2][0];
		tmp[0][1] = are[0][0] * mat[0][1] + are[0][1] * mat[1][1] + are[0][2] * mat[2][1];
		tmp[0][2] = are[0][0] * mat[0][2] + are[0][1] * mat[1][2] + are[0][2] * mat[2][2];
		
		tmp[1][0] = are[1][0] * mat[0][0] + are[1][1] * mat[1][0] + are[1][2] * mat[2][0];
		tmp[1][1] = are[1][0] * mat[0][1] + are[1][1] * mat[1][1] + are[1][2] * mat[2][1];
		tmp[1][2] = are[1][0] * mat[0][2] + are[1][1] * mat[1][2] + are[1][2] * mat[2][2];
		
		tmp[2][0] = are[2][0] * mat[0][0] + are[2][1] * mat[1][0] + are[2][2] * mat[2][0];
		tmp[2][1] = are[2][0] * mat[0][1] + are[2][1] * mat[1][1] + are[2][2] * mat[2][1];
		tmp[2][2] = are[2][0] * mat[0][2] + are[2][1] * mat[1][2] + are[2][2] * mat[2][2];
		
		// printf ("%8.6g	%8.6g\n", tmp[1][0], tmp[1][1]);
		if (ic >= nsz) return;
		xv[ic] = tmp[0][2];
		yv[ic] = tmp[1][2];
		ic ++;

		recur (tmp, lvl, (width-1));
	}
#endif 

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
	int i;
	double mat[3][3];
	
	mat[0][0] = 1.0;
	mat[0][1] = 0.0;
	mat[0][2] = 0.0;
	mat[1][0] = 0.0;
	mat[1][1] = 1.0;
	mat[1][2] = 1.0;
	mat[2][0] = 0.0;
	mat[2][1] = 0.0;
	mat[2][2] = 1.0;

	init ();
#define NMAX 123123
	nsz = NMAX;
	ic = 0;
	xv = (double *) malloc (nsz*sizeof (double));
	yv = (double *) malloc (nsz*sizeof (double));

	recur (mat, 6, 8);

	// sortme ();

#if 0
	for (i=0; i<ic; i++)
	{
		printf ("%8.6g	%8.6g\n", xv[i], yv[i]);
	}
#endif

	// make bin counts 
	int nbins = 400;
	int bin[5000];
	int n;
	for (n=0; n<nbins; n++)
	{
		bin[n] = 0;
	}

	int cnt = 0;
	for (i=0; i<ic; i++)
	{
		// if (yv[i] < 0.09) continue;
		n = (int) (((double)nbins) *xv[i]);
		bin[n] += 1;
		cnt ++;
	}

	double acc = 0.0;
	for (n=0; n<nbins; n++)
	{
		double dis = bin[n];
		dis /= ((double) cnt);
		acc += dis;
		double x = n;
		x /= nbins;
		printf ("%d	%8.6g	%8.6g	%8.6g\n", n, x, dis, acc);
	}

}
