
/* 
 * orbit of vector under 3D representation of modular group
 *
 * Linas Vepstas November  2004
 * Copyright (c) 2004 Linas Vepstas <linas@linas.org>
 */

/* required include files */
#include <GL/gl.h>
#include <GL/glut.h>

#include "vvector.h"

double www = 0.7;

double *pts[3];
int nsz = 1000000;
int npts = 0;

double gee[3][3];
double ref[3][3];

double rescale = 10.0;

void 
geen (double pos[3][3], int depth)
{
	depth --;
	if (0 > depth) return;
	if (npts > nsz)
	{
		printf ("out of space\n");
		return;
	}
	
	double tmp[3][3];
	double mat[3][3];
	double rmat[3][3];
	COPY_MATRIX_3X3 (mat, pos);

	int i;
	for (i=0; i<5; i++)
	{
		MATRIX_PRODUCT_3X3 (tmp, gee, mat);
		COPY_MATRIX_3X3 (mat, tmp);
		pts[0][npts] = rescale * (mat[2][0] - 1.5);
		pts[1][npts] = rescale * mat[2][1];
		pts[2][npts] = rescale * mat[2][2];
		
#if 0
		pts[0][npts] = rescale * (mat[1][0] -0.5 );
		pts[1][npts] = rescale * mat[1][1];
		pts[2][npts] = rescale * mat[1][2];
#endif
		
		// printf ("its %g   %g   %g\n", pts[0][npts], pts[1][npts], pts[2][npts]);
		npts ++;

		
		MATRIX_PRODUCT_3X3 (rmat, ref, mat);
		geen (rmat, depth);
	}
}
	
void 
InitStuff (void) 
{
	double mat[3][3];

	IDENTIFY_MATRIX_3X3 (gee);
	IDENTIFY_MATRIX_3X3 (ref);

	gee[1][1] = 0.5;
	gee[2][2] = www;
	gee[2][1] = 1.0;

	ref[1][0] = 1.0;
	ref[1][1] = -1.0;

	pts[0] = (double *) malloc (nsz * sizeof (double)); 
	pts[1] = (double *) malloc (nsz * sizeof (double)); 
	pts[2] = (double *) malloc (nsz * sizeof (double)); 
	
	IDENTIFY_MATRIX_3X3 (mat);
	geen (mat, 6);

	printf ("made %d points\n", npts);
}

void 
DrawStuff (void) 
{
	glBegin (GL_POINTS);
	// glBegin (GL_LINE_STRIP);

	// glScalef (rescale, rescale, rescale);

	int i;
	for (i=0; i<npts; i++)
	{
		glVertex3d (pts[0][i], pts[1][i], pts[2][i]);
		// glNormal3d (pts[0][i], pts[1][i], pts[2][i]);
	}

	glEnd();
}

/* ------------------------- end of file ----------------- */
