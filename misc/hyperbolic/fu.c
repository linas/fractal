/*
 * fu.c 
 *
 * Draw a binary tree on the hyperbolic disk, using postscript
 * This is an explictly geometric construction. For an algebraic
 * construction of the same thing (which is a tad simpler), see
 * mod-tree.c
 *
 * Linas Vepstas April 2007
 */

#include "cplex.h"
#include "flt.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



/* arc */
void print_arc(mobius_t m, cplex zf, cplex zt)
{
	int i;
	int nseg=1;

	// clip anything taller than yclip
	double yclip = 5;
	if (yclip <zf.im)
	{
		double m = (zt.im-zf.im)/(zt.re-zf.re);
		double xclip = zt.re +(yclip-zt.im)/m;
		zf.re = xclip;
		zf.im = yclip;
	}
	if (yclip <zt.im)
	{
		double m = (zt.im-zf.im)/(zt.re-zf.re);
		double xclip = zf.re +(yclip-zf.im)/m;
		zt.re = xclip;
		zt.im = yclip;
	}

	cplex zstart = zf;
	cplex zdelta = cplex_scale((1.0/(double)nseg), cplex_sub (zt,zf));
	cplex za = mobius_xform (m, zstart);
	cplex zb;

	cplex zs = za;

	for (i=0; i<nseg; i++)
	{
		cplex zend = cplex_add (zstart, zdelta);
		zb = mobius_xform (m, zend);
		double dist = cplex_dist (zs, zb);
		if (dist > 0.03) {
			printf("duuude %Lg %Lg\n", zs.re, zs.im);
			printf("other end= %Lg %Lg\n", zb.re, zb.im);
			zs = zb;
		}
		zstart = zend;
		za = zb;
	}
}

/* ==================================================== */

mobius_t go_to_fork_tip(double sign)
{
	cplex c = cplex_set (-0.5, 0.0);
	mobius_t m = disk_center (c);
	cplex rot = cplex_exp_itheta (sign*2.0*M_PI/6.0);
	m = mobius_scale (m,rot);
// show_mobius (m);
	return m;
}

void draw_fork(mobius_t m, int level)
{
	if (level == 0) return;
	level--;

	cplex za, zb;

	za = cplex_set (0.0, 0.0);

	zb = cplex_set (0.25, 0.25*sqrt(3.0));
	// za = mobius_xform (m,za);
	// zb = mobius_xform (m,zb);
	// printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	print_arc (m, za, zb);

/****
	zb = cplex_set (0.25, -0.25*sqrt(3.0));
	// zb = mobius_xform (m,zb);
	// printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	print_arc (m, za, zb);
****/

	mobius_t tip = go_to_fork_tip(+1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);

/***
	tip = go_to_fork_tip(-1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);
***/
}


void draw()
{
	mobius_t m;
	int level=3;
level=20;

	// Originally drawn with -0.268 -- this gets the vertex centers
	// correctly located.
	double cent = sqrt(3.0) - 2.0;
	cplex z = cplex_set (cent, 0.0);
	// z = cplex_set (-0.25, 0.0); xxxxxx
	z = cplex_set (0.0, 0.0);
	mobius_t off = disk_center (z);

	mobius_t rot = mobius_rotate (-0.5*M_PI);
	// mobius_t rot = mobius_rotate ((-0.5-0.166666)*M_PI);
	off = mobius_mul (rot, off);

#define HALF_PLANE
#ifdef HALF_PLANE
	mobius_t xfm = mobius_set (1,0,0,1);
	xfm = to_half_plane_xform();
	off = mobius_mul (xfm, off);
#endif

	m = go_to_fork_tip(-1.0);
	m = mobius_mul(off,m);
	draw_fork (m, level, 255,0,0);
}

/* ==================================================== */

int
main (int argc, char * argv[]) 
{
	draw();
}
