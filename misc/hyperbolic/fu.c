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
#include "eps.h"
#include "flt.h"
#include "flt-eps.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/* ==================================================== */
int drawcusp = 1;

/* draw fundamental domains */
void draw_cusp (mobius_t m)
{
	cplex za,zb;

	if (!drawcusp) return;
	// Draw cusps in orange, dashed linestyle.
	eps_set_color(240,130,0);
	printf ("[0.02 0.01 0.005 0.01] 0 setdash\n");

	za = cplex_set (0.0, 0.0);
	zb = cplex_set (1.0, 0.0);
	draw_arc (m, za, zb);

	zb = cplex_set (-0.5, 0.5*sqrt(3.0));
	draw_arc (m, za, zb);

	zb = cplex_set (-0.5, -0.5*sqrt(3.0));
	draw_arc (m, za, zb);

	printf ("[] 0 setdash\n");
}

mobius_t go_to_fork_tip(double sign)
{
	cplex c = cplex_set (-0.5, 0.0);
	mobius_t m = disk_center (c);
	cplex rot = cplex_exp_itheta (sign*2.0*M_PI/6.0);
	m = mobius_scale (m,rot);
// show_mobius (m);
	return m;
}

void draw_fork(mobius_t m, int level, int r, int g, int b)
{
	if (level == 0) return;
	level--;

	cplex za, zb;

	if (drawcusp)
	{
		draw_cusp(m);
		eps_set_color(r,g,b);
	}

	za = cplex_set (0.0, 0.0);

	zb = cplex_set (0.25, 0.25*sqrt(3.0));
	// za = mobius_xform (m,za);
	// zb = mobius_xform (m,zb);
	// printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	draw_arc (m, za, zb);

	zb = cplex_set (0.25, -0.25*sqrt(3.0));
	// zb = mobius_xform (m,zb);
	// printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	draw_arc (m, za, zb);

	mobius_t tip = go_to_fork_tip(+1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level, r, g, b);

	tip = go_to_fork_tip(-1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level, r, g, b);
}


void draw()
{
	mobius_t m;
	int level=3;
level=10;

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

eps_set_color_red();
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
