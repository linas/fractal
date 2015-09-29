/*
 * tree.c 
 *
 * draw a binary tree on the hyperbolic disk, using postscript
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


void eps_setup_disk (void)
{
	eps_setup_basic_linstyles();
}
void eps_setup_plane (void)
{
	eps_setup_basic_linstyles();
	printf ("0.0 -0.8 translate\n");
	printf ("0.4 1.6 scale\n");
}

/* ==================================================== */

/* draw three-pointed stick figure */
void draw_tristar (mobius_t m)
{
	cplex za,zb;

	za = cplex_set (0.0, 0.0);

	// printf ("n 0.0 0.0 m 0.25 0.433012702 l s\n");
	zb = cplex_set (0.25, 0.25*sqrt(3.0));

	// za = mobius_xform (m,za);
	// zb = mobius_xform (m,zb);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	draw_arc (m, za, zb);

	// printf ("n 0.0 0.0 m 0.25 -0.433012702 l s\n");
	zb = cplex_set (0.25, -0.25*sqrt(3.0));
	// zb = mobius_xform (m,zb);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	draw_arc (m, za, zb);

	// printf ("n 0.0 0.0 m -0.5 0.0 l s\n");
	zb = cplex_set (-0.5, 0.0);
	// zb = mobius_xform (m,zb);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	draw_arc (m, za, zb);

	// draw a splat 
	printf ("0.0600000 slw\n");
	za = cplex_set (-0.23, 0.0);
	za = mobius_xform (m,za);
	zb = cplex_set (-0.27, 0.0);
	zb = mobius_xform (m,zb);
	printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);
	// printf ("0.0100000 slw\n");
	printf ("0.010000 slw\n");
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
	draw_arc (m, za, zb);

	zb = cplex_set (0.25, -0.25*sqrt(3.0));
	// zb = mobius_xform (m,zb);
	// printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);
	// eps_draw_lineseg (za.re, za.im,zb.re, zb.im);
	draw_arc (m, za, zb);

	mobius_t tip = go_to_fork_tip(+1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);

	tip = go_to_fork_tip(-1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);
}


void draw(int n)
{
	mobius_t m;
	int level=3;

	cplex z = cplex_set (-0.268, 0.0);
	// cplex z = cplex_set (0.0, 0.0);
	mobius_t off = disk_center (z);

	mobius_t rot = mobius_rotate (-0.5*M_PI);
	// mobius_t rot = mobius_rotate ((-0.5-0.166666)*M_PI);
	off = mobius_mul (rot, off);

#define XLATE
#ifdef XLATE
	int a,b,c,d;
	a=1;
	b=n;
	c=0;
	d=1;
	mobius_t xfm = mobius_set (a,b,c,d);
	xfm = to_disk (xfm);
	off = mobius_mul (xfm, off);
#endif

#define HALF_PLANE
#ifdef HALF_PLANE
	xfm = to_half_plane_xform();
	off = mobius_mul (xfm, off);
#endif

	draw_tristar(off);

eps_set_color_green();
	m = go_to_fork_tip(+1.0);
	m = mobius_mul(off,m);
	draw_fork (m, level);

eps_set_color_red();
	m = go_to_fork_tip(-1.0);
	m = mobius_mul(off,m);
	draw_fork (m, level);

eps_set_color_blue();
	m = go_to_fork_tip(-3.0);
	m = mobius_mul(off,m);
	draw_fork (m, level);
}

/* ==================================================== */

int
main (int argc, char * argv[]) 
{

	if (argc < 1)
	{
		fprintf (stderr, "Usage: %s <shift>\n", argv[0]);
		exit(1);
	}

	int n = atoi (argv[1]);

	eps_print_prolog(220,220);
	// eps_print_prolog(400,100);
	eps_setup_disk();
	// eps_setup_plane();
	eps_draw_circle();

	draw(n);

	printf ("showpage\n");
}
