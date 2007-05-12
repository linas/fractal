/*
 * mod-tree.c 
 * draw a binary tree on the hyperbolic upper half-plane
 * (or on poincare disk)
 * Built algebraically from modular xforms
 * as opposed to the geometric construction in tree.c
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

void draw_geo (mobius_t m, cplex a, cplex b)
{
	cplex z0 = mobius_xform (m,a);
	cplex z1 = mobius_xform (m,b);

	/* if the points are above one-another... */
	double dx = z1.re - z0.re;
	if ((1e-4>dx) && (-1e-4 < dx))
	{
		eps_draw_lineseg (z0.re, z0.im, z1.re, z1.im);
		return;
	}
	
	/* else draw geodesic */
	cplex zm = cplex_add (z0,z1);
	zm = cplex_scale (0.5, zm);
	double slope = (z1.im-z0.im)/(z1.re-z0.re);
	double xcenter = slope*zm.im+zm.re;
	
	double t0 = atan2 (z0.im, z0.re - xcenter);
	double t1 = atan2 (z1.im, z1.re - xcenter);
	int npts =20;
	double delta = (t1-t0)/npts;

	double dc = cos(delta);
	double ds = sin(delta);

	double xa = z0.re-xcenter;
	double ya = z0.im;

	int i;
	for (i=0; i<npts; i++)
	{
		double xb = xa*dc - ya*ds;
		double yb = xa*ds + ya*dc;
		eps_draw_lineseg (xa+xcenter, ya, xb+xcenter,yb);
		xa = xb;
		ya = yb;

	}
}

/* ========================================================= */

/**
 * @depth: recursion depth
 * @lr: 0 or 1, indicating that the left or rightt tree is being drawn
 * @draw_fund: 0 or 1, do draw the fundamental domain lines, or not
 */
void recursive_draw_binary_tree (int depth, int lr, int draw_fund, mobius_t m)
{
	if (0 >= depth)
		return;
	depth --;

	cplex tip;
	if (lr) tip = cplex_set(0.5,0.5*sqrt(3.0));
	else tip = cplex_set(-0.5,0.5*sqrt(3.0));

eps_set_color_red();
	mobius_t are;
	if (lr) are = mobius_set (1,1,0,1);
	else are = mobius_set (1,-1,0,1);

	mobius_t mr = mobius_mul (m,are);

	cplex zb = mobius_xform (are, tip);
	draw_geo (m, tip, zb);

	recursive_draw_binary_tree (depth, lr, draw_fund, mr);

	// seems to have the right curves
eps_set_color_green();
	mobius_t ell;
	if (lr) ell = mobius_set (1,0,1,1);
	else ell = mobius_set (1,0,-1,1);

	zb = mobius_xform (ell, tip);
	draw_geo (m, tip, zb);

	mobius_t ml = mobius_mul (m,ell);
	recursive_draw_binary_tree (depth, lr, draw_fund, ml);

	if (draw_fund)
	{
		eps_set_color(240,130,0);
		printf ("[0.02 0.01 0.005 0.01] 1 setdash\n");

		// good
		cplex tap = cplex_set(0.5,0.5*sqrt(3.0));
		cplex top = cplex_set(0,0);
		draw_geo (m, top, tap);

#if 0
		top = cplex_set(0,1160);
		cplex tep = cplex_set(-0.5,0.5*sqrt(3.0));
		draw_arc (m, top, tep);
#endif

		eps_set_color(0,70,220);
		printf ("[] 0 setdash\n");
	}
}

/* ========================================================= */

void draw (int n)
{
	mobius_t ident = mobius_ident();
	mobius_t tee = mobius_set (1,1,0,1);
	mobius_t ess = mobius_set (0,1,-1,0);

	mobius_t vee = mobius_set (0,-1,1,0);
	mobius_t pee = mobius_set (0,-1,1,1);

	mobius_t ell = mobius_set (1,0,1,1);
	mobius_t are = mobius_set (1,1,0,1);

	mobius_t xfm = ident;

// #define DISK_COORDS
#ifdef DISK_COORDS
	/* The following sets up a transform to disk coords */
	xfm = to_disk_xform (xfm);
#endif

#if 0
	/* Assorted translations and rotations */
	xfm = mobius_mul (xfm, ell);
	xfm = mobius_mul (xfm, are);
	xfm = mobius_mul (xfm, are);
	xfm = mobius_mul (xfm, pee);
	xfm = mobius_mul (xfm, pee);
	xfm = mobius_mul (xfm, ess);
	xfm = mobius_mul (xfm, tee);
	xfm = mobius_mul (xfm, ess);
	xfm = mobius_mul (tee,xfm);
	xfm = mobius_mul (ess,xfm);
	xfm = mobius_mul (tee,xfm);
	xfm = mobius_mul (tee,xfm);
	xfm = mobius_mul (ess,xfm);
#endif 

eps_set_color_blue();
// eps_set_color_green();
	recursive_draw_binary_tree (n,1, 1, xfm);
#if 0
// eps_set_color_red();
	recursive_draw_binary_tree (n,0, 1, xfm);
#endif

	cplex ltip = cplex_set(-0.5,0.5*sqrt(3.0));
	cplex rtip = cplex_set(0.5,0.5*sqrt(3.0));

// eps_set_color_blue();
	draw_geo (xfm, ltip, rtip);

	ltip = cplex_set(0.5,0.5*sqrt(3.0));
	rtip = cplex_set(1.5,0.5*sqrt(3.0));

	draw_geo (xfm, ltip, rtip);


	// draw a splat 
	printf ("0.0600000 slw\n");
eps_set_color_red();
	printf ("n %f %f m %f %f l s\n", -0.04, 1.0, 0.04, 1.0);
eps_set_color_green();
	printf ("n %f %f m %f %f l s\n", 0.96, 0.0, 1.04, 0.0);
	printf ("0.010000 slw\n");
}

/* ==================================================== */

int
main (int argc, char * argv[]) 
{
	if (argc < 2)
	{
		fprintf (stderr, "Usage: %s <shift>\n", argv[0]);
		exit(1);
	}

	int n = atoi (argv[1]);

	// eps_print_prolog(220,220);
	eps_print_prolog(800,200);
	// eps_setup_disk();
	eps_setup_plane();
	// eps_draw_circle();

	draw(n);

	printf ("showpage\n");
}
