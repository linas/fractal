/*
 * similar.c
 *
 * Draw a binary tree on the poincare disk
 * Dervied from mod-tree.c
 *
 * Linas Vepstas April 2007
 * Updates Sept 2015
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

static int drawdash=0;

/* draw a formal geodesic */
/* Similar to the draw_arc routine used elsewhwere ... but different ... */
void draw_geo (mobius_t m, cplex a, cplex b)
{
	cplex z0 = mobius_xform (m,a);
	cplex z1 = mobius_xform (m,b);

#if 0
	/* vertical clip */
	if (100.0 < z0.im) {
		z0.im = 100;
		z0.re = z1.re;
	}
	if (100.0 < z1.im) {
		z1.im = 100;
		z1.re = z0.re;
	}
#endif

	/* if the points are above one-another... */
	double dx = z1.re - z0.re;
	if ((1e-4>dx) && (-1e-4 < dx))
	{
		eps_draw_lineseg (z0.re, z0.im, z1.re, z1.im);
		return;
	}
	
	/* else draw geodesic */
// #define HALF_PLANE
#ifdef HALF_PLANE
	cplex zm = cplex_add (z0,z1);
	zm = cplex_scale (0.5, zm);
	double slope = (z1.im - z0.im) / (z1.re - z0.re);
	double xcenter = slope * zm.im + zm.re;
	double ycenter = 0.0;
#endif

#define POINCARE_DISK
#ifdef POINCARE_DISK
	// Solve for small circle that passes through z0 and z1 and
	// is orthogonal to the unit circle.
	double cross = 2.0 * (z0.re*z1.im - z1.re*z0.im);
	double a2 = z0.re*z0.re + z0.im*z0.im + 1.0;
	double b2 = z1.re*z1.re + z1.im*z1.im + 1.0;
	double xcenter = ( z1.im * a2 - z0.im * b2 ) / cross;
	double ycenter = -( z1.re * a2 - z0.re * b2 ) / cross;
#endif

	double t0 = atan2 (z0.im - ycenter, z0.re - xcenter);
	double t1 = atan2 (z1.im - ycenter, z1.re - xcenter);

	if (M_PI < fabs(t1-t0))
	{
		if (t0 < 0.0) t0 += 2.0*M_PI;
		if (t1 < 0.0) t1 += 2.0*M_PI;
	}

#ifdef POSTSCRIPT_ARC
	double radius = sqrt((z0.re - xcenter)*(z0.re - xcenter)+z0.im*z0.im);
	
	if (t0>t1) {
		cplex tmp = z0;
		z0 = z1;
		z1 = tmp;

		double t = t0;
		t0 = t1;
		t1 = t;
	}
	printf ("n %f %f m %f %f %f %f %f a s\n",
		z0.re, z0.im, xcenter, 0.0,  radius, t0*180/M_PI, t1*180/M_PI);
#endif

#define HOMEMADE_ARC
#ifdef HOMEMADE_ARC
	int npts = 20*fabs(t1-t0);
	if (0 == npts) npts = 1;
	double delta = (t1-t0)/ ((double) npts);

	double dc = cos(delta);
	double ds = sin(delta);

	double xa = z0.re-xcenter;
	double ya = z0.im-ycenter;

	int i;
	double len = 0.0;
	for (i=0; i<npts; i++)
	{
		double xb = xa*dc - ya*ds;
		double yb = xa*ds + ya*dc;
		if (drawdash) {
			printf ("[0.02 0.01 0.005 0.01] %f setdash\n", len);
			len += sqrt ((yb-ya)*(yb-ya) + (xb-xa)*(xb-xa));
		}
		eps_draw_lineseg (xa+xcenter, ya+ycenter, xb+xcenter,yb+ycenter);
		xa = xb;
		ya = yb;
	}
#endif
}

/* ========================================================= */

/**
 * @depth: recursion depth
 * @lr: 0 or 1, indicating that the left or right tree is being drawn
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

	mobius_t are;
	if (lr) are = mobius_set (1,1,0,1);
	else are = mobius_set (1,-1,0,1);

	mobius_t mr = mobius_mul (m,are);

	cplex zb = mobius_xform (are, tip);
	draw_geo (m, tip, zb);

	recursive_draw_binary_tree (depth, lr, draw_fund, mr);

	// seems to have the right curves
	mobius_t ell;
	if (lr) ell = mobius_set (1,0,1,1);
	else ell = mobius_set (1,0,-1,1);

	zb = mobius_xform (ell, tip);
	draw_geo (m, tip, zb);

	mobius_t ml = mobius_mul (m,ell);
	recursive_draw_binary_tree (depth, lr, draw_fund, ml);

	if (draw_fund)
	{
		cplex top;
		eps_set_color(240,130,0);
		printf ("[0.02 0.01 0.005 0.01] 0 setdash\n");
		drawdash=1;

		// the ones on the right side ...
		cplex tap;
		if (lr) tap = cplex_set(0.5,0.5*sqrt(3.0));
		else tap = cplex_set(-0.5,0.5*sqrt(3.0));
		top = cplex_set(0,0);
		draw_geo (m, top, tap);

/*
		// the cusps on the left side ...
		top = cplex_set(0,1e8);
		cplex tep;
		if (lr) tep = cplex_set(-0.5,0.5*sqrt(3.0));
		else tep = cplex_set(0.5,0.5*sqrt(3.0));
		draw_geo (m, top, tep);
*/

		eps_set_color(0,70,220);
		printf ("[] 0 setdash\n");
		drawdash=0;
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

	/* The following sets up a transform to disk coords */
	xfm = to_disk_xform ();

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
eps_set_color(0,70,220);
// eps_set_color_green();
	recursive_draw_binary_tree (n, 1, 1, xfm);
// eps_set_color_red();
	recursive_draw_binary_tree (n, 0, 1, xfm);

	// Draw the center-line
	cplex ltip = cplex_set(-0.5,0.5*sqrt(3.0));
	cplex rtip = cplex_set(0.5,0.5*sqrt(3.0));
	draw_geo (xfm, ltip, rtip);
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

	eps_print_prolog(800,800);
	eps_setup_disk();
	eps_draw_circle();

	draw(n);

	printf ("showpage\n");
}
