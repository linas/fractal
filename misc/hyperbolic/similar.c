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
#include <string.h>

void eps_setup_disk (void)
{
	eps_setup_basic_linstyles();
}

static int drawdash=0;

/* draw a formal geodesic */
/* Similar to the draw_arc routine used elsewhwere ... but different ... */
void draw_geo (mobius_t m, cplex a, cplex b)
{
	// So .. why this, and why not draw_arc() ?? The coordinate
	// systems seemto be slightly off, one from the other, I don't
	// know why. Seems strange to me, this should all have worked...
	// draw_arc(m,a,b); return;

	cplex z0 = mobius_xform (m,a);
	cplex z1 = mobius_xform (m,b);

// #define HALF_PLANE
#ifdef HALF_PLANE
	/* vertical clip */
	if (100.0 < z0.im) {
		z0.im = 100;
		z0.re = z1.re;
	}
	if (100.0 < z1.im) {
		z1.im = 100;
		z1.re = z0.re;
	}

	/* If the points are above one-another... */
	double dx = z1.re - z0.re;
	if ((1e-4>dx) && (-1e-4 < dx))
	{
		eps_draw_lineseg (z0.re, z0.im, z1.re, z1.im);
		return;
	}
	
	/* else draw geodesic */
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
	if (fabs(cross) < 1.0e-10) cross = 1.0e-10; // avoid divide-by-zero.
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

		// the ones on the left side ...
		cplex tap;
		if (lr) tap = cplex_set(0.5,0.5*sqrt(3.0));
		else tap = cplex_set(-0.5,0.5*sqrt(3.0));
		top = cplex_set(0,0);
		draw_geo (m, top, tap);

		// the cusps on the right side ...
		top = cplex_set(0,1e8);
		cplex tep;
		if (lr) tep = cplex_set(-0.5,0.5*sqrt(3.0));
		else tep = cplex_set(0.5,0.5*sqrt(3.0));
		draw_geo (m, top, tep);

		eps_set_color(0,70,220);
		printf ("[] 0 setdash\n");
		drawdash=0;
	}
}

/* ========================================================= */

/* draw a small rectangle */
/* Observe that the center of the disk is located at i */
/* This is consistent with upper half-plane coordinates! */
void draw_splat (mobius_t m)
{
#ifdef DRAW_CROSS
	cplex za,zb;

	double xsplat = 0.0;
	double ysplat = 1.3;

	// draw a cross
	printf ("0.00800000 slw\n");
	za = cplex_set (xsplat-0.04, ysplat);
	za = mobius_xform (m,za);
	zb = cplex_set (xsplat+0.10, ysplat);
	zb = mobius_xform (m,zb);
	printf ("n %Lf %Lf m %Lf %Lf l s\n", za.re, za.im, zb.re, zb.im);
	za = cplex_set (xsplat, ysplat);
	za = mobius_xform (m,za);
	zb = cplex_set (xsplat, ysplat+0.18);
	zb = mobius_xform (m,zb);
	printf ("n %Lf %Lf m %Lf %Lf l s\n", za.re, za.im, zb.re, zb.im);
	printf ("0.010000 slw\n");
#endif

	int i;
// eps_set_color(0xbb, 0xdd, 0x0);
	printf ("0.00800000 slw\n");
	for (i=0; i<15; i++)
	{
		cplex za,zb,zc,zd;
		double y = 0.96 + (i+1)*(i+1)*0.03;

		za = cplex_set (-0.5, y);
		za = mobius_xform (m,za);
		zb = cplex_set (-0.1666, y);
		zb = mobius_xform (m,zb);
		zc = cplex_set (0.1666, y);
		zc = mobius_xform (m,zc);
		zd = cplex_set (0.5, y);
		zd = mobius_xform (m,zd);
		printf ("n %Lf %Lf m %Lf %Lf l s\n", za.re, za.im, zb.re, zb.im);
		printf ("n %Lf %Lf m %Lf %Lf l s\n", zb.re, zb.im, zc.re, zc.im);
		printf ("n %Lf %Lf m %Lf %Lf l s\n", zc.re, zc.im, zd.re, zd.im);
	}
	printf ("0.010000 slw\n");

}

/* ========================================================= */

/* Draw splats at each similarity point */
void sim_splat(mobius_t local, mobius_t global)
{
	mobius_t vee = mobius_set (0,-1,1,0);   // v^2 == 1  V=-S
	mobius_t pee = mobius_set (0,-1,1,1);   // P^3 == 1
	mobius_t vin = mobius_set (0, 1, -1, 0); // v^-1 == V^2

// eps_set_color(0xaa, 0xcc, 0x0);
	mobius_t sim = mobius_mul(global, local);
	draw_splat(sim);

	// Now do the Vee similarities
	sim = mobius_mul(local, vee);
	sim = mobius_mul(vin, sim);
	sim = mobius_mul(global, local);
	draw_splat(sim);

	sim = mobius_mul(local, vin);
	sim = mobius_mul(vee, sim);
	sim = mobius_mul(global, local);
	draw_splat(sim);
}

/* ========================================================= */

/**
 * @depth: recursion depth
 * @lr: 0 or 1, indicating that the left or right tree is being drawn
 * @draw_fund: 0 or 1, do draw the fundamental domain lines, or not
 */
void recursive_draw_similar (int n, mobius_t local, mobius_t global)
{
	if (0 >= n)
		return;
	n--;

	mobius_t sim = local;
	sim_splat(sim, global);

	mobius_t ell = mobius_set (1, 0, 1, 1);
	mobius_t ill = mobius_set (1, 0, -1, 1);
	mobius_t are = mobius_set (1, 1, 0, 1);
	mobius_t ari = mobius_set (1, -1, 0, 1);
	mobius_t vee = mobius_set (0, -1, 1, 0); // v^3 == I
	mobius_t vin = mobius_set (0, 1, -1, 0); // v^-1 == V^2

	// Move in the L direction.
	sim = mobius_mul (ell, local);
	sim = mobius_mul (sim, ill);

// eps_set_color_red();
	sim_splat(sim, global);

	recursive_draw_similar(n, sim, global);

	// Move in the Linv direction
	sim = mobius_mul (ill, local);
	sim = mobius_mul (sim, ell);

// eps_set_color_green();
	sim_splat(sim, global);
	recursive_draw_similar(n, sim, global);

	// Move in the R direction.
	sim = mobius_mul (are, local);
	sim = mobius_mul (sim, ari);

// eps_set_color_blue();
	sim_splat(sim, global);
	recursive_draw_similar(n, sim, global);

	// Move in the Rinv direction.
	sim = mobius_mul (ari, local);
	sim = mobius_mul (sim, are);

// eps_set_color(0x0, 0xff, 0xff);
	sim_splat(sim, global);
	recursive_draw_similar(n, sim, global);
}

/* ========================================================= */

void similar (int n, const char *xform, mobius_t global)
{
	mobius_t ident = mobius_ident();
	mobius_t tee = mobius_set (1,1,0,1);    // T = R
	mobius_t ess = mobius_set (0,1,-1,0);

	mobius_t vee = mobius_set (0,-1,1,0);   // v^2 == 1  V=-S
	mobius_t pee = mobius_set (0,-1,1,1);   // P^3 == 1

	mobius_t ell = mobius_set (1,0,1,1);
	mobius_t are = mobius_set (1,1,0,1);
	mobius_t ill = mobius_set (1, 0, -1, 1);
	mobius_t ari = mobius_set (1, -1, 0, 1);

	mobius_t xfm = ident;

	int i;
	for (i=0; i<strlen(xform); i++)
	{
		char x = xform[i];
		if ('L' == x) xfm = mobius_mul(xfm, ell);
		if ('l' == x) xfm = mobius_mul(xfm, ill);

		if ('R' == x) xfm = mobius_mul(xfm, are);
		if ('r' == x) xfm = mobius_mul(xfm, ari);

		if ('S' == x) xfm = mobius_mul(xfm, ess);
		if ('T' == x) xfm = mobius_mul(xfm, tee);
		if ('t' == x) xfm = mobius_mul(xfm, ari);

		if ('P' == x) xfm = mobius_mul(xfm, pee);
		if ('V' == x) xfm = mobius_mul(xfm, vee);
	}

	// recursive_draw_similar(n, xfm, global);

eps_set_color(0x0, 0x0, 0x0);
	recursive_draw_similar(n, ident, global);

eps_set_color(0x99, 0x0, 0xdd);
	recursive_draw_similar(n, vee, global);

eps_set_color(0xdd, 0x0, 0xbb);
	recursive_draw_similar(n, pee, global);
eps_set_color(0xbb, 0x44, 0x77);
	xfm = ident;
	xfm = mobius_mul(xfm, pee);
	xfm = mobius_mul(xfm, pee);
	recursive_draw_similar(n, xfm, global);

eps_set_color(0xbb, 0xdd, 0x0);
	recursive_draw_similar(n, ell, global);
eps_set_color(0xdd, 0xaa, 0x0);
	recursive_draw_similar(n, are, global);

eps_set_color(0x0, 0x77, 0xdd);
	xfm = ident;
	xfm = mobius_mul(xfm, ell);
	xfm = mobius_mul(xfm, ell);
	recursive_draw_similar(n, xfm, global);
eps_set_color(0x0, 0xdd, 0x77);
	xfm = ident;
	xfm = mobius_mul(xfm, are);
	xfm = mobius_mul(xfm, are);
	recursive_draw_similar(n, xfm, global);

eps_set_color(0x0, 0xdd, 0xbb);
	xfm = ident;
	xfm = mobius_mul(xfm, ell);
	xfm = mobius_mul(xfm, ell);
	xfm = mobius_mul(xfm, ell);
	recursive_draw_similar(n, xfm, global);
eps_set_color(0x0, 0xaa, 0xdd);
	xfm = ident;
	xfm = mobius_mul(xfm, are);
	xfm = mobius_mul(xfm, are);
	xfm = mobius_mul(xfm, are);
	recursive_draw_similar(n, xfm, global);
}

/* ========================================================= */

mobius_t draw (int n)
{
	mobius_t ident = mobius_ident();
	mobius_t tee = mobius_set (1,1,0,1);    // T = R
	mobius_t ess = mobius_set (0,1,-1,0);

	mobius_t vee = mobius_set (0,-1,1,0);   // v^2 == 1  V=-S
	mobius_t pee = mobius_set (0,-1,1,1);   // P^3 == 1

	mobius_t ell = mobius_set (1,0,1,1);
	mobius_t are = mobius_set (1,1,0,1);

	mobius_t xfm = ident;

	/* The following sets up a transform to disk coords */
	/* Notice that the coordinates used throughout are upper-half-plane
	 * coordinates, and that this xform merely converts them to a disk
	 * when they are drawn. However, the "raw" z values used throughout
	 * are upper-half-plane z's.  This is a key point!
	 */
	xfm = to_disk_xform ();

#if 0
	// Relocate the center of the thing.
	// This is not as intuitive as you might think:
	// Although changing the real component does a hyperbolic rotation,
	// Changing the imaginary component scales the thing!!
	double cent = sqrt(3.0) - 2.0;
	cplex z = cplex_set (cent, 0.0);
	// z = cplex_set (-0.25, 0.0);
	// cplex z = cplex_set (0.0, 0.0);
	mobius_t off = disk_center (z);

	xfm = mobius_mul(xfm, off);

	// z = cplex_set (0.0, 1.0);
	// off = mobius_scale(off, z);
	// xfm = mobius_mul(off, xfm);
#endif

#if 0
	/* Assorted translations and rotations */
	xfm = mobius_mul (xfm, ell);
	xfm = mobius_mul (xfm, are);
	xfm = mobius_mul (xfm, vee);
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
	// Note that the center-line spans the bottom of the fun domain
	// and that the coordinates given are in upper-half-plane coordinates!
	cplex ltip = cplex_set(-0.5,0.5*sqrt(3.0));
	cplex rtip = cplex_set(0.5,0.5*sqrt(3.0));
	draw_geo (xfm, ltip, rtip);

	return xfm;
}

/* ==================================================== */

int
main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s <recurs> <xform>\n", argv[0]);
		exit(1);
	}

	int n = atoi (argv[1]);

	eps_print_prolog(800,800);
	eps_setup_disk();
	eps_draw_circle();

	mobius_t global = draw(n);
	similar(n, argv[2], global);

	printf ("showpage\n");
}
