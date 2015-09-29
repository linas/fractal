/*
 * similar.c
 *
 * draw a binary tree on the hyperbolic disk, using postscript
 * This is an explictly geometric construction. For an algebraic
 * construction of the same thing (which is a tad simpler), see
 * mod-tree.c
 *
 * Modified for visualizng similarity transforms.
 *
 * Linas Vepstas April 2007
 * Modified September 2015
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
	// printf ("0.0 -0.8 translate\n");
	// printf ("0.9 0.9 scale\n");
}

/* ==================================================== */

/* draw a small rectangle */
void draw_splat (mobius_t m)
{
	cplex za,zb;

	// draw a splat
	printf ("0.0600000 slw\n");
	// printf ("0.0100000 slw\n");
	// za = cplex_set (-0.23, 0.0);
	za = cplex_set (-0.02, 0.0);
	za = mobius_xform (m,za);
	// zb = cplex_set (-0.27, 0.0);
	zb = cplex_set (0.02, 0.0);
	zb = mobius_xform (m,zb);
	printf ("n %Lf %Lf m %Lf %Lf l s\n", za.re, za.im, zb.re, zb.im);
	// printf ("0.0100000 slw\n");
	printf ("0.010000 slw\n");
}

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

	cplex zz, zl, zr, zc;

	zz = cplex_set (0.0, 0.0);

	zl = cplex_set (0.25, 0.25*sqrt(3.0));
	draw_arc (m, zz, zl);

	zr = cplex_set (0.25, -0.25*sqrt(3.0));
	draw_arc (m, zz, zr);

	// Draw cusps in orange, dashed linestyle.
	eps_set_color(240,130,0);
	printf ("[0.02 0.01 0.005 0.01] 0 setdash\n");

	// This one isn't needed except for the tri-star...
	zc = cplex_set (1.0, 0.0);
	draw_arc (m, zz, zc);

	// Back up by one iteration...
	zz = cplex_set (-0.5, 0.0);
	zc = cplex_set (-1.0, 0.0);
	draw_arc (m, zz, zc);
	eps_set_color_black();
	printf ("[] 0 setdash\n");

	// Now recurse.
	mobius_t tip = go_to_fork_tip(+1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);

	tip = go_to_fork_tip(-1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);
}

/* Draw splats at each similarity point */
void sim_splat(mobius_t off)
{
	draw_splat(off);
}

/* Draw similarity point */
void similar(mobius_t off, int n)
{
	n--;
	if (0 == n) return;

	mobius_t ell = mobius_set (1, 0, 1, 1);
	mobius_t ill = mobius_set (1, 0, -1, 1);
	mobius_t are = mobius_set (1, 1, 0, 1);
	mobius_t ari = mobius_set (1, -1, 0, 1);
	mobius_t vee = mobius_set (0, -1, 1, 0);
	mobius_t vin = mobius_set (0, 1, -1, 0);
	ell = to_disk (ell);
	ill = to_disk (ill);
	are = to_disk (are);
	ari = to_disk (ari);

	// Move in the L direction.
	mobius_t sim;
	sim = mobius_mul (ell, off);
	sim = mobius_mul (sim, ill);

eps_set_color_red();
	sim_splat(sim);
	similar(sim, n);

	// Move in the Linv direction
	sim = mobius_mul (ill, off);
	sim = mobius_mul (sim, ell);

eps_set_color_green();
	sim_splat(sim);
	similar(sim, n);

	// Move in the R direction.
	sim = mobius_mul (are, off);
	sim = mobius_mul (sim, ari);

eps_set_color_blue();
	sim_splat(sim);
	similar(sim, n);

	// Move in the Rinv direction.
	sim = mobius_mul (ari, off);
	sim = mobius_mul (sim, are);

eps_set_color(0x0, 0xff, 0xff);
	sim_splat(sim);
	similar(sim, n);

	// Do the V versions
eps_set_color(0xcc, 0x0, 0xff);
	sim = mobius_mul (vee, ell);
	sim = mobius_mul (sim, off);
	sim = mobius_mul (sim, ill);
	sim = mobius_mul (sim, vin);
	sim_splat(sim);

	sim = mobius_mul (vee, ill);
	sim = mobius_mul (sim, off);
	sim = mobius_mul (sim, ell);
	sim = mobius_mul (sim, vin);
	sim_splat(sim);

	sim = mobius_mul (vee, are);
	sim = mobius_mul (sim, off);
	sim = mobius_mul (sim, ari);
	sim = mobius_mul (sim, vin);
	sim_splat(sim);

	sim = mobius_mul (vee, ari);
	sim = mobius_mul (sim, off);
	sim = mobius_mul (sim, are);
	sim = mobius_mul (sim, vin);
	sim_splat(sim);

}

/* Draw the entire tree */
void draw(int n)
{
	mobius_t m;
	int level = n;

	// Originally drawn with -0.268 -- this gets the vertex centers
	// correctly located.
	double cent = sqrt(3.0) - 2.0;
	cplex z = cplex_set (cent, 0.0);
	// z = cplex_set (-0.25, 0.0);
	// z = cplex_set (0.0, 0.0);
	mobius_t off = disk_center (z);

	mobius_t rot = mobius_rotate (-0.5*M_PI);
	// rot = mobius_rotate ((-0.5-0.166666)*M_PI);
	off = mobius_mul (rot, off);

// #define XLATE
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

	draw_splat(off);
	draw_tristar(off);

// eps_set_color_green();
	m = go_to_fork_tip(+1.0);
	m = mobius_mul(off,m);
	draw_fork (m, level);

// eps_set_color_red();
	m = go_to_fork_tip(-1.0);
	m = mobius_mul(off,m);
	draw_fork (m, level);

// eps_set_color_blue();
	m = go_to_fork_tip(-3.0);
	m = mobius_mul(off,m);
	draw_fork (m, level);


#if 0
	off = disk_center (z);
	rot = mobius_rotate (-0.5*M_PI);
	rot = mobius_rotate ((-0.5-0.166666)*M_PI);
	rot = mobius_rotate ((-0.5-0.3333333)*M_PI);
	// WTF? why is this one degenerate?
	// why is the next one doubled up?
	rot = mobius_rotate ((-0.5-0.5)*M_PI);
	rot = mobius_rotate ((-0.5-0.0833333)*M_PI);
// show_mobius(off);
	off = mobius_mul (rot, off);
// show_mobius(off);
#endif
	similar(off, n);
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
