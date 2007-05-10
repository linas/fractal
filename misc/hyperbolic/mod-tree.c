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

// eps_set_color_red();
	mobius_t are;
	if (lr) are = mobius_set (1,1,0,1);
	else are = mobius_set (1,-1,0,1);

	cplex zb = mobius_xform (are, tip);
	draw_seg (m, tip, zb);

	if (draw_fund)
	{
		eps_set_color(240,130,0);
		printf ("[0.02 0.01 0.005 0.01] 1 setdash\n");

		cplex top = cplex_set(-1,0);
		draw_seg (m, top, tip);

		eps_set_color(0,70,220);
		printf ("[] 0 setdash\n");
	}

	mobius_t mr = mobius_mul (m,are);
	recursive_draw_binary_tree (depth, lr, draw_fund, mr);

// eps_set_color_green();
	mobius_t ell;
	if (lr) ell = mobius_set (1,0,1,1);
	else ell = mobius_set (1,0,-1,1);

	zb = mobius_xform (ell, tip);
	draw_seg (m, tip, zb);

	mobius_t ml = mobius_mul (m,ell);
	recursive_draw_binary_tree (depth, lr, draw_fund, ml);
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

#if DISK_COORDS
	/* The following are inappropriate for disk coords */
	// xfm = to_disk (xfm);
	// xfm = disk_center (cplex_set (0,1));

	/* The following sets up a transform to disk coords */
	xfm.a = cplex_set (1,0);
	xfm.b = cplex_set (0,-1);
	xfm.c = cplex_set (1,0);
	xfm.d = cplex_set (0,1);
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
// eps_set_color_red();
	recursive_draw_binary_tree (n,0, 1, xfm);

	cplex ltip = cplex_set(0.5,0.5*sqrt(3.0));
	cplex rtip = cplex_set(-0.5,0.5*sqrt(3.0));

// eps_set_color_blue();
	draw_seg (xfm, ltip, rtip);
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
