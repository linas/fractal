/*
 * tree.c 
 * draw a binary tree on the hyperbolic disk, using postscript
 *
 * Linas Vepstas April 2007
 */

#include "cplex.h"
#include <math.h>
#include <stdio.h>

void eps_print_prolog (void)
{
	printf ("%!PS-Adobe-3.0 EPSF-3.0\n");
	// printf ("%!PS-Adobe-2.0 EPSF-2.0\n");
	printf ("%%Title: blah blah\n");
	printf ("%%Creator: fractal/misc/hyperbolic/tree.c\n");
	printf ("%%CreationDate: Fri Apr 27 22:11:59 2007\n");
	printf ("%%For: linas\n");
	printf ("%%Orientation: Portrait\n");
	printf ("%%Magnification: 1.0000\n");
	printf ("%%BoundingBox: 0 0 220 220\n");
	printf ("%%EndComments\n");
	printf ("%%BeginProlog\n");
	printf ("/cp {closepath} bind def\n");
	printf ("/c {curveto} bind def\n");
	printf ("/f {fill} bind def\n");
	printf ("/a {arc} bind def\n");
	printf ("/ef {eofill} bind def\n");
	printf ("/ex {exch} bind def\n");
	printf ("/gr {grestore} bind def\n");
	printf ("/gs {gsave} bind def\n");
	printf ("/sa {save} bind def\n");
	printf ("/rs {restore} bind def\n");
	printf ("/l {lineto} bind def\n");
	printf ("/m {moveto} bind def\n");
	printf ("/rm {rmoveto} bind def\n");
	printf ("/n {newpath} bind def\n");
	printf ("/s {stroke} bind def\n");
	printf ("/sh {show} bind def\n");
	printf ("/slc {setlinecap} bind def\n");
	printf ("/slj {setlinejoin} bind def\n");
	printf ("/slw {setlinewidth} bind def\n");
	printf ("/srgb {setrgbcolor} bind def\n");
	printf ("/rot {rotate} bind def\n");
	printf ("/sc {scale} bind def\n");
	printf ("/sd {setdash} bind def\n");
	printf ("/ff {findfont} bind def\n");
	printf ("/sf {setfont} bind def\n");
	printf ("/scf {scalefont} bind def\n");
	printf ("/sw {stringwidth pop} bind def\n");
	printf ("/tr {translate} bind def\n");
	printf ("%%EndProlog\n");
	printf ("\n");
}

void eps_setup_misc (void)
{
	printf ("0.0100000 slw\n");
	printf ("[] 0 sd\n");
	printf ("[] 0 sd\n");
	printf ("0 slc\n");
	printf ("0.000000 0.000000 0.000000 srgb\n");
	printf ("110.0 110.0 translate\n");
	printf ("100.0 -100.0 scale\n");
}

void eps_set_color_red (void)
{
	printf ("1.000000 0.000000 0.000000 srgb\n");
}

void eps_set_color_green (void)
{
	printf ("0.000000 1.000000 0.000000 srgb\n");
}

void eps_set_color_blue (void)
{
	printf ("0.000000 0.000000 1.000000 srgb\n");
}

/* =============================================== */

/* draw a circle of unit radius about the origin */
void eps_draw_circle(void)
{
	int i;
	double si, co, dsi, dco;
	double theta = 2.0*M_PI/360.0;
	dsi = sin (theta);
	dco = cos (theta);
	si = dsi;
	co = dco;

	// for example,draw one line:
	// n 4.350000 5.150000 m 10.950000 14.550000 l s
	printf ("n 0.0 1.0 m ");
	for (i=0; i<=360; i++)
	{
		printf("%f %f l \n", si,co);
		double tmp = si*dco + co*dsi;
		co = co*dco - si*dsi;
		si = tmp;
	}
	printf (" s\n");
}

/* draw three-pointed stick figure */
void draw_tristar (void)
{
	printf ("n 0.0 0.0 m -0.5 0.0 l s\n");
	printf ("n 0.0 0.0 m 0.25 0.433012702 l s\n");
	printf ("n 0.0 0.0 m 0.25 -0.433012702 l s\n");
}

/* ==================================================== */

/* fractional linear transform */
typedef struct {
	cplex a,b,c,d;
} mobius_t;

/* perform a multiplicative scaling of the thing */
mobius_t mobius_scale(mobius_t m, const cplex z)
{
	m.a = cplex_mul (m.a, z);
	m.b = cplex_mul (m.b, z);
	return m;
}

/* prooduct of mobius transforms */
mobius_t mobius_mul(const mobius_t l, const mobius_t r)
{
	mobius_t m;

	m.a = cplex_add (cplex_mul (l.a, r.a), cplex_mul(l.b, r.c));
	m.c = cplex_add (cplex_mul (l.c, r.a), cplex_mul(l.d, r.c));
	m.b = cplex_add (cplex_mul (l.a, r.b), cplex_mul(l.b, r.d));
	m.d = cplex_add (cplex_mul (l.c, r.b), cplex_mul(l.d, r.d));

	return m;
}

/* apply mobius xform to z */
cplex mobius_xform (const mobius_t m, const cplex z)
{
	cplex numer = cplex_mul(m.a, z);
	numer = cplex_add (numer, m.b);
	cplex deno = cplex_mul(m.c, z);
	deno = cplex_add (deno, m.d);
	cplex w = cplex_div (numer,deno);
	return w;
}

/* return mobius xform that recenters the disk at w */
mobius_t disk_center (cplex w)
{

	cplex nu = w;
	nu.re += 1.0;
	cplex de = w;
	de.re -= 1.0;

	cplex z = cplex_div (nu,de);
	z = cplex_neg (z);
	z = cplex_times_i (z);
	cplex zb = cplex_conj (z);

	cplex mi = cplex_set(0.0,-1.0);

	mobius_t m;
	m.a = cplex_sub (mi, z);
	m.b = cplex_add (mi, z);
	m.c = cplex_sub (mi, zb);
	m.d = cplex_add (mi, zb);
	return m;
}

void show_mobius(mobius_t m)
{
	printf ("a=%f +i%f    b=%f+i%f\n", m.a.re, m.a.im, m.b.re, m.b.im);
	printf ("c=%f +i%f    d=%f+i%f\n", m.c.re, m.c.im, m.d.re, m.d.im);
}

mobius_t go_to_fork_tip(double sign)
{
	cplex c = cplex_set (-0.5, 0.0);
	mobius_t m = disk_center (c);
	cplex rot = cplex_exp_itheta (sign*2.0*M_PI/6.0);
	m = mobius_scale (m,rot);
	return m;
}

void draw_fork(mobius_t m, int level)
{
	if (level == 0) return;
	level--;

	cplex za, zb;

	za = cplex_set (0.0, 0.0);
	za = mobius_xform (m,za);

	zb = cplex_set (0.25, 0.433012702);
	zb = mobius_xform (m,zb);
	printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);

	zb = cplex_set (0.25, -0.433012702);
	zb = mobius_xform (m,zb);
	printf ("n %f %f m %f %f l s\n", za.re, za.im,zb.re, zb.im);

	mobius_t tip = go_to_fork_tip(+1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);

	tip = go_to_fork_tip(-1.0);
	tip = mobius_mul (m, tip);
	draw_fork (tip, level);
}


void draw(void)
{
	mobius_t m;
	int level=11;

eps_set_color_red();
	m = go_to_fork_tip(+1.0);
	draw_fork (m, level);

eps_set_color_green();
	m = go_to_fork_tip(-1.0);
	draw_fork (m, level);

eps_set_color_blue();
	m = go_to_fork_tip(-3.0);
	draw_fork (m, level);
}

/* ==================================================== */

main () 
{
	eps_print_prolog();
	eps_setup_misc();
	eps_draw_circle();
	draw_tristar();

	draw();

	printf ("showpage\n");
}
