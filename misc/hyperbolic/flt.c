/*
 * tree.c 
 * draw a binary tree on the hyperbolic disk, using postscript
 *
 * Linas Vepstas April 2007
 */

#include "cplex.h"
#include <math.h>
#include <stdio.h>

void eps_print_prolog (int width, int height)
{
	printf ("%%!PS-Adobe-3.0 EPSF-3.0\n");
	// printf ("%!PS-Adobe-2.0 EPSF-2.0\n");
	printf ("%%%%Title: blah blah\n");
	printf ("%%%%Creator: fractal/misc/hyperbolic/tree.c\n");
	printf ("%%%%CreationDate: Fri Apr 27 22:11:59 2007\n");
	printf ("%%%%For: linas\n");
	printf ("%%%%Orientation: Portrait\n");
	printf ("%%%%Magnification: 1.0000\n");
	// printf ("%%%%BoundingBox: 0 0 220 220\n");
	printf ("%%%%BoundingBox: 0 0 %d %d\n", width, height);
	printf ("%%%%EndComments\n");
	printf ("%%%%BeginProlog\n");
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
	printf ("%%%%EndProlog\n");
	printf ("\n");
	printf ("%d %d translate\n", width/2, height/2);
	printf ("%d %d scale\n", (45*width)/100, (45*height)/100);
}

void eps_setup_disk (void)
{
	printf ("0.0100000 slw\n");
	printf ("[] 0 sd\n");
	printf ("[] 0 sd\n");
	printf ("0 slc\n");
	printf ("0.000000 0.000000 0.000000 srgb\n");
}
void eps_setup_plane (void)
{
	printf ("0.0100000 slw\n");
	printf ("[] 0 sd\n");
	printf ("[] 0 sd\n");
	printf ("0 slc\n");
	printf ("0.000000 0.000000 0.000000 srgb\n");
	printf ("0.0 -0.8 translate\n");
	printf ("0.4 1.6 scale\n");
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
	int nseg = 120;
	double si, co, dsi, dco;
	double theta = 2.0*M_PI/((double) nseg);
	dsi = sin (theta);
	dco = cos (theta);
	si = dsi;
	co = dco;

	// for example,draw one line:
	// n 4.350000 5.150000 m 10.950000 14.550000 l s
	printf ("n 0.0 1.0 m ");
	for (i=0; i<=nseg; i++)
	{
		printf("%f %f l \n", si,co);
		double tmp = si*dco + co*dsi;
		co = co*dco - si*dsi;
		si = tmp;
	}
	printf (" s\n");
}

void eps_draw_lineseg (double fx, double fy, double tx, double ty)
{
	printf ("n %f %f m %f %f l s\n", fx, fy, tx, ty);
}

/* ==================================================== */

/* fractional linear transform */
typedef struct {
	cplex a,b,c,d;
} mobius_t;

/* create an element of sl(2,z). Its up to user to ensure
 * that ad-bc=1 */
mobius_t mobius_set(int a, int b, int c, int d)
{
	mobius_t m;
	m.a = cplex_set (a, 0);
	m.b = cplex_set (b, 0);
	m.c = cplex_set (c, 0);
	m.d = cplex_set (d, 0);
	return m;
}

/* return a transformation that simply rotates by theta radians */
mobius_t mobius_rotate (double theta)
{
	mobius_t m;
	m.a = cplex_exp_itheta (theta);
	m.b = cplex_zero();
	m.c = cplex_zero();
	m.d = cplex_one();
	return m;
}

/* perform a multiplicative scaling of the thing */
mobius_t mobius_scale(mobius_t m, const cplex z)
{
	m.a = cplex_mul (m.a, z);
	m.b = cplex_mul (m.b, z);
	return m;
}

/* product of mobius transforms (just a matrix multiply) */
mobius_t mobius_mul(const mobius_t l, const mobius_t r)
{
	mobius_t m;

	m.a = cplex_add (cplex_mul (l.a, r.a), cplex_mul(l.b, r.c));
	m.c = cplex_add (cplex_mul (l.c, r.a), cplex_mul(l.d, r.c));
	m.b = cplex_add (cplex_mul (l.a, r.b), cplex_mul(l.b, r.d));
	m.d = cplex_add (cplex_mul (l.c, r.b), cplex_mul(l.d, r.d));

	return m;
}

/* apply mobius xform to z, return (az+b)/(cz+d) */
cplex mobius_xform (const mobius_t m, const cplex z)
{
	cplex numer = cplex_mul(m.a, z);
	numer = cplex_add (numer, m.b);
	cplex deno = cplex_mul(m.c, z);
	deno = cplex_add (deno, m.d);
	cplex w = cplex_div (numer,deno);
	return w;
}

/* Return mobius xform that recenters the disk at w */
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

/* Return a mobius transform that maps the point z=+1 to the 
 * center of the Poincare disk */
mobius_t to_half_plane_xform(void)
{
	mobius_t m;
	m.a = cplex_set (0.0, -1.0);
	m.b = cplex_set (0.0, -1.0);
	m.c = cplex_set (1.0, 0.0);
	m.d = cplex_set (-1.0, 0.0);
	return m;
}

/* convert a half-plane xform to a disk xform */
mobius_t to_disk(mobius_t m)
{
	cplex apd = cplex_add (m.a, m.d);
	cplex amd = cplex_sub (m.a, m.d);
	cplex bpc = cplex_times_i (cplex_add (m.b, m.c));
	cplex bmc = cplex_times_i (cplex_sub (m.b, m.c));

	mobius_t n;
	n.a = cplex_scale (0.5, cplex_add (apd, bmc));
	n.b = cplex_scale (0.5, cplex_sub (amd, bpc));
	n.c = cplex_scale (0.5, cplex_add (amd, bpc));
	n.d = cplex_scale (0.5, cplex_sub (apd, bmc));
	
	return n;
}

/* convert a disk xform to a plane xform */
mobius_t to_half_plane(mobius_t m)
{
	cplex apd = cplex_add (m.a, m.d);
	cplex amd = cplex_sub (m.a, m.d);
	cplex bpc = cplex_add (m.b, m.c);
	cplex bmc = cplex_sub (m.b, m.c);

	mobius_t n;
	n.a = cplex_scale (0.5, cplex_add (apd, bpc));
	n.b = cplex_scale (0.5, cplex_times_i (cplex_sub (amd, bmc)));
	n.b = cplex_neg (n.b);
	n.c = cplex_scale (0.5, cplex_times_i (cplex_add (amd, bmc)));
	n.d = cplex_scale (0.5, cplex_sub (apd, bpc));
	
	return n;
}

void show_mobius(mobius_t m)
{
	printf ("a=%f +i%f    b=%f+i%f\n", m.a.re, m.a.im, m.b.re, m.b.im);
	printf ("c=%f +i%f    d=%f+i%f\n", m.c.re, m.c.im, m.d.re, m.d.im);
}

/* ==================================================== */

/* Draw straight line segment */
void draw_seg(mobius_t m, cplex zf, cplex zt)
{
	cplex za = mobius_xform (m, zf);
	cplex zb = mobius_xform (m, zt);
	eps_draw_lineseg (za.re, za.im, zb.re, zb.im);
}

/* draw line segment in the klein model */
void draw_klein_seg(mobius_t m, cplex zf, cplex zt)
{
	cplex za = mobius_xform (m, zf);
	cplex zb = mobius_xform (m, zt);

	double ma = 2.0/(1.0 + za.re * za.re + za.im * za.im);
	double mb = 2.0/(1.0 + zb.re * zb.re + zb.im * zb.im);
	za = cplex_scale (ma,za);
	zb = cplex_scale (mb,zb);
	eps_draw_lineseg (za.re, za.im, zb.re, zb.im);
}

/* draw statically-tesselated arc */
void draw_tesselated_arc(mobius_t m, cplex zf, cplex zt)
{
	int i;
	// int nseg=20;
	int nseg=5;

	cplex zstart = zf;
	cplex zdelta = cplex_scale((1.0/(double)nseg), cplex_sub (zt,zf));
	cplex za = mobius_xform (m, zstart);

	for (i=0; i<nseg; i++)
	{
		cplex zend = cplex_add (zstart, zdelta);
		cplex zb = mobius_xform (m, zend);
		eps_draw_lineseg (za.re, za.im, zb.re, zb.im);
		zstart = zend;
		za = zb;
	}
}

#if 0
/* draw dynamically-tesselated arc */
void draw_arc(mobius_t m, cplex zf, cplex zt)
{
	draw_klein_seg(m,zf,zt);
}
#endif

void draw_arc(mobius_t m, cplex zf, cplex zt)
{
	int i;
	int nseg=20;

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
			eps_draw_lineseg (zs.re, zs.im, zb.re, zb.im);
			zs = zb;
		}
		zstart = zend;
		za = zb;
	}

	double dist = cplex_dist (zs, zb);
	if (dist > 0.001) 
		eps_draw_lineseg (zs.re, zs.im, zb.re, zb.im);
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
	int n = atoi (argv[1]);

	eps_print_prolog(220,220);
	// eps_print_prolog(400,100);
	eps_setup_disk();
	// eps_setup_plane();
	eps_draw_circle();

	draw(n);

	printf ("showpage\n");
}
