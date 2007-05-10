/*
 * eps.c 
 * Encapsulated PostScript utilities
 *
 * Linas Vepstas April 2007
 */

#include <math.h>
#include <stdio.h>
#include "eps.h"

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

void eps_setup_basic_linstyles (void)
{
	printf ("0.0100000 slw\n");
	printf ("[] 0 sd\n");
	printf ("[] 0 sd\n");
	printf ("0 slc\n");
	printf ("0.000000 0.000000 0.000000 srgb\n");
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

void eps_set_color (int r, int g, int b)
{
	printf ("%f %f %f srgb\n", r/255.0, g/255.0, b/255.0);
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

/* ======================= END OF FILE =================== */
