/*
 * tree.c 
 * draw a binary tree on the hyperbolic disk, using postscript
 *
 * Linas Vepstas April 2007
 */

void eps_print_prolog (void)
{
	printf ("%!PS-Adobe-3.0 EPSF-3.0\n");
	printf ("%!PS-Adobe-2.0 EPSF-2.0\n");
	printf ("%%Title: blah blah\n");
	printf ("%%Creator: fractal/misc/hyperbolic/tree.c\n");
	printf ("%%CreationDate: Fri Apr 27 22:11:59 2007\n");
	printf ("%%For: linas\n");
	printf ("%%Orientation: Portrait\n");
	printf ("%%Magnification: 1.0000\n");
	printf ("%%BoundingBox: 0 0 192 271\n");
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

%%28.346000 -28.346000 scale
%%-4.280348 -14.619652 translate


0.100000 slw
[] 0 sd
[] 0 sd
0 slc
0.000000 0.000000 0.000000 srgb
n 4.350000 5.150000 m 10.950000 14.550000 l s
showpage
