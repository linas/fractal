/*
 * brat.h
 *
 * FUNCTION:
 * Explore Hausdorf measure of mandelbrot set.
 * And other stuff.
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 */

void MakeHisto (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
	double 	renorm);

void 
MakeHistoCB (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm,
	double   (*cb)(double, double, int));


#define MAKE_HISTO(cb)  \
void MakeHisto (        \
   float  	*glob,      \
   int 		sizex,      \
   int 		sizey,      \
   double	re_center,  \
   double	im_center,  \
   double	width,      \
   double	height,     \
   int		itermax,    \
	double 	renorm)     \
{                       \
   MakeHistoCB (glob, sizex, sizey, re_center, im_center,  \
       width, height, itermax, renorm, cb);                \
}
