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
MakeHistoWrap (
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


#define DECL_MAKE_HISTO(cb)  \
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
   MakeHistoWrap (glob, sizex, sizey, re_center, im_center,  \
       width, height, itermax, renorm, cb);                  \
}

/* Bifurcation diagram callback, does one row at a time */
typedef void
MakeBifurCB (
	float *array, 
	int array_size, 
	double x_center,
	double x_width,
	double y_paramter, 
	int itermax,
	double renorm);

void
MakeBifurWrap (
	float    *glob,
	int      sizex,
	int      sizey,
	double   re_center,
	double   im_center,
	double   width,
	double   height,
	int      itermax,
	double   renorm,
	MakeBifurCB cb);

#define DECL_MAKE_BIFUR(cb)  \
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
   MakeBifurWrap (glob, sizex, sizey, re_center, im_center,  \
       width, height, itermax, renorm, cb);                  \
}

