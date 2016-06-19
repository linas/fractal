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

/**
 * MakeHisto -- callback for making a generic scatterplot.
 *
 * Implement this callback to draw a generic scatterplot. Upon return,
 * the callback should fill in the 2D array with values.
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

/**
 * MakeHeightCB - callback for making a plain height-map.
 *
 * To graph a simple height map, implement this callback, and then use
 * DECL_MAKE_HEIGHT() to run it.
 *
 * The callback should return a single real number, given, as input,
 * a fixed point (x,y) on the 2D plane.
 */
typedef double MakeHeightCB
	(double x, double y, int itermax, double param);

void
MakeHeightWrap (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm,
	MakeHeightCB cb);


#define DECL_MAKE_HEIGHT(cb)  \
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
   MakeHeightWrap (glob, sizex, sizey, re_center, im_center,  \
       width, height, itermax, renorm, cb);                  \
}

/**
 * MakeBifurCB- Bifurcation diagram callback, does one row at a time.
 *
 * To draw a bifurcation diagram, implement this callback,
 * and then declare DECL_MAKE_BIFUR. The callback will be called
 * with a steadily-incremented y_parameter each time. For each
 * y_parameter, the callback should compute a 1D distribution, and
 * return that distribution in "array".
 */

typedef void
MakeBifurCB (
	float *array,
	int array_size,
	double x_center,
	double x_width,
	double y_parameter,
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
