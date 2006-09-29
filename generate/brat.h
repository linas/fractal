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

/** callback for making a generic scatterplot.
 *  Just implement this callback to draw a generic
 *  scatterplot.
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

/** callback for making a plain height-map. 
 *  To graph is simple height map i.e. a distinct real 
 *  for a point (x,y), implement this callback, 
 *  and then use DECL_MAKE_HEIGHT()
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

/* Bifurcation diagram callback, does one row at a time.
 * To draw a bifurcation diagram, implement this callback, 
 * and then declare DECL_MAKE_BIFUR */

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

