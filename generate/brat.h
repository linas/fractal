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
 * godawful hackery - October 2016
 */

#include <string.h>

/**
 * MakeHisto -- callback for making a generic scatterplot.
 *
 * Implement this callback to draw a generic scatterplot. Upon return,
 * the callback should fill in the 2D array with values.
 */
void MakeHisto (
   char     *name,     /* value of argv[0] */
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

#define MAX_NUM_NAMES 8000
extern int num_names;
extern MakeHeightCB* callbacks[MAX_NUM_NAMES];
extern const char* main_names[MAX_NUM_NAMES];

// Declare a named callback.
extern void decl_height(const char* name, MakeHeightCB* cb);

#define DECL_HEIGHT(name,cb) \
	main_names[num_names] = name; \
	callbacks[num_names] = &cb; \
	num_names ++;

#define MAKE_HEIGHT     \
void MakeHisto (        \
   char     *name,      \
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
	int i;                                                        \
	for (i=0; i<num_names; i++) {                                 \
      if (0 == strcmp(name, main_names[i])) break;               \
	}                                                             \
	if (num_names <= i) {                                         \
		fprintf(stderr, "Unable to find callback for %s\n", name); \
		exit(1);                                                   \
	}                                                             \
   MakeHeightWrap (glob, sizex, sizey, re_center, im_center,     \
       width, height, itermax, renorm, *callbacks[i]);           \
}



#define DECL_MAKE_HEIGHT(cb)  \
void MakeHisto (        \
   char     *name,      \
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
   char     *name,      \
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
