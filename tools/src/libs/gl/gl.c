
/* 
 * FUNCTION:
 * Emulate GL library
 *
 * HISTORY:
 * Linas Vepstas February 1994
 */

#include "gl.h"
#include <X11/Xlib.h>

/* ------------------------------------------------------------ */
/* defines */

Display *dpy = NULL;

typedef struct GL_gc {
   struct GL_gc *next;
   int wid;

   Display *dpy;
   Window win;
   GC gc;

   /* window dimensions */
   int width;
   int height;

   /* current matrix */
   float offset_x;
   float offset_y;
   float scale_x;
   float scale_y;

} GL_gc_t;

int nextid =0;

GL_gc_t *win_list = NULL;
GL_gc_t *curr_gc;

int width = 100;
int height = 100;

/* ------------------------------------------------------------ */
/* emulate the GL prefsize routine */

void prefsize (w, h)
int w, h;
{
   width = w;
   height = h;
}

/* ------------------------------------------------------------ */

int
winopen () 
{
   Window win;
   GC gc;
   XGCValues gcv;
   XSetWindowAttributes wa;
   GL_gc_t *GLwin;

   if (!dpy) dpy = XOpenDisplay ("unix:0");

   win = XCreateSimpleWindow (dpy, RootWindow (dpy, 0), 
                             10, 10, width, height, 1, 0, 1);

   XMapWindow (dpy, win);
   
   gcv.foreground = 0;
   gcv.background = 1;
   gc = XCreateGC (dpy, win, GCForeground |GCBackground, gcv);

   wa.event_mask = ExposureMask;
   XChangeWindowAttributes (dpy, win, CWEventMask, &wa);
   
   GLwin = (GL_gc_t *) malloc (sizeof (GL_gc_t));

   GLwin->next = win_list;
   win_list = GLwin;

   /* copy important stuff into the gc */
   GLwin->dpy = dpy;
   GLwin->win = win;
   GLwin->gc = gc;
   GLwin->width = width;
   GLwin->height = height;

   /* initialize the pipe */
   GLwin -> scale_x = 1.0;
   GLwin -> scale_y = 1.0;

   nextid ++;
   GLwin->wid = nextid;

   curr_gc = GLwin;

   gsync ();

   return (GLwin->wid);
}

/* ------------------------------------------------------------ */
/* emulate the GL winset routine */

void
winset (id)
int id;
{
   GL_gc_t *GLwin;

   for (GLwin=win_list; GLwin !=0x0; GLwin=GLwin->next) {
      if (GLwin -> wid = id) {
         break;
      }
   }

   curr_gc = GLwin;
}

/* ------------------------------------------------------------ */
/* emulate orhto2 */

void ortho2 (xlow, xhi, ylow, yhi)
float xlow, xhi, ylow, yhi;
{
   curr_gc->scale_x = ((float) curr_gc->width) / (xhi-xlow);
   curr_gc->scale_y = ((float) curr_gc->height) / (yhi-ylow);
   curr_gc->offset_x = xlow;
   curr_gc->offset_y = ylow;
}
  
/* ------------------------------------------------------------ */
/* emulate the v2 routines */

void v2d (x, y) 
double x, y;
{
   int wx, wy;

   wx = curr_gc->scale_x * (x - curr_gc->offset_x);
   wy = curr_gc->scale_y * (y - curr_gc->offset_y);

   XDrawPoint (curr_gc->dpy, curr_gc->win, curr_gc->gc, 
               wx, curr_gc->height - wy);
}

/* ------------------------------------------------------------ */
/* emulate the v2 routines */

void v2f (x, y) 
float x, y;
{
   v2d ((double) x, (double) y);
}

/* ------------------------------------------------------------ */
/* emulate the v2 routines */

void v2i (x, y) 
int x, y;
{
   v2d ((double) x, (double) y);
}

/* ------------------------------------------------------------ */
/* emulate the GL clear routine */

void clear ()
{
   XClearWindow (curr_gc->dpy, curr_gc->win);
}

/* ------------------------------------------------------------ */

void gsync ()
{
   XSync (curr_gc->dpy, False);
}

/* ------------------------------------------------------------ */

/* ---------------------- END OF FILE ------------------------- */
