
/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas Januery 16 1994
 */

#include <X11/Xlib.h>

#include "Farey.h"
#include <stdio.h>
#include <math.h>

/* ------------------------------------------------------------ */
/* returns the farey-number" */

double f_sub_x (f, nume, deno, z) 
struct Farey *f;
int nume, deno;
double z;
{
   double r;
   RatioToContinuedFraction (f, nume, deno);
   r = ContinuedFractionToEFarey (f, z);
   return (r);
}

/* ------------------------------------------------------------ */
/* returns the "z-number" */

double r_sub_x (f, nume, deno, z) 
struct Farey *f;
int nume, deno;
double z;
{
   double r;
   RatioToContinuedFraction (f, nume, deno);
   r = ContinuedFractionToZReal (f, z);
   return (r);
}

/* ------------------------------------------------------------ */
/* returns a symmetric version of the z-number.
 * (the gaps are symmetrically arranged)
 */

double sym_sub_x (f, nume, deno, z) 
struct Farey *f;
int nume, deno;
double z;
{
   double r;
   r = r_sub_x (f, nume, deno, z);
   r -= r_sub_x (f, deno-nume, deno, z);
   r +=1.0;
   r *= 0.5;
   return (r);
}

/* ------------------------------------------------------------ */
/* Returns a symmetric z-number, 
 * where the gap at x=0.5 is directly given by the argument z
 */

double sym_gap_x (f, nume, deno, z) 
struct Farey *f;
int nume, deno;
double z;
{
   double r;
   r = sym_sub_x (f, nume, deno, (1.0 - 2.0*z) /(1.0 + 2.0*z));
   return (r);
}

/* ------------------------------------------------------------ */

Display *dpy = NULL;
Window win;
GC gc;

typedef struct GL_gc {
   struct GL_gc *next;
   int wid;
   Window w;
   GC g;
} GL_gc_t;

int nextid =0;

GL_gc_t *win_list = NULL;

/* ------------------------------------------------------------ */

winopen () 
{
   XGCValues gcv;
   XSetWindowAttributes wa;
   GL_gc_t *GLwin;

   if (!dpy) dpy = XOpenDisplay ("unix:0");

   win = XCreateSimpleWindow (dpy, RootWindow (dpy, 0), 
                             10, 10, 500, 500, 1, 0, 1);

   XMapWindow (dpy, win);
   
   gcv.foreground = 0;
   gcv.background = 1;
   gc = XCreateGC (dpy, win, GCForeground |GCBackground, gcv);

   wa.event_mask = ExposureMask;
   XChangeWindowAttributes (dpy, win, CWEventMask, &wa);
   
   GLwin = (GL_gc_t *) malloc (sizeof (GL_gc_t));

   GLwin->next = win_list;
   win_list = GLwin;
   GLwin->w = win;
   GLwin->g = gc;

   nextid ++;
   GLwin->wid = nextid;

   XSync (dpy, False);

   return (GLwin->wid);
}

/* ------------------------------------------------------------ */

void
winset (id)
int id;
{
   GL_gc_t *GLwin;

   for (GLwin=win_list; GLwin !=0x0; GLwin=GLwin->next) {
      if (GLwin -> wid = id) {
         win = GLwin->w;
         gc = GLwin->g;
         break;
      }
   }
}

/* ------------------------------------------------------------ */

void v2i (x, y) 
int x, y;
{
   XDrawPoint (dpy, win, gc, x, 500 - y);
}

/* ------------------------------------------------------------ */

struct Farey *f;

redraw (n, z) 
int n;
double z;
{
   double x, y, t;
   int i;
   int nume, deno;
   int ix, iy;
   double r, theta;

   XClearWindow (dpy, win);

   deno = n* 6899 +1;
   for (i=0; i<n; i++){
      nume = 6899*i;


/* ----------- */

#ifdef INTERESTING
      t = - log (0.5 * (1.0+z));
      x = f_sub_x (f, nume, deno, t);
      y = sym_gap_x (f, nume, deno, -z);
      /*
      t = - log (0.5 * (1.0-z));
      y = f_sub_x (f, nume, deno, t);
      */

      theta = 2.0 * M_PI * x;
      r = y;
      r = sin (M_PI * r);

      x = r * cos (theta);
      y = r * sin (theta);

      ix = 250 * (x+1.0);
      iy = 250 * (y+1.0);
      v2i (ix, iy);
#endif 

      t = - log (0.5 * (1.0+z));
      x = f_sub_x (f, nume, deno, t);
      y = sym_gap_x (f, nume, deno, -z);

      theta = M_PI * x;
      r = y;
      r -=  ((double) nume) / ((double) deno);
      r *= 2.0;

      x = r * cos (theta);
      y = r * sin (theta);

      ix = 250 * (x+1.0);
      iy = 250 * (y+1.0);
      v2i (ix, iy);
   }

   XSync (dpy, False);
}

/* ------------------------------------------------------------ */

main (argc, argv)
int argc;
char *argv[];
{
   double r, x, y, z, t;
   int i, n;
   int nume, deno;
   int ix, iy;

   if (argc <3) {
      printf ("Usage: %s <number of terms> <base> \n", argv[0]);
      exit (1);
   }

   winopen ();

   f = CreateFarey();

   n = atoi (argv[1]);
   z = atof (argv[2]);

   while (TRUE) {
      XEvent ev;
      XNextEvent (dpy, &ev);

      switch (ev.type) {
         case ConfigureNotify:
         case Expose:
            printf (" Expose \n");
            redraw (n, z);
            break;
         default:
            break;
      }
   }

}

/* ---------------------- END OF FILE ------------------------- */
