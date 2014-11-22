
/*
 * FUNCTION:
 * GL library emulation include header file
 *
 * HISTORY:
 * Linas Vepstas February 1994
 */

#include <X11/Xlib.h>

extern Display *dpy;

extern void clear ();
extern void gsync ();
extern void ortho2 ();
extern void prefsize ();
extern void v2d ();
extern void v2f ();
extern void v2i ();
extern int winopen ();
extern void winset ();

/* ------------------- END OF FILE ------------------------ */


