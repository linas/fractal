
/*
 * Copyright (c) 1987,1989 IBM Corporation
 * Author: Bruce Lucas (with help from Daniel Hu, Mickey Coggins,
 * and Bob Pearson for the X11 conversion)
 *
 * Revision History:
 *
 */

#include <X11/Xlib.h>
#include <stdio.h>
#include "debug.h"
#include "extern.h"
/*
 * globals
 */

float gcp = 0.5;                /* gamma correction parameter */

int height=0, width=0;          /* main window dimensions */
char *dpy_name = NULL;          /* X server to connect to */

char title[200];                /* window title */
char table[256];                /* map from dither indices to X pixel values */
int progress;                   /* rows availiable so far */

Display *dpy;                   /* display we're talking to */
Window win, iwin;               /* main, icon windows */
Pixmap pixmap, ipixmap;         /* pixmaps for the image and for the icon */
GC gc;                          /* a graphics context to use */

