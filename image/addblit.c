
/*
 * Copyright (c) 1987,1989 IBM Corporation
 * Author: Bruce Lucas (with help from Daniel Hu, Mickey Coggins,
 * and Bob Pearson for the X11 conversion)

Revision History:
15 June 1989 by Linas Vepstas -- insert an XInstallColormap 

24 October 1989 by Linas Vepstas -- if all colormap entries already
                   taken, search for color closest to desired color 

11 January 1990 by Linas Vepstas -- above change defined out,
                   appearently this was an X server bug.  Since some
                   servers may still have this bug, the code was not
                   totally thrown away.

 */

#include <math.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include "debug.h"
#include "extern.h"

typedef struct rgb {
	char r, g, b;
        };

/*-------------------------------------------------------------------*/

/*
 * Initialize X, create windows, compute color table
 */

initialize_RGB_window (argc, argv, pix_width, pix_height)
int	argc;
char 	**argv;
unsigned int	pix_width, pix_height;	/* size of window/pixmap */
{
    char comp[256];            /* gamma correction table */
    struct rgb map[256];       /* map from dither indices to colors */
    Colormap cmap;             /* a colormap to use */
    int black, white;          /* black, white pixel values */
    Window root;
    XWMHints wm_hints;
    XSizeHints size_hints;
    XColor cur_map[256];
    int i, j, n, bomb;

    width = pix_width;
    height = pix_height;

    /* open display etc. */
    dpy = XOpenDisplay(dpy_name);
    if (!dpy)
        Exit("can't open display");
    root = RootWindow(dpy, DefaultScreen(dpy));
    cmap = DefaultColormap(dpy, DefaultScreen(dpy));
    XInstallColormap (dpy, cmap);    /*  Server problem ???? */
    black = BlackPixel(dpy, DefaultScreen(dpy));
    white = WhitePixel(dpy, DefaultScreen(dpy));

    gc = XCreateGC(dpy, root, 0, NULL);
    XSetGraphicsExposures(dpy, gc, 0);

    /* iwin */
    iwin = XCreateSimpleWindow(dpy, root, 0, 0, IWIDTH, IHEIGHT, 1,
        black, white);
    size_hints.flags = USSize;
    XSelectInput(dpy, iwin, ExposureMask);

    /* win */
    win = XCreateSimpleWindow(dpy, root, 0, 0, width, height, 3, black, white);
    size_hints.flags = PSize;
    size_hints.width = width;
    size_hints.height = height;
    XSetStandardProperties(dpy, win, title, title, NULL,
        argv, argc, &size_hints);
    wm_hints.flags = IconWindowHint;
    wm_hints.icon_window = iwin;
    XSetWMHints(dpy, win, &wm_hints);
    XSelectInput(dpy, win, ExposureMask|ResizeRedirectMask);
    XMapRaised(dpy, win);

    {
        /* Need to do this to get window resize to work */
        XSetWindowAttributes swin_attr; 

        swin_attr.override_redirect = TRUE;
        XChangeWindowAttributes (dpy, win, CWOverrideRedirect, &swin_attr);
    }

    /* fill in the color table */
    fprintf(stderr, "creating color table ...");
    n = ColorMap(map);
    /* Set up a compensating gamma correction map */
    for (i=0; i<256; i++)
        comp[i] = sqrt((i)/256.0)*256.0*gcp + i*(1.0-gcp);
   /* allocate as many colors as possible */
   bomb = FALSE;
   for (i=0; i<n; i++) {
        XColor color;
        color.red = comp[map[i].r] << 8;
        color.green = comp[map[i].g] << 8;
        color.blue = comp[map[i].b] << 8;
        TRACE(("request color %d %d %d ", color.red, color.blue, color.green));
        /* If XAllocColor finds an empty colormap slot, it fills it and
         * returns true.  If it cannot find an empty slot, it returns
         * the nearest slot and returns FALSE. Even if it returns false,
         * it should have done a "lower manhattan" search for best fit.
         */
        if (!XAllocColor(dpy, cmap, &color)) {
            TRACE(("XAllocColor bombed "));
            table[i] = -1;
            bomb = TRUE;
        } else {
            table[i] = color.pixel;
        }
        TRACE(("got color %d %d %d %d\n", color.pixel, color.red,
        color.blue, color.green));
    }

    /* if unable to satisfy all color requests */
    /* perform a least lower manhattan distance search */
    if (bomb) fprintf(stderr, "\nWarning: sharing non-optimal color table ...");
    fprintf(stderr, "\n");

    /* allocate server-side pixmaps */
    pixmap = XCreatePixmap(dpy, win, width, height, 8);
    ipixmap = XCreatePixmap(dpy, iwin, IWIDTH, IHEIGHT, 8);

}

/*-------------------------------------------------------------------*/

void add_block_to_window (block, block_width, block_height,
			x_offset, y_offset)
char	block [];
unsigned int	block_width;
unsigned int	block_height; 
unsigned int	x_offset, y_offset; /* where to start writing */

{
    char *ipixels, *ip;
    Visual *visual = DefaultVisual(dpy, DefaultScreen(dpy));
    XImage *image, *iimage;
    register int i, j;
    register char *p;

    image = XCreateImage(dpy, visual, 8, ZPixmap, 0, block,
            block_width, block_height, 8, 0);
    ip = ipixels = (char *) Malloc(IWIDTH*IHEIGHT);
    iimage = XCreateImage(dpy, visual, 8, ZPixmap, 0, ipixels,
            IWIDTH, IHEIGHT, 8, 0);

    /* copy from block to icon */
    for (j=0, p=block; j<block_height; j+=I, p+=I*block_width)
        for (i=0; i<block_width; i+=I)
            *ip++ = p[i];

    /* put our image data into the server-maintained pixmap */
    XPutImage(dpy, pixmap, gc, image, 0, 0, x_offset, y_offset, 
		block_width, block_height);
    /* copy from server-maintained pixmap into the window */
    XCopyArea(dpy, pixmap, win, gc, x_offset, y_offset, 
		block_width, block_height, x_offset, y_offset);

    /* put our image data into the server-maintained pixmap */
    XPutImage(dpy, ipixmap, gc, iimage, 0, 0, 0, 0, IWIDTH, IHEIGHT);

    free (ip);	/* free icon pixels */
    XDestroyImage (iimage); /* free the icon image data structure */
    XDestroyImage (image);	/* free image data structure */
}

/*-------------------------------------------------------------------*/
void fix_damage() {

    /* fix damage */
    while (XPending(dpy))
        handle_event();
}

/*-------------------------------------------------------------------*/
/*
 * X event loop
 */

io_error()
{
    exit(0);
}

event_loop()
{
    if (fork()) return;    /* detatch child process - parent returns */
    setpgrp();     /* so interrupts from terminal don't get to us */
    XSetIOErrorHandler(io_error);

    for (;;)
        handle_event();
}

/*-------------------------------------------------------------------*/

handle_event()
{
    XEvent event;
    int x, y, wid, hig;

    XNextEvent(dpy, &event);

    switch (event.type) {
    case Expose:
        x = event.xexpose.x;
        y = event.xexpose.y;
        wid = event.xexpose.width;
        hig = event.xexpose.height;
        if (event.xexpose.window==win) {
            XCopyArea(dpy, pixmap, win, gc, x, y, wid, hig, x, y);
            XFlush (dpy);

        } 
        /* else icon was exposed */
        else if (event.xexpose.window==iwin)
            XCopyArea(dpy, ipixmap, iwin, gc, x, y, wid, hig, x, y);

        break;

    default:
        break;
    }
}

/*-------------------------------------------------------------------*/

/*
 * Constants
 */

#define R 5                /* number of red shades */
#define G 9                /* number of green shades */
#define B 5                /* number of blue shades */

/*
 * MIX(r,g,b) calculates the color table index for r,g,b;
 * 0<=r<R, 0<=g<G, 0<=b<B.
 */

#define MIX(r,g,b) (((r)*G+(g))*B+(b))


/*
 * Fill in a color map
 */

ColorMap(map)
struct rgb map[256];
{
    register int r, g, b, i;
    for (r=0; r<R; r++) {
        for (g=0; g<G; g++) {
            for (b=0; b<B; b++) {
                i = MIX(r,g,b);
                map[i].r = r*255/(R-1);
                map[i].g = g*255/(G-1);
                map[i].b = b*255/(B-1);
            }
        }
    }
    return R*G*B;
}

/*-------------------------------------------------------------------*/

Malloc(n)
{
    int m = malloc(n);
    if (!m)
        Exit("can't malloc");
    return m;
}

/*-------------------------------------------------------------------*/

Exit(s)
char *s;
{
    fprintf(stderr, "%s\n", s);
    exit(1);
}

