
/*
 * globals
 */

#define I 6                     /* icon is 1/Ith - even I is best */

#define IWIDTH ((width+I-1)/I)  /* size of icon */
#define IHEIGHT ((height+I-1)/I)

extern float gcp;           /* gamma correction parameter */

extern int height, width;   /* main window dimensions */
extern char *dpy_name;      /* X server to connect to */

extern char title[200];     /* window title */
extern char table[256];     /* map from dither indices to X pixel values */

extern Display *dpy;        /* display we're talking to */
extern Window win, iwin;    /* main, icon windows */
extern Pixmap pixmap, ipixmap;     /* pixmaps for the image and for the icon */
extern GC gc;                      /* a graphics conntext to use */

