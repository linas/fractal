#include  <X11/Xlib.h>
#include  <X11/Intrinsic.h>
#include  <X11/Xutil.h>


void blit_pixels_to_new_window 
			(argc, argv, my_pix_data,
				my_pix_width,
				my_pix_height)
int	 	argc;		/* passed from main */
char 		**argv;		/* passed from main */ 
char		my_pix_data [];	/* array contianing pixel data */
unsigned int 	my_pix_width, my_pix_height; 	/* pixmap size */

/*
 * This routine is very simple in conception.  It accepts eight-bit
 * pixels as input, opens a new window, puts the raw pixels into the
 * window (using the default colormap ???), and then enters an infinite
 * loop to handle the redraw events.  The inifinite loop then runs as a
 * detached process, so nothing new can be added to the window after it has
 * been opened.
 */

{
Display		*display;	/* Connection to X server */
int		screen;		/* cathode ray tube id */
Window		root_win;	/* root attributes */
Window 		child_win;
XSizeHints	my_size_hints;	/* hints for the window manager */
XEvent		event_buffer;	/* buffer for holding events */
GC		context;	/* graphics context */
XImage		*image;		/* this is my pixmap */
Pixmap		pixmap_id;	/* pixmap */

unsigned long 	border_color, ground_color;
int 		display_width, display_height, display_depth;
int		my_pix_padding;			/* pixmap junk */
int		bytes_per_line;	/* number of bytes in a scanline */

display = XOpenDisplay (NULL);
screen = DefaultScreen (display);	/* get CRT id number */ 
root_win = RootWindow (display, screen);	/* get default attributes */
context = XCreateGC (display, root_win, 0, NULL);

my_pix_padding = 32;	/* scanlines must start on word boundaries */
bytes_per_line = 0;	/* let X figure out the length of a scanline */
/* put my pixmap data into the client side X image data structure */
image = XCreateImage (display, NULL, 8, ZPixmap, 0, 
	my_pix_data,
	my_pix_width, my_pix_height, my_pix_padding,
	bytes_per_line);

/* tell server that start managing my pixmap */
pixmap_id = XCreatePixmap (display, root_win, my_pix_width,
	my_pix_height, 8);

/* copy from client to server */
(void) XPutImage (display, pixmap_id, context, image, 0,0, 0, 0,
	 	my_pix_width, my_pix_height);

/* free up the client side pixmap data area */
XDestroyImage (image);

/* set up window size hints for the window manager */
my_size_hints.flags = PSize;
my_size_hints.height = my_pix_width;
my_size_hints.width = my_pix_height;

child_win = XCreateSimpleWindow (display, root_win, 
	10, 10, my_pix_width, my_pix_height,
	0, BlackPixel(display, screen), WhitePixel (display, screen));

/* pass my size hints to the window manager, along with window
	and icon names */
(void) XSetStandardProperties (display, child_win, "Color Pixmap",
	"Pixmap", None, argv, argc, &my_size_hints);

/* copy the server copy of the pixmap into the window */
(void) XCopyArea (display, pixmap_id, child_win, context, 0,0,
		my_pix_width, my_pix_height, 0, 0);

/* set masks for event reporting */
(void) XSelectInput (display, child_win, ExposureMask);

/* get the window up onto the screen */
(void) XMapWindow (display, child_win);
(void) XSync (display, 0);	/* make it a squeaky clean map */

/* detach a new process */
if (fork()) return;

while (1) { 		/* begin event loop */
   XNextEvent (display, &event_buffer);
   if (event_buffer.type == Expose) {
      (void) XCopyArea (display, pixmap_id, child_win, context, 0,0,
		my_pix_width, my_pix_height, 0, 0);
   }
}

}


/*-------------------------------------------------------------------*/

void blit_pixmap_to_root 
			(argc, argv, my_pix_data,
				my_pix_width,
				my_pix_height,
				bytes_per_scanline)
int	 	argc;		/* passed from main */
char 		**argv;		/* passed from main */ 
char		my_pix_data [];	/* array contianing pixel data */
unsigned int 	my_pix_width, my_pix_height; 	/* pixmap size */
unsigned int	bytes_per_scanline;

/* This is the version to use with the Megapel -- it blits a one byte
* deep pixmap to the root window 				*/

{
Display		*display;	/* Connection to X server */
int		screen;		/* cathode ray tube id */
Window		root_win;	/* root attributes */
XImage		*image;
Visual		visual;
GC		gc;
Pixmap		pixmap;		/* pixmap */

display = XOpenDisplay (NULL);
screen = DefaultScreen (display);	/* get CRT id number */ 
root_win = RootWindow (display, screen);	/* get default attributes */
gc = XCreateGC (display, root_win, 0, NULL);

/* create client side data structure */
image = XCreateImage (display, NULL, 8, ZPixmap, 0, my_pix_data,
	my_pix_width, my_pix_height, 8, bytes_per_scanline); 

/* image = XCreateImage (display, &visual, 8, ZPixmap, 0, my_pix_data,
	my_pix_width, my_pix_height, 8, bytes_per_scanline); */

/* create server side data structure */
pixmap = XCreatePixmap (display, root_win, my_pix_width,
	my_pix_height, 8);

/* copy from client to server */
XPutImage (display, pixmap, gc, image, 0, 0, 0, 0, 
	my_pix_width, my_pix_height);

/* free client side data structures */
XDestroyImage  (image);

(void) XSetWindowBackgroundPixmap (display, root_win, pixmap); 
/* Use clear area to fake out window manager with an exposure event */
XClearArea (display, root_win, 0, 0, 1, 1, TRUE);
/* flush out the protocol buffer */
(void) XSync (display, 0);

XFreePixmap (display, pixmap);
XFreeGC (display, gc);
XCloseDisplay (display);
printf ("Done Blitting\n");

}
/*-------------------------------------------------------------------*/

void blit_bitmap_to_root 
			(argc, argv, my_pix_data,
				my_pix_width,
				my_pix_height,
				bytes_per_scanline)
int	 	argc;		/* passed from main */
char 		**argv;		/* passed from main */ 
char		my_pix_data [];	/* array contianing pixel data */
unsigned int 	my_pix_width, my_pix_height; 	/* pixmap size */
unsigned int	bytes_per_scanline;

/* This is the version to use with the APA -- it blits a one bit
* deep bitmap to the root window 				*/

{
Display		*display;	/* Connection to X server */
int		screen;		/* cathode ray tube id */
Window		root_win;	/* root attributes */
Pixmap		pixmap;		/* pixmap */

display = XOpenDisplay (NULL);
screen = DefaultScreen (display);	/* get CRT id number */ 
root_win = RootWindow (display, screen);	/* get default attributes */

/* allocate space for the pixmap */
pixmap = XCreatePixmapFromBitmapData
 	(display, root_win, 
	my_pix_data, 8*bytes_per_scanline, my_pix_height, 0, 1, 1); 

(void) XSetWindowBackgroundPixmap (display, root_win, pixmap); 
/* Use clear area to fake out window manager with an exposure event */
XClearArea (display, root_win, 0, 0, 1, 1, TRUE);
/* flush out the protocol buffer */
(void) XSync (display, 0);

XCloseDisplay (display);
printf ("Done Blitting\n");

}

#define	DUTY_CYCLE	0.5

/*-------------------------------------------------------------------*/

void add_float_data_to_pixmap 
			(data_array, data_len,
			data_x_offset, data_y_offset,
			pix_array, pix_width, pix_height)
float		data_array [];
unsigned int	data_len;		/* length of data_array */
unsigned int	data_x_offset, data_y_offset; /* where to start writing */
char		pix_array [];		/* bitmap to be written into */
unsigned int	pix_width, pix_height;	/* size of bitmap */

/* 
add_float_data_to_pixmap --- Linas Vepstas --- 3 May 1989

This routine takes the array of floating point numbers in data_array
(which should all be normalized -- i.e. run between zero and one)
and stuffs them into a 8 bit deep pixmap.  Each float represents one byte.
The location written to is given by the two offsets.
Writing is performed in the x (horizontal) direction.

The code has been optimized for human readability, not for fast
execution.  
*/
{
int		i;
unsigned int 	byte;
char		off;
unsigned int	byte_off, bytes_per_line;

byte_off = pix_width * data_y_offset + data_x_offset;
/* loop over x directions */
for (i=0; i<data_len; i++) {
   pix_array [i+byte_off] = (int) (255.0 * data_array [i]);
   } 
}

/*-------------------------------------------------------------------*/

void add_float_data_to_bitmap 
			(data_array, data_len,
			data_x_offset, data_y_offset,
			pix_array, pix_width, pix_height,
			bytes_per_line_return)
float		data_array [];
unsigned int	data_len;		/* length of data_array */
unsigned int	data_x_offset, data_y_offset; /* where to start writing */
char		pix_array [];		/* bitmap to be written into */
unsigned int	pix_width, pix_height;	/* size of bitmap */
unsigned int 	*bytes_per_line_return;
/* 
add_float_data_to_bitmap --- Linas Vepstas --- 8 April 1989

This routine takes the array of floating point numbers in data_array
(which should all be normalized -- i.e. run between zero and one)
and stuffs them into a bitmap.  Each float represents one bit.
Each bit is turned on if its floating value is greater than #define
DUTY_CYCLE.  The location written to is given by the two offsets.
Writing is performed in the x (horizontal) direction.

The code has been optimized for human readability, not for fast
execution.  
*/
{
int		i,j;
unsigned int 	byte;
char		off;
unsigned int	bit_off, bytes_per_line;

/* scan lines must start on word boundaries */
bytes_per_line = pix_width/32;
if (pix_width%32 != 0) bytes_per_line = bytes_per_line +1;
bytes_per_line = 4 * bytes_per_line;
*bytes_per_line_return = bytes_per_line;

/* loop over x directions */
bit_off = bytes_per_line * data_y_offset;
for (i=0; i<data_len; i++) {
   if (data_array [i] > DUTY_CYCLE)  {
      j = i+data_x_offset;
      byte = bit_off + j/8;	
#ifdef BYTE_SWAP
      off = 0x80 >> (j%8); 
#else
      off = 0x1 << (j%8); 
#endif
      /* set that bit ! */
      pix_array [byte] = pix_array [byte] | off;
      }
   } 
}
