/*-------------------------------------------------------------------*/

void filled_circle (glob, sizex, sizey)
float  		glob [];
unsigned int	sizex, sizey;
{
unsigned int	i,j, globlen;
int		icen, jcen, ip, jp;
float		globby;

globlen = sizex*sizey;
for (i=0; i<globlen; i++) glob [i] = 0.0;

icen = 0.65* sizex;
jcen =  0.8 * sizey;
for (i=0; i<sizey; i++)
   {
   for (j=0; j<sizex; j++) 
      {
      ip = i - icen; 
      jp = j - jcen;
      globby = (float) ((ip*ip+jp*jp) % 2500);
      globby = globby / 2500.0;
      glob [i*sizex +j] = globby; 
      }
   }
}

/*-------------------------------------------------------------------*/

void mandelbrot_out (glob, sizex, sizey,
               re_center, im_center, width, itermax)
float  		glob [];
unsigned int	sizex, sizey;
double		re_center, im_center, width;
int		itermax;
/* this routine fills in the exterior of the mandelbrot set using */
/* the classic algorithm */
{
unsigned int	i,j, globlen;
int		ip, jp;
double		re_start, im_start, delta;
double		re_position, im_position;
double		re, im, tmp;
int		loop;

delta = width / (double) sizex;
re_start = re_center - width / 2.0;
im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);

globlen = sizex*sizey;
for (i=0; i<globlen; i++) glob [i] = 0.0;

im_position = im_start;
for (i=0; i<sizey; i++) {
   if (i%10==0) printf(" start row %d\n", i);
   re_position = re_start;
   for (j=0; j<sizex; j++) {
      re = re_position;
      im = im_position;
      for (loop=1; loop <itermax; loop++) {
         tmp = re*re - im*im + re_position;
         im = 2.0*re*im + im_position;
         re = tmp;
         if ((re*re + im*im) > 4.0) break;
      }    
      glob [i*sizex +j] = ((float) (loop%10)) / 10.0; 
      re_position += delta;
      }
   im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/*-------------------------------------------------------------------*/

void main (argc, argv) 
int argc;
char **argv; 
{
static float	data [1201000];		/* my data array */
unsigned int	data_width, data_height;/* data array dimensions */
float		zmax;
static char	my_pix_data [160000];	/* my pixmap */
unsigned int 	pix_width, pix_height; 	/* pixmap size */
double		re_center, im_center, width;
int		itermax;

printf ("whoops!\n");
data_width = 264;
data_height = 250;
(void) filled_circle (data, data_width, data_height); 

#ifdef MANDEL
re_center = -1.483;
im_center = 0.004;
width = 0.0125;
itermax = 400;
(void) mandelbrot_out (data, data_width, data_height,
               re_center, im_center, width, itermax);
#endif

printf ("filled in the circle \n");

pix_width = data_width;
pix_height = data_height;
(void) add_float_data_to_pixmap (data, data_width*data_height, 0, 0, 
			my_pix_data, pix_width, pix_height)

printf ("created the pixels \n");

(void) blit_pixels_to_new_window 
			(argc, argv, my_pix_data,
			pix_width, pix_height);
}
