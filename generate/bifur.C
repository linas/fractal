/*-------------------------------------------------------------------*/
/* This routine does bifurcation diagrams */

void 
MakeBifurWrap (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm
	MakeBifurCB cb)
{
   int		i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   
   delta = width / (double) sizex;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;

   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      if (i%10==0) printf(" start row %d\n", i);

		cb (&glob[i*sizex], sizex, re_center, width, im_position, itermax, renorm);
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
