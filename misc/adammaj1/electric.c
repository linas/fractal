
/* 


   http://linas.org/art-gallery/escape/ray.html
   image : 
   http://linas.org/art-gallery/escape/electric.gif
   
   
   "Gradient of the potential

The gradient is more interesting. The gradient 2Df always points 'uphill':  

2Df = f zn Dzn / |zn|2 log |zn|

The picture below shows the modulus of the gradient.  Of course, it resembles a coronal discharge.  
If one were to cut out a metal plate in the shape of an M-set, apply a high voltage to it, 
it is exactly this that you would see, because it is the areas of high electric field (gradient) where the coronal ionization will occur. 
Maybe make a chic lightning rod for the barn. " Linas Vepstas
   
  
   c console program:
   
   
   --------------------------------
   1. draws Mandelbrot set for Fc(z)=z*z +c
   using Mandelbrot algorithm ( boolean escape time )
   -------------------------------         
   2. technique of creating ppm file is  based on the code of Claudio Rocchini
   http://en.wikipedia.org/wiki/Image:Color_complex_plot.jpg
   create 24 bit color graphic file ,  portable pixmap file = PPM 
   see http://en.wikipedia.org/wiki/Portable_pixmap
   to see the file use external application ( graphic viewer)
   
   
   complex point c -> virtual 2D array -> memory 1D array -> ppm file on the disc -> png file 
   
   C -> pixel (iX,iY)  -> index k  -> 24bit color 
   
   -----
   https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
   complex numbers are built in type 
 
 --------------
 formated with emacs
   -------------
   to compile : 

 
 
   gcc electric.c -lm -Wall 
 
 
   ./a.out
   
   
   to convert to png using ImageMagic

   convert electric.ppm electric.png  



   ----------------------
   git add  electric.png electric.c
   git commit -m " electric"
   push -u origin master



 
 
*/
#include <stdio.h>
#include <stdlib.h>		// malloc
#include <math.h>
#include <complex.h> // https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
 
 

 

 
 
/* screen ( integer) coordinate */

const int iWidth  = 1000; 
const int iHeight = 1000;


/* world ( double) coordinate = parameter plane*/
// double complex C =  Cx + Cy*I ;
const double CxMin=-2.3;
const double CxMax= 0.7;
const double CyMin=-1.5;
const double CyMax= 1.5;

/* */
double PixelWidth; //=(CxMax-CxMin)/iWidth;
double PixelHeight; // =(CyMax-CyMin)/iHeight;


/* color component ( R or G or B) is coded from 0 to 255 */
/* it is 24 bit color RGB file */
int ColorBytes = 3; // 3*8 = 24 bit color        

/* iterations  */
const int IterationMax=18;

/* bail-out value , radius of circle ;  */
const double EscapeRadius=3.0;
double log_2; // = log(2.0);
       
// memmory virtual 1D array 
unsigned char *data;       
size_t MemmorySize;       
    
    
void GiveLinasColor(double position , int k, unsigned char c[])
{
  /* based on the code by Linas Vepstas January 16 1994 : void make_cmap (void) */

   
  int i;
  int iMax = 240;
  i=(int)(iMax-1)*position;  
  c[0] = c[1] = c[2] = 0;
  /* set up a default look up table */
  /* ramp up to blue */
  if (i<60) {
    c[k] = 0;
    c[k+1] = 0;
    c[k+2] = (unsigned char) i*3;
  }
  /* ramp down from blue, up to green */
  if (i>=60 && i<120) {
    c[k] = 0;
    c[k+1] = (unsigned char) (i-60)*3;
    c[k+2] = (unsigned char) (120-i)*3;
  }
  /* ramp from green to yellow */
  if (i>=120 && i<180) {
    /* vlt[i].r = (char) (((i-120)*7) / 2); */
    c[k] = (unsigned char) (210 - (7*(180-i)*(180-i)) / 120);
    c[k+1] = (unsigned char) (210 -i/4);
    c[k+2] = 0;
  }
  /* ramp from yellow to red (pink) */
  if (i>=180 && i<iMax) {
    c[k] = (unsigned char) (210 + (3*(i-180))/4);
    c[k+1] = (unsigned char) (510 - 2*i);
    c[k+2] = (unsigned char) (i-180)/3;
  }
   
}

       
        
/* 
   gives position ( index) in 1D virtual array  of 2D point (iX,iY) from ; uses also global variable iWidth 
   without bounds check !!
*/
int f(int ix, int iy)
{ return ColorBytes*(ix + iy*iWidth); }
        
        
        
double complex give_c(int iX, int iY){
  double Cx,Cy;
  Cy=CyMin + iY*PixelHeight;
  
  Cx=CxMin + iX*PixelWidth;
   
  return Cx+Cy*I;
 
 
}
 

double complex GiveGradient(double complex C ){

  double g; 
  int i; // iteration
  double complex Z  = 0.0; // initial value for iteration Z0
  double complex dZ = 1.0; // derivative with respect to c
 // https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/def_cqp#Derivative_with_respect_to_c 
  
  
  // iteration
  for(i=0;i<IterationMax;i++)
    {
      
      dZ = 2 * Z * dZ +1.0; 
      Z=Z*Z+C; 
      if(cabs(Z)>EscapeRadius) break; // bailout test 
    }
   
     
  
    
  // the renormalized, fractional iteration count
  // m(R) = n+1 - log(log |zn|) / log 2
  // http://linas.org/art-gallery/escape/escape.html
  double cabsZ = cabs(Z);  
  
  double m = i + 1.0 - log(log(cabsZ))/log_2;
  //m = m/IterationMax; // normalize = map to [0,1]
  
  
  // the Douady-Hubbard potential is just f = e-m log2 = 2-m
  double potential = pow(2,-m);
  // 2Df = m zn Dzn / |zn|2 log |zn|
 
  g = potential*Z*dZ/(cabsZ*cabsZ*log(cabsZ));
  
  return g;
}
 


int ComputeAndSavePixelColor(int iX, int iY){
 
  
  complex double C;
  double complex g;
  double gm;
  int k; // index of the 1D array
   
  C = give_c(iX, iY);
   
  g = GiveGradient(C);
  gm = cabs(g); 
  
  // index of 1D memory array
  k = f(iX, iY);  
  
   
  //Apply this method to both exterior and interior 
  GiveLinasColor(gm , k,  data); // 
      
    
 
   
  return 0;
}
 
 
 
int setup(){
 
  //
  PixelWidth=(CxMax-CxMin)/iWidth;
  PixelHeight=(CyMax-CyMin)/iHeight;
  //
  log_2 = log(2.0);  
  //
  MemmorySize = iWidth * iHeight * ColorBytes * sizeof (unsigned char);	// https://stackoverflow.com/questions/492384/how-to-find-the-sizeof-a-pointer-pointing-to-an-array
        
  /* create dynamic 1D arrays for colors   */
  data = malloc (MemmorySize);
  if (data == NULL )
    { fprintf (stderr, " Error: Could not allocate memory !!! \n"); return 1; }

  printf (" No errors. End of setup \n");
  return 0;

}
 
 
 
 
 
// save dynamic "A" array to pgm file 
int SaveArray_2_PPM_file (unsigned char A[])
{

  FILE *fp;
  const unsigned int MaxColorComponentValue = 255;	/* color component is coded from 0 to 255 ;  it is 8 bit color file */
  char *filename = "electric.ppm";
  char *comment = "# ";		/* comment should start with # */

  /* save image to the pgm file  */
  fp = fopen (filename, "wb");	/*create new file,give it a name and open it in binary mode  */
  fprintf (fp, "P6\n %s\n %u %u\n %u\n", comment, iWidth, iHeight, MaxColorComponentValue);	/*write header to the file */
  fwrite (A, MemmorySize, 1, fp);	/*write A array to the file in one step */
  printf ("File %s saved. \n", filename);
  fclose (fp);

  return 0;
}


 
 
void CreateImage(){
  int iX,iY; // screen = integer coordinate in pixels       

  // fill the array = render image = scanline 2D  of virtual 2D array 
  for(iY=0;iY<iHeight;iY++)
    for(iX=0;iX<iWidth;iX++)
      ComputeAndSavePixelColor(iX, iY); 
      	
      	
  SaveArray_2_PPM_file (data);     	  
} 
 
 
 
void info(){

  printf(" Parameter plane ( c plane) with Mandelbrot set for complex quadratic polynomial fc(z) = z^2 + c\n ");
  printf(" Rectangle part of 2D parameter plane: corners: \n CxMin = %f;   CxMax = %f;  CyMin = %f; CyMax = %f \n ", CxMin, CxMax, CyMin, CyMax);
  printf(" center and radius: \n CenterX = %f;   CenterY = %f;  radius = %f\n ", (CxMax+CxMin)/2.0, (CyMax+CyMin)/2.0, fabs(CyMax-CyMin)/2.0);
  printf(" Mag = zoom = %f\n ",  2.0/fabs(CyMax-CyMin));
  printf("PixelWidth = %f and PixelHeight =%f\n", PixelWidth, PixelHeight);
  printf(" Escape Radius = %f\n ", EscapeRadius);
  printf(" Iteration Max = %d\n ", IterationMax);



} 
 
 
 
void close(){
 
  info(); 
  free(data); 
}
 
 
 
int main()
{
 
  setup();      
  CreateImage();     
  close();
  
        
  return 0;
}
