
// ray.C
//
// Billards ray tracer
//
// put sphere in center of cube.
// cube is 2 units on a side.
// sphere has radius < 1.0


#include <math.h>
#include <stdio.h>

#include "vvector.h"

class Ray
{
   public:
      void Set (double x[3], double t[3]);
   public:
      double position[3];   
      double direction[3];
};

void Ray::Set (double x[3], double t[3])
{
   VEC_COPY(position, x);
   VEC_COPY(direction,  t);
   VEC_NORMALIZE(direction);
}

class SinaiRay
  : public Ray
{
   public:
      int bounce;
};


class SinaiView
{
   public:
      SinaiView (int, int);
      void ToPixels (void);
      void WriteMTV (const char *filename);
   public:
      int nx;  // dimension of pixel grid
      int ny;

      double eye[3]; // position of eye

   private:
      SinaiRay **sr;
      unsigned int *abgr;
      unsigned int *pack;

};

SinaiView::SinaiView (int px, int py)
{
   int i, j;

   nx = px;
   ny = py;

   eye[0] = 0.0;
   eye[1] = 0.0;
   eye[2] = -5.0;

   sr = new (SinaiRay *) [nx];
   for (i = 0; i<nx; i++) sr[i] = new SinaiRay[ny];

   // initialize ray field
   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         double pixel[3];
         double dir[3];
         pixel[0] = 2.0 * (((double) i) + 0.5)/ ((double) nx) - 1.0;
         pixel[1] = 2.0 * (((double) j) + 0.5)/ ((double) ny) - 1.0;
         pixel[2] = -1.0;
         VEC_DIFF (dir, eye, pixel);
         sr[i][j].Set (pixel, dir);
      }
   }

   abgr = new unsigned int [nx*ny];
   pack = new unsigned int [3*nx*ny/4+1];
   for (i=0; i<nx*ny; i++) abgr[i] = 0;
   for (i=0; i<3*nx*ny/4; i++) pack[i] = 0;
}

/* ==================================== */

void 
SinaiView::ToPixels (void)
{
   int i,j;

   for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
         abgr[nx*j+i] = i%255;
      }
   }

   unsigned int *p = pack;
   for (i=0; i<nx*ny; i++) {
      *p = abgr[i] & 0xffffff;
      ((char *) p) += 3;
   }
}

void 
SinaiView::WriteMTV (const char * filename)
{
   FILE *fh = fopen (filename, "w");
   fprintf (fh, "%d %d\n", nx, ny);
   fwrite (pack, 4, 3*nx*ny/4+1, fh);
   fclose (fh);
}

main ()
{
   SinaiView v (400,400);


   v.ToPixels();
   v.WriteMTV ("junk.mtv");
    

}


