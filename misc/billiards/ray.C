
// ray.C
//
// Billards ray tracer
//
// put sphere in center of cube.
// cube is 2 units on a side.
// sphere has radius < 1.0


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
   public:
      int nx;  // dimension of pixel grid
      int ny;

      double eye[3]; // position of eye

   private:
      SinaiRay **sr;

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



}

main ()
{
    SinaiView v (400,400);


    

}


