
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
   nx = px;
   ny = py;


   eye[0] = 0.0;
   eye[1] = 0.0;
   eye[2] = -5.0;

   sr = new (SinaiRay *) [nx];
   for (int i = 0; i<nx; i++) sr[i] = new SinaiRay[ny];



}

main ()
{
    SinaiView v (400,400);


    

}


