
// ray.C
//
// Billards ray tracer
//
// put sphere in center of cube.
// cube is 2 units on a side.
// sphere has radius < 1.0



class Ray
{
   public:
      double position[3];   
      double direction[3];
};

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

   private:
      SinaiRay **sr;

};

SinaiView::SinaiView (int px, int py)
{
   nx = px;
   ny = py;

   sr = new (SinaiRay *) [nx];
   for (int i = 0; i<nx; i++) sr[i] = new SinaiRay[ny];
   

}

main ()
{
    SinaiView v (400,400);


}


