//
// FILE:
// ray.h
//
// FUNCTION:
// Sinai's Billards ray tracer.  Sinai's Billiards is a well-known
// example of a chaotic system.
//
// Design:
// put sphere in center of cube.
// cube is 2 units on a side.
// sphere has radius < 1.0
//
// HISTORY:
// Linas Vepstas
// November 2001

#ifndef SINAI_RAY_H
#define SINAI_RAY_H

/* ==================================== */
// Class Ray is just a plain ray

// Intersect() determines the intersection between the ray and the plane.
// It returns the distance to the intersection if the intersection
// is forward along the ray.  Otherwise it returns minus the distance.
//
// Bounce() determines if a forward trace of the ray will intersect the
// plane.  If it does, then the ray will be traced forward and reflected
// off the plane.  
//
// Forward() determines if a forward trace of the ray will intersect the
// plane.  If it does, then the ray will be traced forward to the plane.
//
// Sphere() returns a positive number if the ray intersects with a sphere
// of radius R centered at the origin; otherwise negative.  The returned
// value is the distance the ray travelled to hit the sphere.

class Ray
{
   public:
      void Set (double x[3], double t[3]);  // assign value
      void Set (double, double, double, double, double, double);
      double Intersect (Ray& plane);
      void Bounce (Ray& plane);
      void Forward (Ray& plane);
      double Sphere (double r);
   public:
      double position[3];   
      double direction[3];
};

/* ==================================== */
// class inherits from ray; collects statistics about
// the path of the ray as its traced, including distance 
// traveled, number of bounces, etc.

class SinaiRay
  : public Ray
{
   public:
      SinaiRay (void);
      double distance;  // distance travelled by ray   
      int last_wall;    // last wall of intersection
      int bounces[6];   // how many bounces have occured off each wall
      int sphere_hits;  // how many times ray bounced off sphere
};

/* ==================================== */
// The class SinaiBox defines the boundry conditions, the max
// iterations, the shape of the actual box.
//
// The Trace() method performs ray-tracing using the classic Sinai
// reflective boundary conditions.
//
// The TraceToroid() method performs ray-tracing using periodic
// toroidial boundary conditions.

class SinaiBox
{
   public:
      SinaiBox (void);
      void Trace (SinaiRay &);
      void TraceToroid (SinaiRay &);
   public:

      int niterations;
      double max_distance;
      int max_manhattan;
      double radius;

   private:
      Ray walls[6]; // left, right, top, bottom, front, back;
};

/* ==================================== */
// The Pixels class handles colors
//

class Pixels
{
   public:
      Pixels (int, int);
      void TestPattern (void);
      void WriteMTV (const char *filename);

   public:
      int nx;  // dimension of pixel grid
      int ny;

      unsigned int *abgr;
   private:
      unsigned int *pack;
};

#endif /* SINAI_RAY_H */
/* ===================== end of file ====================== */
