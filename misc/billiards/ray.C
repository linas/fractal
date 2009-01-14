//
// FILE:
// ray.C
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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ray.h"
#include "vvector.h"
#include "intersect.h"

/* ==================================== */

void 
Ray::Set (double x[3], double t[3])
{
   VEC_COPY(position, x);
   VEC_COPY(direction,  t);
   VEC_NORMALIZE(direction);
}

void
Ray::Set (double x0, double x1, double x2, 
          double t0, double t1, double t2)
{
   position[0] = x0;
   position[1] = x1;
   position[2] = x2;
   direction[0] = t0;
   direction[1] = t1;
   direction[2] = t2;
   VEC_NORMALIZE(direction);
}

/* ==================================== */

double 
Ray::Intersect (Ray& plane)
{
   double p2[3];
   VEC_SUM (p2, position, direction);

   int valid;
   double result[3];
   INTERSECT (valid, result, 
             plane.position, plane.direction,
             position, p2);

   if (FALSE == valid) return -1.0/DEGENERATE_TOLERANCE;
   
   // The 'direction' should be parallel to the result, and thus,
   // the dot product should be the length.
   
   double dot;
   VEC_DIFF (p2, result, position);
   VEC_DOT_PRODUCT (dot, p2, direction);

   return dot;
}

/* ==================================== */

void
Ray::Bounce (Ray& plane)
{
   double p2[3];
   VEC_SUM (p2, position, direction);

   int valid;
   double result[3];
   INTERSECT (valid, result, 
             plane.position, plane.direction,
             position, p2);

   if (FALSE == valid) return;
   VEC_COPY (position, result);
   VEC_REFLECT (direction, direction, plane.direction);  
}

/* ==================================== */

void
Ray::Forward (Ray& plane)
{
   double p2[3];
   VEC_SUM (p2, position, direction);

   int valid;
   double result[3];
   INTERSECT (valid, result, 
             plane.position, plane.direction,
             position, p2);

   if (FALSE == valid) return;

   VEC_COPY (position, result);
}

/* ==================================== */

double
Ray::Sphere (double radius)
{
   double b[3], pact, r2, s;
   double radial[3], dist, oldpos[3];

   VEC_PERP (b, position, direction);  
   VEC_LENGTH (pact, b);  // impact paramter

   r2 = radius*radius;
   if (r2 < pact) return -1.0;
   s = - sqrt (r2 - pact);

   VEC_COPY (oldpos, position);
   VEC_BLEND (radial, 1.0, b, s, direction);
   VEC_COPY (position, radial);
   VEC_NORMALIZE (radial);
   VEC_REFLECT (direction, direction, radial);

   VEC_DIFF (oldpos, position, oldpos);
   VEC_LENGTH (dist, oldpos);

   return dist;
}

/* ==================================== */

void
SinaiRay::Init (void)
{
   sphere_hits = 0;
   last_wall = -1;
   distance = 0.0;
   for (int i=0; i<6; i++) bounces[i] = 0;
}


SinaiRay::SinaiRay (void)
{
   Init();
}

/* ==================================== */

SinaiBox::SinaiBox (void)
{
   radius = 0.5;
   niterations = 10000;
   max_manhattan = 10000;
   max_distance = 10000.0;
	max_sphere_hits = 10000;

   walls[0].Set(-1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
   walls[1].Set(1.0, 0.0, 0.0, -1.0, 0.0, 0.0);
   walls[2].Set(0.0, 1.0, 0.0, 0.0, -1.0, 0.0);
   walls[3].Set(0.0, -1.0, 0.0, 0.0, 1.0, 0.0);
   walls[4].Set(0.0, 0.0, 1.0, 0.0, 0.0, -1.0);
   walls[5].Set(0.0, 0.0, -1.0, 0.0, 0.0, 1.0);
}

/* ==================================== */

#define NNN (-1.0 + 100.0*DEGENERATE_TOLERANCE)
#define PPP (1.0 - 100.0*DEGENERATE_TOLERANCE)

void 
SinaiBox::Trace (SinaiRay &sr)
{
   for (int n=0; n<niterations; n++)
   {
      double dist = 0.0;

      // first, see if ray bounces off sphere
      dist = sr.Sphere (radius);
      if (0.0 < dist) 
      {
         sr.sphere_hits ++;
         sr.distance += dist;
         sr.last_wall = -1;
      }

      // next, trace ray to wall.
      int next_wall = -1;
      double nearest = 1000000.0;
      for (int iwall=0; iwall<6; iwall ++)
      {
         if (iwall == sr.last_wall) continue;

         double dist = sr.Intersect (walls[iwall]);
         if (0.0 > dist) continue;
         if (nearest > dist) { nearest = dist; next_wall = iwall; }
      }
      sr.bounces[next_wall] ++;
      sr.last_wall = next_wall;
      sr.distance += nearest;

      sr.Bounce (walls[next_wall]);

#if BROKEN
      // check for corner conditions
      if (((NNN >= sr.position[0]) ||
           (PPP <= sr.position[0])) &&
           (0 != sr.last_wall) &&
           (1 != sr.last_wall))
      {
         sr.direction[0] = -  sr[i].direction[0];
      }

      if (((NNN >= sr.position[1]) ||
           (PPP <= sr.position[1])) &&
           (2 != sr.last_wall) &&
           (3 != sr.last_wall))
      {
         sr.direction[1] = -  sr[i].direction[1];
      }

      if (((NNN >= sr.position[2]) ||
           (PPP <= sr.position[2])) &&
           (4 != sr.last_wall) &&
           (5 != sr.last_wall))
      {
         sr.direction[2] = -  sr[i].direction[2];
      }
#endif

      if (sr.distance > max_distance) break;
      if (sr.bounces[next_wall] >= max_manhattan) break;
      if (sr.sphere_hits >= max_sphere_hits) break;
   }
}

/* ==================================== */

void 
SinaiBox::TraceToroid (SinaiRay &sr)
{
   for (int n=0; n<niterations; n++)
   {
      double dist = 0.0;

      // first, see if ray bounces off sphere
      dist = sr.Sphere (radius);
      if (0.0 < dist) 
      {
         sr.sphere_hits ++;
         sr.distance += dist;
         sr.last_wall = -1;
      }

      // next, trace ray to wall.
      int next_wall = -1;
      double nearest = 1000000.0;
      for (int iwall=0; iwall<6; iwall ++)
      {
         if (iwall == sr.last_wall) continue;

         double dist = sr.Intersect (walls[iwall]);
         if (0.0 > dist) continue;
         if (nearest > dist) { nearest = dist; next_wall = iwall; }
      }
      sr.bounces[next_wall] ++;
      sr.last_wall = next_wall;
      sr.distance += nearest;

      sr.Forward (walls[next_wall]);

      // enforce periodic boundary conditions
      sr.position[0] += 2.0 * walls[next_wall].direction[0];
      sr.position[1] += 2.0 * walls[next_wall].direction[1];
      sr.position[2] += 2.0 * walls[next_wall].direction[2];

      // make sure we set the flags correctly
      switch (sr.last_wall) {
         case 0: sr.last_wall = 1; break;
         case 1: sr.last_wall = 0; break;
         case 2: sr.last_wall = 3; break;
         case 3: sr.last_wall = 2; break;
         case 4: sr.last_wall = 5; break;
         case 5: sr.last_wall = 4; break;
         default: break;
      }
         
      if (sr.distance > max_distance) break;
      if (sr.bounces[next_wall] >= max_manhattan) break;
      if (sr.sphere_hits >= max_sphere_hits) break;
   }
}

/* ==================================== */

Pixels::Pixels (int px, int py)
{
   int i;

   nx = px;
   ny = py;

   abgr = new unsigned int [nx*ny];
   pack = new unsigned int [3*nx*ny/4+1];
   for (i=0; i<nx*ny; i++) abgr[i] = 0;
   for (i=0; i<3*nx*ny/4; i++) pack[i] = 0;
}

void 
Pixels::TestPattern (void)
{
   int i,j;

   for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
         abgr[nx*j+i] = i%255 | (j%255 << 8);
      }
   }
}

/* ==================================== */

void 
Pixels::WriteMTV (const char * filename)
{
   // first, pack the pixels byte by byte
   unsigned char *p = (unsigned char *) pack;
   for (int i=0; i<nx*ny; i++) {
      *((unsigned int *) p) = abgr[i] & 0xffffff;
      p += 3;
   }

   // write out MTV format pixmap.
   FILE *fh = fopen (filename, "w");
   fprintf (fh, "%d %d\n", nx, ny);
   fwrite (pack, 4, 3*nx*ny/4+1, fh);
   fclose (fh);
}

/* ===================== end of file ====================== */
