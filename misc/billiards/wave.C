//
// FILE:
// wave.C
//
// FUNCTION:
// Sinai's Billards ray tracer.  Sinai's Billiards is a well-known
// example of a chaotic system.
//
// This has a twist: we do a feynmann path integral over rays.
//
// HISTORY:
// Linas Vepstas
// November 2001


#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ray.h"
#include "vvector.h"

/* ==================================== */
// The PathIntegral class generates the pretty pictures, converts 
// colors, etc.
//
// The Trace() method performs ray-tracing using the classic Sinai
// reflective boundary conditions.
//
// The TraceToroid() method performs ray-tracing using periodic
// toroidial boundary conditions.

class PathIntegral
   : public SinaiBox, 
     public Pixels
{
   public:
      PathIntegral (int, int);
      void Init (void);
      void Trace (void);
      void TraceToroid (void);
      void SumRays (double);
      void AccumIntensity (void);
      void ToPixels (void);

   public:
      int nx;  // dimension of pixel grid
      int ny;

      double omega;  // energy
      double oversample;   // samples per pixel

      double shooter[3]; // direction of source
      double phi;        // direction angle

      complex<double> *side_amplitude;
      complex<double> *front_amplitude;
      complex<double> *back_amplitude;
      int *side_count;
      int *front_count;
      int *back_count;

      double *side_intensity;
      double *front_intensity;
      double *back_intensity;

      double **side_raylens;
      double **front_raylens;
      double **back_raylens;

      int nintense;
};

/* ==================================== */

PathIntegral::PathIntegral (int px, int py)
   : Pixels (px, py)
{
   int i;

   oversample = 1.0;
   omega = 1.0;

   nx = px;
   ny = py;

   shooter[0] = 0.0;
   shooter[1] = 0.0;
   shooter[2] = -1.0;

   side_amplitude = new complex<double> [nx*ny];
   front_amplitude = new complex<double> [nx*ny];
   back_amplitude = new complex<double> [nx*ny];

   side_count = new int [nx*ny];
   front_count = new int [nx*ny];
   back_count = new int [nx*ny];

   side_intensity = new double [nx*ny];
   front_intensity = new double [nx*ny];
   back_intensity = new double [nx*ny];

   side_raylens = new double* [nx*ny];
   front_raylens = new double* [nx*ny];
   back_raylens = new double* [nx*ny];

   Init ();
   for (i=0; i<nx*ny; i++)
   {
      side_intensity[i] = 0.0;
      front_intensity[i] = 0.0;
      back_intensity[i] = 0.0;

      side_raylens[i] = (double *) malloc (8*sizeof (double));
      front_raylens[i] = (double *) malloc (8*sizeof (double));
      back_raylens[i] = (double *) malloc (8*sizeof (double));
   }

   nintense = 0;
}

/* ==================================== */

void
PathIntegral::Init (void)
{
   int i;

   for (i=0; i<nx*ny; i++)
   {
      side_amplitude[i] = 0.0;
      front_amplitude[i] = 0.0;
      back_amplitude[i] = 0.0;
      side_count[i] = 0;
      front_count[i] = 0;
      back_count[i] = 0;
   }
}

/* ==================================== */

complex<double> myexp (double const &ph)
{
    complex<double> retval (cos(ph), sin (ph));
    return retval;
}

/* ==================================== */
void
PathIntegral::TraceToroid(void)
{
   SinaiRay sr;
   int i,j;

   int ox = (int) (oversample * ((double) nx));
   int oy = (int) (oversample * ((double) ny));

   double co = cos (phi);
   double si = sin (phi);

   // rotate by phi radians around y axis
   double tmp =  co*shooter[0] + si*shooter[2];
   shooter[2] = -si*shooter[0] + co*shooter[2];
   shooter[0] = tmp;

   // project rays from source
   for (i=0; i<ox; i++) 
   {
      for (j=0; j<oy; j++) 
      {

         // define initial ray direction and position.
         // we will use an orthonormal projection: the
         // shooter is at infinity, and all rays are 
         // parallel and in phase.  'In phase' means the rays
         // start on a plane that is perpendicular to the
         // ray direction.   This means that the ray trace will
         // be that of a coherent wavefront.
         double pixel[3];
         pixel[0] = 2.0 * (((double) i) + 0.51333)/ ((double) ox) - 1.0;
         pixel[1] = 2.0 * (((double) j) + 0.51666)/ ((double) oy) - 1.0;
         pixel[2] = 13.0;

         // rotate by phi radians around y axis
         tmp      =  co*pixel[0] + si*pixel[2];
         pixel[2] = -si*pixel[0] + co*pixel[2];
         pixel[0] = tmp;

         sr.Init();
         sr.Set (pixel, shooter);
         sr.last_wall = 4;

         // ray trace
         SinaiBox::TraceToroid (sr);

         // map final ray direction to pixel
         // do this by projection
         // treat the four side walls as identical,
         // ignore the front and back walls for now
         int use_ray = 0;
         if ((sr.direction[0] > 0.0)  && 
             (fabs (sr.direction[1]) < sr.direction[0]) &&
             (fabs (sr.direction[2]) < sr.direction[0]))
         {
            use_ray = 1;
         }
#define ROTATED_FACES
#ifdef ROTATED_FACES
         else // 90 degrees
         if ((sr.direction[1] > 0.0)  && 
             (fabs (sr.direction[0]) < sr.direction[1]) &&
             (fabs (sr.direction[2]) < sr.direction[1]))
         {
            double tmp = sr.direction[0];
            sr.direction[0] = sr.direction[1];
            sr.direction[1] = -tmp;
            use_ray = 1;
         }
         else   // 180 degrees
         if ((sr.direction[0] < 0.0)  && 
             (fabs (sr.direction[1]) < -sr.direction[0]) &&
             (fabs (sr.direction[2]) < -sr.direction[0]))
         {
            sr.direction[0] = - sr.direction[0];
            sr.direction[1] = - sr.direction[1];
            use_ray = 1;
         }
         else // 270 degrees
         if ((sr.direction[1] < 0.0)  && 
             (fabs (sr.direction[0]) < -sr.direction[1]) &&
             (fabs (sr.direction[2]) < -sr.direction[1]))
         {
            double tmp = sr.direction[1];
            sr.direction[1] = sr.direction[0];
            sr.direction[0] = -tmp;
            use_ray = 1;
         }
#endif

         if (use_ray)
         {

             // first convert ray direction to grid coords
             double x = sr.direction[1] / sr.direction[0];
             double y = sr.direction[2] / sr.direction[0];

             int px = (int) (((double) nx) * 0.5 * (x+1.0));
             int py = (int) (((double) ny) * 0.5 * (y+1.0));

             if ((0 > px) || (px >=nx) || (0 > py) || (py >=ny))
             {
                printf ("duude out of bounds !! %d %d \n", px, py);
             }

#if DO_PHASE
             // next perform phase summation
             double phase = omega * sr.distance;
             side_amplitude [nx*py+px] += myexp (phase);
             side_count [nx*py+px] ++;
#endif

             // make sure we have enough room to store the ray results
             int n = side_count [nx*py+px];
             int bits = 1;
             while (n >>= 1) bits++;
             if (3<bits)
             {
                size_t sz = 1<< bits;
                side_raylens[nx*py+px] = 
                     (double *) realloc (side_raylens[nx*py+px],
                          sz*sizeof (double));
             }

             // store the ray length
             n = side_count [nx*py+px];
             side_raylens[nx*py+px][n] = sr.distance;
             side_count [nx*py+px] ++;
          
#if 0
printf ("duude ph= %f ", phase);
printf ("ampd= %f %f ", real(amplitude [nx*py+px]), imag(amplitude [nx*py+px]));
printf ("duude amp= %f\n", abs<double>(amplitude [nx*py+px]));
#endif
         }
      }

      if (1.0 < oversample)
      {
         if (0 == i % (int)oversample) { printf ("."); fflush (stdout); }
      } else
      {
         if (0 == i) { printf ("."); fflush (stdout); }
      }
   }

}

/* ==================================== */

void 
PathIntegral::AccumIntensity (void)
{

   for (int i=0; i<nx*ny; i++)
   {
      if (0 < side_count[i])
      {
         double re = real (side_amplitude[i]) / ((double) side_count[i]);
         double im = imag (side_amplitude[i]) / ((double) side_count[i]);
         side_intensity[i] += re*re+im*im;

         re = real (front_amplitude[i]) / ((double) front_count[i]);
         im = imag (front_amplitude[i]) / ((double) front_count[i]);
         front_intensity[i] += re*re+im*im;

         re = real (back_amplitude[i]) / ((double) back_count[i]);
         im = imag (back_amplitude[i]) / ((double) back_count[i]);
         front_intensity[i] += re*re+im*im;
      }

   }

   nintense ++;

   // zero out amplitudes
   Init ();
}

/* ==================================== */

void 
PathIntegral::SumRays (double k)
{
   nintense = 0;

   for (int i=0; i<nx*ny; i++)
   {
      for (int nr= side_count[i]; nr>0; nr++)
      {
         double phase = k * side_raylens[i][nr];
         side_amplitude [i] += myexp (phase);
      }
   }
}

/* ==================================== */

void 
PathIntegral::ToPixels (void)
{

   for (int i=0; i<nx*ny; i++)
   {
      double red = 0.0;
      double green = 0.0;
      double blue = 0.0;

      if (0 < side_count[i])
      {
         // red = 255.0 * abs<double> (side_amplitude[i]) / ((double) side_count[i]);
         red = 255.0 * side_intensity[i] / ((double) nintense);
      }

      abgr[i] = 0xff & ((unsigned int) red);
      abgr[i] |= (0xff & ((unsigned int) green)) << 8;
      abgr[i] |= (0xff & ((unsigned int) blue)) << 16;
   }
}

/* ==================================== */

main (int argc, char * argv[])
{
   PathIntegral v (400,400);

   if (6 > argc) {
      printf ("Usage: %s <fileout> <radius> <omega> <samples> <maxdist> [<niterations> [<max manhattan>]]\n", argv[0]);
      exit (1);
   }

   char * outfile = argv[1];
   double radius = atof (argv[2]);
   double omega = atof (argv[3]);
   double samples = atof (argv[4]);
   double maxdist= atof (argv[5]);

   int niter = 1000000;
   if (7 == argc) niter = atoi (argv[6]);

   int manhat = 1000000;
   if (8 == argc) manhat = atoi (argv[7]);

   v.radius = radius;
   v.omega = omega;
   v.oversample = samples;
   v.max_distance = maxdist;
   v.niterations = niter;
   v.max_manhattan = manhat;

   // v.Trace();
   v.TraceToroid();

   int nframes = 240;
   double delta_k = pow (10.0, 6.0 / ((double) nframes));
   double k = 0.001;
   for (int i=0; i<nframes; i++)
   {
      char buff [200];
      v.SumRays (k);
      v.AccumIntensity ();
      v.ToPixels ();
      sprintf (buff, "%s-%d.mtv", outfile, i);
      v.WriteMTV (buff);
      k *= delta_k;
   }
}

/* ===================== end of file ====================== */
