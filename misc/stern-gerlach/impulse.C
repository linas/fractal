
/*
 * impulse.C
 *
 * Simulate classical magnetic dipole absorbing soft photon impulses.
 *
 * Linas Vepstas July 2001
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main () 
{
   double gb0 = 2.0e-5;  // gyromagnetic * magnetic field 

   long long imax = 48023;  // max iterations
   imax *= 1024*1024;

   int nprec_max = 20000;

   printf ("# \n");
   printf ("# random start \n");
   printf ("# nprec_max = %d\n", nprec_max);
   printf ("# gb0 = %g\n", gb0);
   printf ("# \n");
   printf ("# nprec	i	delta_i	initial_theta	theta\n");
   fflush (stdout);

   double initial_theta;
   for (int loop=1; loop<100; loop++)
   {
      long long last_i = 0;
      double theta, d_theta;  // dipole angle to z axis
      double phi, d_phi;    // dipole angle to x axis
      double sin_theta, cos_theta;
      double sin_phi, cos_phi;
      int nprec = 0;        // number of precessions

      initial_theta = 1.8 * loop;
      initial_theta -= M_PI * ((int) (((double)1.8*loop)/M_PI));
      theta = initial_theta;
      phi = 0.0;
      sin_theta = sin(theta);
      cos_theta = cos(theta);
      sin_phi = sin(phi);
      cos_phi = cos(phi);

      nprec = 0;

      for (long long i=0; i<imax; i++)
      {
        double tmp;
        double magnet_theta;  // photon mag field orientation
        double magnet_phi;    // photon mag field orientation
   
        // theta in upper hemisphere only
        magnet_theta = (M_PI * rand()) / (2.0 * RAND_MAX);
        magnet_phi = (2.0 * M_PI * rand()) / (1.0 * RAND_MAX);
   
        double mag_sin_theta = sin (magnet_theta);
        double mag_cos_theta = cos (magnet_theta);
        double mag_sin_phi = sin (magnet_phi);
        double mag_cos_phi = cos (magnet_phi);
   
   #if 0
        // we don't really need to calculate these
        double impulse_x;
        double impulse_y;
        double impulse_z;
   
        impulse_z = cos_phi * mag_sin_phi - sin_phi * mag_cos_phi;
        impulse_z *= sin_theta * mag_sin_theta;
        impulse_z *= gb0;
   
        impulse_x = mag_cos_theta * sin_theta * sin_phi;
        impulse_x -= cos_theta * mag_sin_theta * mag_sin_phi;
        impulse_x *= gb0;
   
        impulse_y = cos_theta * mag_sin_theta * mag_cos_phi;
        impulse_y -= mag_cos_theta * sin_theta * cos_phi;
        impulse_y *= gb0;
   #endif
   
        // d_theta is defined as 
        // d_theta = -impulse_z / sin_theta;
        // but skip the multiply and divide by sin_theta
        d_theta = cos_phi * mag_sin_phi - sin_phi * mag_cos_phi;
        d_theta *= -gb0 * mag_sin_theta;
   
        theta += d_theta;
   
        // d_phi is defined as
        // d_phi = (impulse_y - d_theta * cos_theta * sin_phi ) / (sin_theta * cos_phi);
        // but we can avoid division problems by using the long-hand
        if (0.0 != sin_theta)
        {
           d_phi = cos_phi * mag_cos_phi + sin_phi * mag_sin_phi;
           d_phi *= cos_theta * mag_sin_theta / sin_theta;
           d_phi -= mag_cos_theta;
           d_phi *= gb0;
   
           phi += d_phi;
        }
        else 
        {
           d_phi = 0.0;
           phi = 0.0;
        }
   
   
        // compute new sine and cosine of phi
        if (2.0*M_PI < phi)
        {
           phi -= 2.0 * M_PI;
           sin_phi = sin(phi);
           cos_phi = cos(phi);
           nprec ++;
           if (0 == nprec%1000) 
           {
              long long di = i - last_i;
              int idi = di;
              long long ni = i/(1024*1024);
              long ini = ni;
              printf ("%d	%d	%lld	%g	%g\n",
                 nprec, ini, di, initial_theta, theta);
              fflush (stdout);
              last_i = i;
           }
        }
        else
        if (-2.0*M_PI > phi)
        {
           phi += 2.0 * M_PI;
           sin_phi = sin(phi);
           cos_phi = cos(phi);
           nprec --;
           if (0 == nprec%1000) 
           {
              long long di = i - last_i;
              int idi = di;
              long long ni = i/(1024*1024);
              long ini = ni;
              printf ("%d	%d	%lld	%g	%g\n",
                 nprec, ini, di, initial_theta, theta);
              fflush (stdout);
              last_i = i;
           }
           if (nprec_max < -nprec) break;
        }
        else 
        {
          tmp = cos_phi;
          cos_phi -= sin_phi * d_phi;
          sin_phi += tmp * d_phi;
        }
   
        // compute new sine and cosine of theta
        if (0.0 > theta)
        {
           theta = -theta;
           if (0.0 < d_phi) { phi -= M_PI; }
           else { phi += M_PI; }
           cos_theta = cos(theta);
           sin_theta = sin(theta);
           sin_phi = sin(phi);
           cos_phi = cos(phi);
        }
        else 
        if (M_PI < theta)
        {
           theta = 2.0 * M_PI - theta;
           if (0.0 < d_phi) { phi -= M_PI; }
           else { phi += M_PI; }
           cos_theta = cos(theta);
           sin_theta = sin(theta);
           sin_phi = sin(phi);
           cos_phi = cos(phi);
        }
        else 
        {
           tmp = cos_theta;
           cos_theta -= sin_theta * d_theta;
           sin_theta += tmp * d_theta;
        }
   
      }
   }
}
