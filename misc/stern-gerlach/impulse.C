
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
   double theta, d_theta;  // dipole angle to z axis
   double phi, d_phi;    // dipole angle to x axis
   double sin_theta, cos_theta;
   double sin_phi, cos_phi;
   int nprec = 0;        // number of precessions
   long long last_i = 0;
   double acc_theta = 0.0;

   double gb0 = 1.0e-3;  // gyromagnetic * magnetic field 

   long long imax = 48023;  // max iterations
   imax *= 1024*1024;

   theta = 0.6;
   nprec = 0;
   phi = 0.0;
   sin_theta = sin(theta);
   cos_theta = cos(theta);
   sin_phi = sin(phi);
   cos_phi = cos(phi);

   printf ("# \n");
   printf ("# initial theta = %g\n", theta);
   printf ("# \n");
   printf ("# nprec	i	delta_i	theta	phi\n");

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

     // compute new sine and cosine of theta
     tmp = cos_theta;
     cos_theta -= sin_theta * d_theta;
     sin_theta += tmp * d_theta;


     // compute new sine and cosine of phi
     if (2.0*M_PI < phi)
     {
        phi -= 2.0 * M_PI;
        sin_phi = sin(phi);
        cos_phi = cos(phi);
        nprec ++;
        // if (0 == nprec%10) 
        {
           long long di = i - last_i;
           int idi = di;
           long long ni = i/(1024*1024);
           long ini = ni;
           printf ("%d	%d	%d	%g	%g\n",
              nprec, ini, idi, theta, phi);
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
        // if (0 == nprec%10) 
        {
           long long di = i - last_i;
           int idi = di;
           long long ni = i/(1024*1024);
           long ini = ni;
           printf ("%d	%d	%d	%g	%g\n",
              nprec, ini, idi, theta, phi);
           last_i = i;
        }
     }
     else 
     {
       tmp = cos_phi;
       cos_phi -= sin_phi * d_phi;
       sin_phi += tmp * d_phi;
     }

   }


}
