
/*
 * sim.C
 *
 * simulate classical magnetic dipole in stern gerlach magnet
 *
 */

#include <math.h>
#include <stdio.h>

main () 
{
   double pos_x, pos_z;  // position
   double vel_x, vel_z;  // veolcity
   double acc_x, acc_z;  // acceleration
   double theta, d_theta;  // angle to z axis
   double phi, d_phi;    // angle to x axis
   double sin_theta, cos_theta;
   double sin_phi, cos_phi;
   int nprec = 0;        // number of precessions
   long long last_i = 0;
   double acc_theta = 0.0;

   double gb0 = 1.0e10;  // gyromagnetic * magnetic field 
   double gb1 = 1.0e9;   // gyromagnetic * magnetic field gradient
   double mub1 = 1.0e4;  // mag moment *mag gradient / mass


   double delta_t = 1.0e-13;  // time step

   long long imax = 48023;  // max iterations
   imax *= 1024*1024;

   theta = 0.6;
   nprec = 0;
   phi = 0.0;
   pos_x = 0.0;
   pos_z = 0.0;
   vel_x = 0.0;
   vel_z = 0.0;
   sin_theta = sin(theta);
   cos_theta = cos(theta);
   sin_phi = sin(phi);
   cos_phi = cos(phi);

   printf ("# \n");
   printf ("# initial theta = %g\n", theta);
   printf ("# delta_t = %g\n", delta_t);
   printf ("# \n");
   printf ("# nprec	i	delta_i	theta	acc_theta	x	z\n");

   for (long long i=0; i<imax; i++)
   {
     double tmp;

     vel_x += mub1 * sin_theta * cos_phi * delta_t;
     vel_z -= mub1 * cos_theta * delta_t;

     pos_x += vel_x * delta_t;
     pos_z += vel_z * delta_t;

     d_theta = gb1 * pos_x * sin_phi * delta_t;
     d_phi = -gb0;
     d_phi += gb1 * (pos_z + pos_x * cos_theta * cos_phi / sin_theta);
     d_phi *= delta_t;

     // theta accumulated only once per precession, below
     // theta += d_theta;
     acc_theta += d_theta;
     phi += d_phi;

#if 0
     if (2.0*M_PI < phi)
     {
        phi -= 2.0 * M_PI;
        sin_phi = sin(phi);
        cos_phi = cos(phi);
        nprec ++;
        printf ("%d	%d	%d	%g	%g	%g	%g\n",
           nprec, i, i-last_i, theta, acc_theta, pos_x, pos_z);
        last_i = i;
        theta += acc_theta;
        acc_theta = 0.0;
     }
     else
#endif

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
           printf ("%d	%d	%d	%20.16g	%g	%g	%g\n",
              nprec, ini, idi, theta, acc_theta, pos_x, pos_z);
           last_i = i;
        }
        theta += acc_theta;

        tmp = cos_theta;
        cos_theta -= sin_theta * acc_theta;
        sin_theta += tmp * acc_theta;

        acc_theta = 0.0;

        if (0.001 > sin_theta && -0.001 < sin_theta)
        {   
           printf ("Error: flirt with divide by zero\n"
                   "theta=%g sin_theta=%g cos_theta=%g\n",
                    theta, sin_theta, cos_theta);
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
