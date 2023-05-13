
/*
 * pauli.C
 *
 * Simulate an SU(2) group element absorbing soft photon impulses.
 * The photons are crudely distributed so that thier average creates
 * a magnetic field along z direction.  
 *
 * Linas Vepstas May 2002
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>

#define I (complex<double>(0,1));



class su2
{
	public:
		su2 (void);
		su2 (complex<double>, 
		          complex<double>, 
					 complex<double>, 
					 complex<double>);

		/* initializer. Assumes x,y,z are tangent vector */
		su2 (double x, double y, double z);

		/* multiply */
		su2& operator *= (const su2&);

		/* move along tangent vector; assumes x,y,z are normed group elt */
		su2& delta (double x, double y, double z);

		/*  conjugate with tangent vector; assumes x,y,z are normed group elt */
		su2& conj (double x, double y, double z);
		
		void print (void);
		
	private:
		/* SU2 as point on sphere in R4 */
		double a0, a1, a2, a3;
};


su2::su2 (void)
{
	a1 = a2 = a3 = 0.0;
	a0 = 1.0;
}

/* take any matrix, project onto SU2 */
su2::su2 (complex<double> m11, 
          complex<double> m12, 
			 complex<double> m21, 
			 complex<double> m22)
{
	complex<double> b0, b1, b2, b3;
	b0 = m11 + m22;
	b1 = (m12 + m21) * I;
	b2 = - (m12 - m21);
	b3 = (m11 - m22) * I;
	
	a0 = real(b0);
	a1 = real(b1);
	a2 = real(b2);
	a3 = real(b3);

	double norm = 1.0 / sqrt (a0*a0+a1*a1+a2*a2+a3*a3);
	a0 *= norm;
	a1 *= norm;
	a2 *= norm;
	a3 *= norm;
}

su2 :: su2 (double x, double y, double z)
{
	double theta = sqrt (x*x+y*y+z*z);

	if (0.0 != theta)
	{
		double s = sin (theta);
		a1 = x * s / theta;
		a2 = y * s / theta;
		a3 = z * s / theta;
		a0 = sqrt (1.0 - (a1*a1+a2*a2+a3*a3));
	}
	else
	{
		a1 = a2 = a3 = 0.0;
		a0 = 1.0;
	}
}

su2&
su2 :: operator *= (const su2& b)
{
	complex<double> ma11=0.0, ma12=0.0, ma21=0.0, ma22=0.0;
	complex<double> mb11=0.0, mb12=0.0, mb21=0.0, mb22=0.0;
	complex<double> mc11=0.0, mc12=0.0, mc21=0.0, mc22=0.0;

	/* convert projection to matrixes */
	ma11 += a0; ma22 += a0;
	ma12 += -a1 * I; ma21 += -a1 * I;
	ma12 -= a2; ma21 += a2;
	ma11 += -a3 * I; ma22 -= -a3 * I;
	
	mb11 += b.a0; mb22 += b.a0;
	mb12 += -b.a1 * I; mb21 += -b.a1 * I;
	mb12 -= b.a2; mb21 += b.a2;
	mb11 += -b.a3 * I; mb22 -= -b.a3 * I;

	/* matrix multiply */
	mc11 = ma11 * mb11 + ma12 * mb21;
	mc12 = ma11 * mb12 + ma12 * mb22;
	mc21 = ma21 * mb11 + ma22 * mb21;
	mc22 = ma21 * mb12 + ma22 * mb22;

	*this = su2 (mc11, mc12, mc21, mc22);
	return *this;
}

su2&
su2 :: delta (double x, double y, double z)
{
	su2 that;
	that.a1 = x;
	that.a2 = y;
	that.a3 = z;
	that.a0 = sqrt (1.0 - x*x+y*y+z*z);
	*this *= that;

	return *this;
}

su2&
su2 :: conj (double x, double y, double z)
{
	su2 that;
	that.a1 = x;
	that.a2 = y;
	that.a3 = z;
	that.a0 = sqrt (1.0 - x*x+y*y+z*z);
	*this *= that;
	
	that.a1 = -x;
	that.a2 = -y;
	that.a3 = -z;
	that *= *this;
	*this = that;

	return *this;
}

void 
su2::print (void)
{
	printf ("A=(%g %g %g %g)\n", a0, a1, a2, a3);
	
	complex<double> ma11=0.0, ma12=0.0, ma21=0.0, ma22=0.0;

#if 0
	/* convert projection to matrix */
	ma11 += a0; ma22 += a0;
	ma12 += a1; ma21 += a1;
	ma12 -= a2 * I; ma21 += a2 *I;
	ma11 += a3; ma22 -= a3;

	printf (" = (%g+%gi   %g+%gi)\n", ma11.real(), ma11.imag(), ma12.real(), ma12.imag());
	printf ("   (%g+%gi   %g+%gi)\n", ma21.real(), ma21.imag(), ma22.real(), ma22.imag());
#endif

}
	

main () 
{
   double gb0 = 2.0e-3;  // gyromagnetic * magnetic field 
	
   long long imax = 48023;  // max iterations
   imax *= 1024*1024;
	imax = (long long) (10.0 / gb0);

	imax = 10000000;

   printf ("# \n");
   printf ("# su2 driven by soft photons \n");
   printf ("# gb0 = %g\n", gb0);
   printf ("# \n");
   printf ("# nprec	i	delta_i	initial_theta	theta\n");
   fflush (stdout);

	/* trun simulation several times */
	int loopmax = 100;
   for (int loop=0; loop<loopmax; loop++)
   {

		double theta = loop * 2.0 * M_PI /  (double) loopmax;
		
		// su2 dipole (0.0, 0.0, theta);
		su2 dipole (0.0, 1.6, 0.0);
		
      for (long long i=0; i<imax; i++)
      {
			double tmp;
			double magnet_theta;  // photon mag field orientation
			double magnet_phi;    // photon mag field orientation
   
			// theta uniform in upper hemisphere only
			// magnet_theta = (M_PI * rand()) / (2.0 * RAND_MAX);

#define HEMI_PHOTONS
#ifdef HEMI_PHOTONS
			// theta in upper hemisphere only, trailing off at equator
			magnet_theta = asin ((double) rand() / (double) RAND_MAX);
			magnet_phi = (2.0 * M_PI * rand()) / (1.0 * RAND_MAX);
   
			double mag_sin_theta = sin (magnet_theta);
			double mag_cos_theta = cos (magnet_theta);
			double mag_sin_phi = sin (magnet_phi);
			double mag_cos_phi = cos (magnet_phi);
   
			double z = gb0 * mag_cos_theta;
			double x = gb0 * mag_sin_theta * mag_cos_phi;
			double y = gb0 * mag_sin_theta * mag_sin_phi;
#endif

#ifdef RAND_PHOTONS
			double x = ((double) rand() / (double) RAND_MAX);
			double y = ((double) rand() / (double) RAND_MAX);
			double z = ((double) rand() / (double) RAND_MAX);

			x = gb0 * (2.0*x - 1.0);
			y = gb0 * (2.0*y - 1.0);
			// z = gb0 * (2.0*z - 1.0);
			z = gb0 * z;
#endif
			
			dipole.conj (x,y,z);
			printf ("i=%d ", i);
			dipole.print ();
		}
		printf ("l=%d ", loop);
		dipole.print ();
	}
}

/* ============================ END OF FILE ================ */
