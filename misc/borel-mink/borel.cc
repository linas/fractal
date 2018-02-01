
/**
 * borel.cc
 *
 * Borel transform, or Riesz–Herglotz transform on the Minkowski measure
 *
 * Linas Vepstas February 2018
 */

#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include "question.h"

#ifndef DEBUG
#include "brat.h"
#endif

// Goddamnd fucked up C++ complex numbers
#define COMPLEX std::complex<double>

/**
 * The  Riesz–Herglotz transform on the minkowski measure.
 * `level` is the iteration level to recurse to.
 *
 * Raw version, uncahced, slow.
 */
COMPLEX riesz_mink_nocache(COMPLEX z, int level)
{
	unsigned __int128 p, q, pm, qm, pmid, qmid;

	// Golden ratio
	double phi = 0.5 * (1.0 + sqrt(5.0));

	// Growth rate
	double grow = pow(0.5 * phi * phi, level) * 0.2 * phi * phi * phi;

	COMPLEX acc = 0.0;
	double quest = 0.0;

	double norm = 1.0 / (double)(1<<level);
	for (int i=0; i< (1<<level); i++)
	{
		double x = ((double) i + 0.5) * norm;
		stern_brocot_tree128(i, level, &p, &q);
		stern_brocot_tree128(i+1, level, &pm, &qm);

		// a and b are the upper and lower endpoints of the interval
		double a = ((double) p) / (double) q;
		double b = ((double) pm) / (double) qm;
		stern_brocot_tree128(2*i+1, level+1, &pmid, &qmid);

		// The midpoint of the interval.
		double y = ((double) pmid) / (double) qmid;

		// unsigned long det = pm * q - p * qm;

		// measure: the value of the measure at the midpoint.
		double measure = norm * q * qm;
		measure /= grow;

		// The integral of the masure: should be the question mark.
		// up to a scale factor
		quest += measure * (b-a);

		// You stupid plick: measure*(b-a) == const

		// theta ranges -pi to pi
		double theta = (2.0*y - 1.0) * M_PI;
		COMPLEX it = I*theta;
		COMPLEX eit = exp(it);

		// The kernel
		COMPLEX krn = (eit + z) / (eit - z);

		acc += krn * measure * (b-a);

		printf("%d	%g	%g	%g	%g	%g	%g\n", i, x, y, measure, quest, real(acc), imag(acc));
	}

	// Normalize
	acc /= quest;

	acc /= 2.0 * M_PI;

	return acc;
}

/**
 * The  Riesz–Herglotz transform on the minkowski measure.
 * `level` is the iteration level to recurse to.
 */
COMPLEX riesz_mink(COMPLEX z, int level)
{
	int npts = (1<<level);
	static bool init = false;
	static double* midpoint = nullptr;
	static double* measure = nullptr;
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	
	if (not init)
	{
		pthread_mutex_lock(&mutex);
		if (not init) {

		unsigned __int128 p, q, pm, qm, pmid, qmid;

		// Golden ratio
		double phi = 0.5 * (1.0 + sqrt(5.0));

		// Growth rate
		double grow = pow(0.5 * phi * phi, level) * 0.2 * phi * phi * phi;

		midpoint = (double*) malloc(npts*sizeof(double));
		measure = (double*) malloc(npts*sizeof(double));
		double norm = 1.0 / (double) npts;
		double quest = 0.0;
		for (int i=0; i<npts; i++)
		{
			stern_brocot_tree128(i, level, &p, &q);
			stern_brocot_tree128(i+1, level, &pm, &qm);

			// a and b are the upper and lower endpoints of the interval
			double a = ((double) p) / (double) q;
			double b = ((double) pm) / (double) qm;
			stern_brocot_tree128(2*i+1, level+1, &pmid, &qmid);

			// The midpoint of the interval.
			double y = ((double) pmid) / (double) qmid;

			// measure: the value of the measure at the midpoint.
			double meas = norm * q * qm;
			meas /= grow;
			// meas *= (b-a);

			// The integral of the masure: should be the question mark.
			// up to a scale factor
			quest += meas;
			midpoint[i] = y;
			measure[i] = meas;
		}

		// Normalize
		for (int i=0; i<npts; i++)
		{
			measure[i] /= quest;
		}

		printf("Done with init\n");
		init = true;
		}
		pthread_mutex_unlock(&mutex);
	}

	COMPLEX acc = 0.0;
	double norm = 1.0 / (double) npts;
	for (int i=0; i< (1<<level); i++)
	{
		double x = ((double) i + 0.5) * norm;
		// double y = midpoint[i];
		// theta ranges -pi to pi
		// double theta = (2.0*y - 1.0) * M_PI;
		double theta = (2.0*x - 1.0) * M_PI;
		COMPLEX it = I*theta;
		COMPLEX eit = exp(it);

		// The kernel
		COMPLEX krn = (eit + z) / (eit - z);

		acc += krn * measure[i];
	}

	// Normalize
	acc /= 2.0 * M_PI;

	return acc;
}

#ifndef DEBUG
double xform(double re_q, double im_q, int itermax, double parm)
{
	COMPLEX z = re_q + I * im_q;
	COMPLEX g = riesz_mink(z, itermax);

#define MAG
#ifdef MAG
	double rv = abs(g);
#else
	double rv = arg(g);
	rv += M_PI;
	rv /= 2.0 * M_PI;
#endif

	return rv;
}

DECL_MAKE_HEIGHT(xform)
#endif // DEBUG

#ifdef DEBUG
int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s level\n", argv[0]);
		exit(1);
	}

	int level = atoi(argv[1]);

	COMPLEX z = I*1.0;
	z = 0.0;
	COMPLEX g = riesz_mink_nocache(z, level);
	g *= 2.0*M_PI;

	printf("# its %g + i%g\n", real(g), imag(g));
}
#endif
