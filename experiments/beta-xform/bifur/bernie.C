/*
 * bernie.C
 *
 * Draw Feigenbaum-style bifurcation diagrams for the beta-map
 * and lots of other maps. These diagrams collect up the "invariant
 * measure" as they are iterated.
 *
 * Some of the maps replace multiplication by multiplication-like
 * bit-mangling algorithms that are simpler (e.g. don't propagate
 * a carry bit).
 *
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bitops.h"
#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 */

#if 0
static double xprod(double x, double y, int itermax, double param)
{
	// return mult_xor(x, y);
	return x*y;
}

DECL_MAKE_HEIGHT (xprod);
#endif

double downshift(double x, double K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (x - 0.5);
	}
	return K*x;
}

// A zig-zag map, with diagonal joining the two branches
double island(double x, double K, double epsilon)
{
	double beta = 2.0 * K;
	if (0.5+epsilon <= x)
	{
		return beta * (x - 0.5);
	}
	if (x < 0.5-epsilon)
	{
		return beta*x;
	}
	return 0.25*beta* (1.0 + (1.0 - 4.0*epsilon) * (0.5 - x) / epsilon);
}

// Like above, but with an S-curve for the middle segment.
double lost_island(double x, double K, double epsilon)
{
	double beta = 2.0 * K;
	if (0.5+epsilon <= x)
	{
		return beta * (x - 0.5);
	}
	if (x < 0.5-epsilon)
	{
		return beta*x;
	}

	// om runs from  -1 to 1;
	double om = (x - 0.5) / epsilon;

	// kos runs from -1 to 1
	double kos = sin(M_PI * 0.5 * om);
	// double kos = -sin(M_PI * 0.5 * om);

	// S-curve interpolation between top and bottom.
	return 0.25*beta - beta*(0.25-epsilon) * kos;
}

// Like above, but with soft top and hard-edged bottom
double hard_island(double x, double K, double epsilon)
{
	double beta = 2.0 * K;
	if (0.5+epsilon <= x)
	{
		return beta * (x - 0.5);
	}
	if (x < 0.5-epsilon)
	{
		return beta*x;
	}

	// om runs from  -1 to 1;
	double om = (x - 0.5) / epsilon;

	// kos still runs from 1 to -1
	double kos = cos(M_PI * 0.25 * (om+1.0));
	kos = 2.0 * (kos - 0.5);

	// First-quarter cosine interpolation between top and bottom.
	return 0.25*beta + beta*(0.25-epsilon) * kos;
}

double sign(double x)
{
	if (x < 0.0) return -1.0;
	if (0.0 < x) return 1.0;
	return 0.0;
}

// Like above, but with half of an S-curve for the middle segment
double ess_island(double x, double K, double epsilon)
{
	double beta = 2.0 * K;
	if (0.5+epsilon <= x)
	{
		return beta * (x - 0.5);
	}
	if (x < 0.5-epsilon)
	{
		return beta*x;
	}

	// om runs from  -1 to 1;
	double om = (x - 0.5) / epsilon;

	// kos runs from -1 to 1
	// double kos = pow(om, 3);

	// With this sign convention, the map is continuous, and
	// has an S-curve shape.
	// double kos = om;
	// double kos = sign(om) * om*om;
	// double kos = om*om*om;
	// double kos = sign(om) * om*om*om*om;
	double kos = om*om*om*om*om;

	// With this sign convention, the middle segment is increasing,
	// and the map consists of three disjoint segments.
	// double kos = -om;
	// double kos = -sign(om) * om*om;
	// double kos = -om*om*om;
	// double kos = -sign(om) * om*om*om*om;
	// double kos = -om*om*om*om*om;

	return 0.25*beta - beta*(0.25-epsilon) * kos;
}

// Like above, but the location of the kink adjustable
double kink_island(double x, double K, double epsilon, double alpha)
{
	double beta = 2.0 * K;
	if (0.5+epsilon <= x)
	{
		return beta * (x - 0.5);
	}
	if (x < 0.5-epsilon)
	{
		return beta*x;
	}

	// om runs from  -1 to 1;
	double om = (x - 0.5) / epsilon;

	// kos runs from -1 to 1
	// double kos = pow(om, 3);

	// With this sign convention, the map is continuous, and
	// has an S-curve shape.
	// double kos = om;
	// double kos = sign(om) * om*om;
	double kos = om*om*om;
	// double kos = sign(om) * om*om*om*om;
	// double kos = om*om*om*om*om;

	// With this sign convention, the middle segment is increasing,
	// and the map consists of three disjoint segments.
	// double kos = -om;
	// double kos = -sign(om) * om*om;
	// double kos = -om*om*om;
	// double kos = -sign(om) * om*om*om*om;
	// double kos = -om*om*om*om*om;

	kos = fabs(kos);
	if (x<0.5) kos = alpha + (1.0-alpha)*kos;
	else kos = alpha - (1.0+alpha)*kos;

	return 0.25*beta - beta*(0.25-epsilon) * kos;
}

double nocarry(double x, double K)
{
	// Notes:
	// 1) The mult_xor function is technically incorrect;
	//    its too large by a factor of two.
	//    That's why we can skip multiplying K by two.
	// 2) The mult_xor function just shuffles, so we post-multiply by
	//    two, to get a return value between zero and one.
	if (0.5 <= x)
	{
		return 2.0 * mult_xor (K, (x - 0.5));
	}
	return 2.0*mult_xor(K, x);
}

int clamp(int n)
{
	// return n+1; // mang-carry-1-more

	if (0 == n) return n;
	return n-1;  // mang-carry-1-less

	// return n%2;     // mang-carry-mod-2.png aka shift-XOR
	// return 1-n%2;   // mang-flip-mod.png

	// return n%3;

	// int max = 1;    // mang-carry-1
	// int max = 2;    // mang-carry-2
	int max = 3;    // mang-carry-3
	// int max = 6;
	if (max < n) return max;
	return n;
}

double mangle_carry(double x, double K)
{
	if (0.5 <= x)
	{
		x -= 0.5;
	}
	return 2.0* mangle_multiply(K, x, clamp);
}

double tent(double x, double K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (1.0 - x);
	}
	return K*x;
}

double notent(double x, double K)
{
	if (0.5 <= x)
	{
		return 2.0*mult_xor (K, (1.0 - x));
	}
	return 2.0*mult_xor(K, x);
}

double feig(double x, double K)
{
	K *= 4.0;
	return K * x * (1.0 - x);
}

// Smoothly interpolates between tent (for epsi=0) and feigenbaum
// (for epsi=1)
double bulgetent(double x, double K, double epsi)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (1.0 - x + epsi * (-1.0 + 3.0*x - 2.0*x*x));
	}
	return K*x * (1.0 + epsi * (1.0-2.0*x));
}

double nofeig(double x, double K)
{
	// Below is intersting/crazy
	// return mult_xor(K, 4.0 * x * (1.0 - x));

	// Below is same as above. Multiplying by 2 commutes with mult_xor
	// return 2.0 * mult_xor(K, 2.0 * x * (1.0 - x));
	return 4.0 * mult_xor(K, x * (1.0 - x));

	// Below is totally uniform
	// return 2.0* mult_xor(K, 2.0 * mult_xor(x, (1.0 - x)));
}

// Superimpose a sine-wave like bulge on the straight line.
double bulge(double x, double K, double epsi)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		double y = epsi * (1.0-x) * (0.5-x);
		x += y;
		return K * (x - 0.5);
	}
	double y = epsi * x * (0.5-x);
	x += y;
	return K*x;
}


static void bifurcation_diagram (float *array,
                                 int array_size,
                                 double x_center,
                                 double x_width,
                                 double K,
                                 int itermax,
                                 double omega)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int cnt=0;

	double Korg = 0.5 + 0.5*K;
	// double eps = 0.5* (1.0-K);
	// double eps = 0.04;
	// double eps = omega;
	// double alpha = omega;
	for (int j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		t = rand();
		t /= RAND_MAX;
		t -= 0.5;
		t /= 800.0; // 800 pixels tall
		t *= 0.5; // K runs from 0.5 to 1.0
		t *= x_width; // in case its zoomed.
		K = Korg + t;

		// For Feigenbaum only, runs from 1.75 to 2.0 at top
		t *= 0.25;
		K = 1.0 - 0.25*(1.0 - Korg) + t;

		/* OK, now start iterating the benoulli map */
		for (int iter=0; iter < 1250; iter++)
		{
			// x = downshift(x, K);
			// x = nocarry(x, K);
			// x = tent(x, K);
			// x = notent(x, K);
			x = feig(x, K);
			// x = nofeig(x, K);
			// x = mangle_carry(x, K);
			// x = island(x, K, eps);
			// x = lost_island(x, K, eps);
			// x = hard_island(x, K, eps);
			// x = ess_island(x, K, eps);
			// x = kink_island(x, K, eps, alpha);
			// x = bulge(x, K, eps);
			// x = bulgetent(x, K, eps);

			double en = array_size * (x-floor(x));
			int n = en;
			if (0 > n) n = 0;
			if (n >= array_size) n = array_size-1;
			array[n] += 1.0;
			cnt ++;
		}
	}

	double wi = 0;
	double cut = 0.02 * cnt / ((double) array_size);
	for (int j=0; j<array_size; j++)
	{
		if (cut < array[j]) wi += 1.0;
	}

	// double norm = ((double) array_size) / ((double) cnt);
	double norm = wi / ((double) cnt);
	for (int j=0; j<array_size; j++)
		array[j] *= norm;

#ifdef MARKUP_RETICULE
	double beta = 2.0 * Korg;
	double pix = 0.45 / array_size;
	int left = (0.5 + 0.5* (beta-1)) * array_size;
	left -= 54;
	double phi;
/*
	phi = sqrt(2.0);
	if (fabs(beta-phi)<pix)
		for (int j=left; j<array_size; j++) {array[j] = 1.5; }
*/

	phi = sqrt(2.0/(1.0-alpha));
	if (fabs(beta-phi)<pix)
		for (int j=left; j<array_size; j++) {array[j] = 1.5; }

	phi = pow(2.0/(1.0-alpha), 0.33333);
	if (fabs(beta-phi)<pix)
		for (int j=left; j<array_size; j++) {array[j] = 1.5; }

	phi = pow(2.0/(1.0-alpha), 0.25);
	if (fabs(beta-phi)<pix)
		for (int j=left; j<array_size; j++) {array[j] = 1.5; }
#endif

#ifdef GOLDEN_ISLAND
	double beta = 2.0 * Korg;
	double pix = 0.45 / array_size;
	int left = (0.5 + 0.5* (beta-1)) * array_size;
	// left -= 234; // For golden at epsilon=0.15
	// left -= 145; // For golden at epsilon=0.10
	// left -= 54; // For golden at epsilon=0.04
	left -= 54;

	double corner = (1 + 2.0*eps) / (1 - 2.0*eps);

	double phi;
	phi = 1.618034;
	phi = corner + (2-corner)*(phi-1);
	if (fabs(beta-phi)<pix)
		for (int j=left; j<array_size; j++) {array[j] = 1.5; }

#if 0
	phi = 1.4655;
	// phi = corner + (2-corner)*(phi-1);
	if (fabs(beta-phi)<pix)
		for (int j=left; j<array_size; j++) {array[j] = 1.5; }

	phi = 1.3803;
	// phi = corner + (2-corner)*(phi-1);
	if (fabs(beta-phi)<pix)
		for (int j=left; j<array_size; j++) {array[j] = 1.5; }

	phi = 1.3247;
	// phi = corner + (2-corner)*(phi-1);
	if (fabs(beta-phi)<pix)
		for (int j=left; j<array_size; j++) {array[j] = 1.5; }
#endif
#endif

#ifdef KORNER
	// Draw lines indicating the left corner of the
	// island and hard-island maps.
	int cleft = eps *(1+2*eps) / (1-2*eps) * array_size;
	array[cleft] = 1.5;

	phi = (1 + 2.0*eps) / (1 - 2.0*eps);
	if (fabs(beta-phi)<pix)
		for (int j=0; j<cleft; j++) {array[j] = 1.5; }
#endif
}

DECL_MAKE_BIFUR(bifurcation_diagram)
