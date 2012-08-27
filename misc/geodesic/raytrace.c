/*
 * raytrace.c
 *
 * Explore geodesics on the fundamental domain of the modular group
 * by ray-tracing.  Arbitrary-precision version need to make sure we 
 * get the right number of digits.  Geodesics are circles with origin
 * on the imaginary axis.  The path that a geodesic takes is a double-
 * sides string of S and T values, which can be converted to a binary
 * string via Minkowski ? function.
 *
 * Linas Vepstas August 2012
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef int bool;

typedef struct 
{
	double x;   // position
	double y;
	double vx;  // velocity
	double vy;
} ray_t;

/* Return S or T or N depending on whther the ray exited on the bottom,
 * or on the right, or on the left of the fundamental domain bounded by
 * two vertical lines from +/- 1/2 +i rho where rho=
 * and a circular arc stretching between them, of radius r=
 */
char bounce (const ray_t in, ray_t* out)
{
	static double rho = 0.5 * sqrt(3.0);

	// Calculate center of circle.  Its at y=0, x=center.
	double center = in.x + in.y * in.vy / in.vx;
printf("--------\ncenter at %g\n", center);
	// Square of the radius of the circle
	double radius_sq = (in.x - center)*(in.x-center) + in.y*in.y;

	char exit_code = 'S';
	double x_exit = (1.0 - radius_sq + center*center) / (2.0 * center);

	// intersect is true if the geodesic and the unit circle interest.
	// If intersect, then can take sqrt to get valid intersection point.
	bool intersect = (x_exit-center)*(x_exit-center) < radius_sq;

	// If intersect is true, then the dot product tells us whether the
	// trajectory is incoming, or outgoing.  If its outgoing, then its
	// a bottom exit, else its a side exit.  dot is true for bottem
	// entry.
	bool dot = 0.0 < (in.x* in.vx + in.y *in.vy);
	if (!intersect || (intersect && dot))
	{
		if (in.vx > 0.0)
		{
			// Right side exit, not bottom exit
			x_exit = 0.5;
			exit_code = 'T';
		}
		else 
		{
			// Left side exit, not bottom exit
			x_exit = -0.5;
			exit_code = 'N';
		}
	}

	double y_exit = sqrt(radius_sq - (x_exit-center)*(x_exit-center));
	if (('S' != exit_code) && (y_exit < rho))
	{
		fprintf(stderr, "Error: bad y value for side exit\n");
		exit(1);
	}

printf(" %c exit x=%g y=%g\n", exit_code, x_exit, y_exit);

	double tan_exit = -(x_exit-center)/y_exit;
	double vx_exit = sqrt(1.0 / (1.0 + tan_exit*tan_exit));
	double vy_exit = sqrt (1.0 - vx_exit*vx_exit);
	if (tan_exit < 0.0) vy_exit = -vy_exit;

printf("velc=%g %g %g\n", vx_exit, vy_exit, vx_exit*vx_exit+vy_exit*vy_exit);

	// now reflect
	switch (exit_code)
	{
		case 'T':
			out->x = -0.5;
			out->y = y_exit;

			out->vx = vx_exit;
			out->vy = vy_exit;
			break;
		case 'N':
			out->x = 0.5;
			out->y = y_exit;

			out->vx = vx_exit;
			out->vy = vy_exit;
			break;
		case 'S':
		{
			out->x = -x_exit;
			out->y = y_exit;

			double deno = 2.0 * center * vx_exit;
			out->vx = -vx_exit + x_exit * deno;
			out->vy = vy_exit - y_exit * deno;
			break;
		}
	}
vx_exit = out->vx;
vy_exit = out->vy;
printf("final velc=%g %g %g\n", vx_exit, vy_exit, vx_exit*vx_exit+vy_exit*vy_exit);

	return exit_code;
}

void sequence()
{
	ray_t in, out;
	in.x = 0.0;
	in.y = 2.0;
	double theta = 0.1;
	in.vx = cos(theta);
	in.vy = sin(theta);	

printf("start %g %g\n", in.vx, in.vy);
	int i;
	for (i=0; i<10; i++)
	{
		bounce(in, &out);
		in = out;
	}
}

int main(int argc, char * argv[]) 
{
	sequence();
}




