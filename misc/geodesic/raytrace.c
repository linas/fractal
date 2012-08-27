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

typedef struct 
{
  double x;   // position
  double y;
  double xv;  // velocity
  double yv;
} ray_t;

/* Return S or T or N depending on whther the ray exited on the bottom,
 * or on the right, or on the left of the fundamental domain bounded by
 * two vertical lines from +/- 1/2 +i rho where rho=
 * and a circular arc stretching between them, of radius r=
 */
char bounce (const ray_t in, ray_t* out)
{
    return 'S';
}

int main(int argc, char * argv[]) 
{
}




