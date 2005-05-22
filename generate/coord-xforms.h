/*
 * coord-xforms.h
 *
 * FUNCTION:
 * provide coordinate transforms from upper half plane to the 
 * poincare disk and the q-series disk, and thier inverses.
 *
 * HISTORY:
 * New, May 2005
 */

void poincare_disk_to_plane_coords (double x, double y, 
                                    double *px, double *py);

void plane_to_poincare_disk_coords (double x, double y, 
                                    double *px, double *py);

void plane_to_q_disk_coords (double tau_re, double tau_im, 
                             double *px, double *py);

/* map from q-series coords to the upper half-plane */
void q_disk_to_plane_coords (double qre, double qim,
                             double *px, double *py);

/* Apply mobius transform (ab)(cd) = (az+b)/(cz_d) */
void mobius_xform (double a, double b, double c, double d,
              double tau_re, double tau_im, 
              double *px, double *py);

