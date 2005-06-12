
/*
 * geo-lib.h
 *
 * Common routines for computation of the lengths and energies of 
 * geodesics on the riemann surface/fundamental domain. 
 */

/* return number of translations befor circular arc hits the
 * scallopped inversion boundary */
int get_n_of_rho (double rho);

/* x coordinate where geodesic hits the scalloped boundary */
double geo_x (double rho);

/* length and energy of geodesic between scalloped boundaries */
double geo_length (double rho);
double geo_energy (double rho);
