
main () {

   float x, y;
   int n;

   x = 2.2; n = (int) x; if (0>x) n--; y = x-n; printf ("its %f %d %f \n", x, n, y);
   x = 1.2; n = (int) x; if (0>x) n--; y = x-n; printf ("its %f %d %f \n", x, n, y);
   x = 0.2; n = (int) x; if (0>x) n--; y = x-n; printf ("its %f %d %f \n", x, n, y);
   x = -0.2; n = (int) x; if (0>x) n--; y = x-n; printf ("its %f %d %f \n", x, n, y);
   x = -1.2; n = (int) x; if (0>x) n--; y = x-n; printf ("its %f %d %f \n", x, n, y);
   x = -2.2; n = (int) x; if (0>x) n--; y = x-n; printf ("its %f %d %f \n", x, n, y);
}

