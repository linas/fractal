extern double drand48();

void main () 
{
   unsigned int	i;
   double		y;
   int   		big, z;

   printf ("whoops!\n");
   
   big = 0;
   for (i=0; i<3000; i++) {
      y = drand48 ();
      z = rand ();
      if ( z > 16382) big++;
      printf ("rand %d %d %f\n",i-big, big, y);
   }
   
   printf ("filled in the circle \n");
}
