#! /usr/bin/perl

# coherent eigenstates of bernoulli map

$nsteps = 50;
$l = 0;

$nmax = 100;
$z = 0.75;

$twopi= 2.0 * 3.14159265358979;
$tl = 2*$l +1;

for ($i=0; $i<=$nsteps; $i ++)
{
   $x = $i/$nsteps;

   $e = 0;
   $tpn = 1;
   $zn = 1;
   for ($n=0; $n < $nmax; $n++)
   {
      $c = cos ($twopi * $tl * $tpn * $x);
      $e += $zn *$c;
      # print "          $tpn    $c    $e \n";

      $tpn *= 2;
      $zn *= $z;
   }
   print "$x     $e\n";

}
