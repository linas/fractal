#! /usr/bin/perl

$x=0;
$y=0;
$omega = 0.1;
$K= 0.3;

for ($i=0; $i<130; $i++) {
   $x = $x + $omega + $K * sin (2.0*3.14159265358979 *$x);
   $y = $y + $omega - $K * sin (2.0*3.14159265358979 *$y);
   $s = $x + $y;
   print "$i $x $y $s\n";

}
