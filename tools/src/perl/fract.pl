#! /usr/bin/perl
# some track from rational expressions 
# (a+b*sin)/(c+d*sin) and (a+b*cos)/(c+d*cos)

$x=0;
$y=0;

$a = 13;
$b = 3;
$c = 30;
$d = 7;

$a = 61;
$c = 142;

$a = 3;
$b = 1;
$c = 7;
$d = 2;

$n = 130;

for ($i=0; $i<$n; $i++) {
   $th = 2.0*3.14159265358979 *$i/$n;
   $co = cos ($th);
   $si = sin ($th);
   $x = $a + $b * $co;
   $x /= $c + $d * $co;
   $y = $a + $b * $si;
   $y /= $c + $d * $si;
   print "$i $th $x $y\n";

}
