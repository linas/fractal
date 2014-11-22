#! /usr/bin/perl

for ($z=0; $z>-1; $z+=-0.05) {
   $f = 0;
   for ($i=0; $i<120; $i++) {
   
      $f = $z /  (1+$f);
   
   
   }
   print "$z $f\n";
}
