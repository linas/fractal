#! /usr/bin/env perl

$intervals = 10;

$binsz = 3.141592653 / $intervals;

for ($i=0;$i<$intervals; $i++) {
   $bin[$i] = 0;
}

while (<>) {
   chop;
   ($nprec, $i, $j, $init, $theta) = split /	/;

   if (($j < 50000000) && ($j > 0)) { next; }   

   $n = $theta / $binsz;
   $n = int ($n);

   $bin[$n] ++;
}

for ($i=0;$i<$intervals; $i++) {
   $theta = ($i+0.5) * $binsz;
   print "$theta	$bin[$i]\n";
}

