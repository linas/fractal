#! /usr/bin/env perl

$binsz = 0.2;

while (<>) {
   chop;
   ($nprec, $i, $j, $init, $theta) = split /	/;
   
   $n = $theta / $binsz;
   $n = int ($n);

   $bin[$n] ++;
}

for ($i=0;$i<3.14/$binsz; $i++) {
   $theta = $i * $binsz;
   print "$theta	$bin[$i]\n";
}

