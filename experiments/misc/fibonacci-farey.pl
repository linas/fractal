#! /usr/bin/perl

# compute the limit of the mandelbrot set


$b[0] = 0;
$a[0] = 1;

$na = 1;
$nb = 1;

# compte the fibonacia L-system sequence
for ($i=0; $i<10; $i++) {

   for ($j=0; $j<$na; $j++) {
      $c[$j] = $a[$j];
   }

   for ($j=0; $j<$nb; $j++) {
      $c[$na+$j] = $b[$j];
   }

   $nc = $na+$nb;

   for ($j=0; $j<$na; $j++) {
      $b[$j] = $a[$j];
   }
   $nb = $na;

   for ($j=0; $j<$nc; $j++) {
      $a[$j] = $c[$j];
   }
   $na = $nc;
}

# construct the fraction

$t = 0.5;
$sum = 0;
for ($j=0;$j<56;$j++) {
   $k=$j+1;
   if ($a[$j] > 0) { $sum += $t; print "$k, ";}
   $t *= 0.5;
}
print "\n";

print "The fibonacci sequence number is $sum\n";


# count the number of bits in the binary expasion
$cnt[0] = 1;
$k=0;
if ($a[0] == 0) { $cnt[0] ++; }
for ($j=1;$j<$na;$j++) {
   if ($a[$j] != $a[$j-1]) { $k++; $cnt[$k] =1; } else { $cnt[$k]++; }
}

$kmax = $k;

# print the bit-count
for ($k=0; $k<=$kmax; $k++) {
   print "$cnt[$k]";
}
print "\n";


$r = 1;
for ($k=$kmax; $k>=0; $k--) {
   $r = $cnt[$k] + 1/$r;
}

$r = 1/$r;
print "The inverse farey nuer is $r\n";
