#! /usr/bin/perl

# compute the limit of the mandelbrot set

# first, compte the limit of the power-of-2 sequence.
$a = 1;
for ($i=1; $i<11; $i++) {
   $tp = exp (($i-1)*log(2));
   $tpp = exp ($tp*log(2));
   $deno = $tpp +1;

   $r = $a/$deno;
   print "$i  $a / $deno =$r\n";

   $a = $a*($tpp-1) +1;
}
$t = $r;

$r = 2;
$r = 2 + 1/$r;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 1/$r;

print "inverse fairy of thue-morse is $r\n";

$r = 2;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 1 + 1/$r;
$r = 1 + 1/$r;
$r = 2 + 1/$r;
$r = 2 + 1/$r;
$r = 1/$r;

print "inverse fairy of thue-morse is $r\n";

# $t = 0.412454033640107;

# get the binary expansion of the thue-morse number
$max = 48;
for ($i=0;$i<$max;$i++) {
    if ($t<0.5) {print "0"; } else {print "1"; }
    if ($t<0.5) {$b[$i]=0; } else { $b[$i]=1; }

    $t *=2;
    if ($t > 1.0) {$t -=1; }
}
print "\n";

# count the number of digits in the binary expansion
$c[0] = 1;
$j=0;
if ($b[0] == 0) { $c[0] ++; }
for ($i=1;$i<$max;$i++) {
   if ($b[$i] != $b[$i-1]) { $j++; $c[$j] =1; } else { $c[$j]++; }
}

$jmax = $j;

# substitue in 3/7, 4/7 for 0.1
$bud = 0;
$scale=0.125;
for ($j=0; $j<=$jmax; $j++) {
   if ($b[$j]==0) { $bud += 3*$scale; } else { $bud +=4*$scale; }
   $scale *=0.125;
   print "Buddy $bud\n";
}
print "\n";

#print the number of digits
for ($j=0; $j<=$jmax; $j++) {
   print "$c[$j]";
}
print "\n";


# convert to continued fraction
$r = 1;
for ($j=$jmax; $j>=0; $j--) {
   $r = $c[$j] + 1/$r;
}

$r = 1/$r;
print " high-precision its $r\n";
