#/usr/bin/env perl
#
# simple hyperconvergent series -- Euler's series 
# (exponential integral)

$z = -0.1

$prod= $z*$z;
$fact = 1;
$sum = $z + $z*$z;

for ($k=2; $k<40, $k++)
{
	$fact *= $k;
	$prod *= $z;
	$sum += $fact * $prod;
	print "$k $sum\n";

}

