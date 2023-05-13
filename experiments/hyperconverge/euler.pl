#!/usr/bin/env perl
#
# simple hyperconvergent series -- Euler's series 
# (exponential integral)
# abramowitz & stegun eqn 5.1.51

$z = - 0.05;

$prod= $z*$z;
$fact = 1;
$sum = $z + $z*$z;

for ($k=2; $k<40; $k++)
{
	$fact *= $k;
	$prod *= $z;
	$sum += $fact * $prod;
	print "$k $sum\n"; 

}

$prod= $z*$z;
$fact = 1;
$sum = $z + $z*$z;

for ($k=2; $k<20; $k+=2)
{
	$fact *= $k;
	$prod *= $z;
	$sum += $fact * $prod * (1 + $k*$z);

	$tmp = (1 + $k*$z);
	$term = $fact * $prod * $tmp;
	print "$k $sum $tmp $term \n"; 

	$fact *= $k;
	$prod *= $z;
}


$e = &pow(2.7, 1.0/$z);

print "its $e \n";
