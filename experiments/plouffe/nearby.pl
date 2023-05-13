#! /usr/bin/env perl
#
# Linas Novemver 2010
#

use POSIX qw/floor/;

$fp = 3.14159265358979;
$fp *= $fp;
$fp *= 4;

$best = 0.99;

for ($n=5; $n < 2000000; $n++)
{
	$m = exp($fp / log($n));
	$fm = $m - floor($m);
	if ($fm < $best)
	{
		$best = $fm;

		$m = floor($m);
		$x = log($n)/(2*3.14159265358979);
		$y = (2*3.14159265358979) / log($m);
		$d = sqrt($x*$y);
		$en = exp((2*3.14159265358979)*$d);
		$em = exp((2*3.14159265358979)/$d);
		print "$n  $m  $fm $d \t$en $em \n";
	}
}

