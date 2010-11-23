#! /usr/bin/env perl
#
# Linas Novemver 2010
#
$fp = 3.14159265358979;
$fp *= $fp;
$fp *= 4;

$best = 0.99;

for ($n=2; $n < 1000000; $n++)
{
	$m = exp($fp / log($n));
	$fm = $m - floor($m);
	if ($fm < $est)
	{
		$best = $fm;
		print "$n  $m  $fm\n";
	}
}

