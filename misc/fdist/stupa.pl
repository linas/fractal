#! /usr/bin/env perl

$n=10;
while(1)
{
	$xlo = 0.5 ** $n;
	$xhi = 2 * $xlo;
	$hgt = 4 * $xlo;
	print "$xlo	0.0\n";
	print "$xlo	$hgt\n";
	print "$xhi	$hgt\n";
	print "$xhi	0.0\n";

	$n--;
	if ($n == 1) { last; }
}

$n=2;
while(1)
{
	$xlo = 0.5 ** $n;
	$xhi = 2 * $xlo;
	$hgt = 4 * $xlo;

	$ylo = 1.0 - $xhi;
	$yhi = 1.0 - $xlo;
	print "$ylo	0.0\n";
	print "$ylo	$hgt\n";
	print "$yhi	$hgt\n";
	print "$yhi	0.0\n";

	$n++;
	if ($n == 10) { last; }
}
