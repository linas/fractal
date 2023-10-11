#!/usr/bin/perl

for ($r=1.75; $r>0; $r -=0.05)
{
	$r4 = $r*$r*$r*$r;


	$m = exp (1-$r4*$r4) / $r4;

	print "$r 	$m\n"


}
