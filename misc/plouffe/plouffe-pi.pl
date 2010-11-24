#! /usr/bin/env perl
#
# quick-n-dirty sanity check of some pluffe sums
#
# Linas Vepstas November 2010
#

$pi = 3.14159265358979;
$q = exp(2*$pi);
$qn = $q;

$acc = 7*$pi*$pi*$pi / 180;

for($n=1; $n<20; $n++)
{
	$term = $n*$n*$n *($qn -1);
	$acc -= 2/$term;
	$qn *= $q;
}

print "Its $acc\n";
