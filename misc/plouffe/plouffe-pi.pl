#! /usr/bin/env perl
#
# quick-n-dirty sanity check of some pluffe sums
#
# Linas Vepstas November 2010
#

$pi = 3.14159265358979;
$q = exp($pi);
$qn = $q * $q;

$acc = 7*$pi*$pi*$pi / 180;

for($n=1; $n<20; $n++)
{
	$term = $n*$n*$n *($qn -1);
	$acc -= 2/$term;
	$qn *= $q*$q;
}

print "Its $acc\n";

# -------------------------------------
$qn = $q;
$qfn = $q*$q*$q*$q;

$acc = 9*$pi*$pi*$pi / 900;

for($n=1; $n<20; $n++)
{
	$term = 4/($qn -1) + 1/($qfn -1);
	$term /= $n*$n*$n;
	$acc -= 0.4 * $term;
	$qn *= $q;
	$qfn *= $q*$q*$q*$q;
}

print "Its $acc\n";



