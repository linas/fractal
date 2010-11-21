#! /usr/bin/env perl
#
# a few sums
#

$gkw = 0.0;
$dyadic = 0.0;

for ($n=1; $n<40; $n++)
{
	$dyadic += $n / pow(2.0, $n);
}
$dyadic *= log(2.0);
print "Dyadic entropy is $dyadic\n";

for ($n=1; $n<80000; $n++)
{
	$gkw += log($n) / $n;
}

print "GKW entropy is $gkw\n";
