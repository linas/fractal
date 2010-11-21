#! /usr/bin/env perl
#
# a few sums
#

$gkw = 0.0;
$dyadic = 0.0;

my $tn = 1.0;
for ($n=1; $n<40; $n++)
{
	$dyadic += $n / $tn;
	$tn *= 2.0;
}
$dyadic *= log(2.0);
print "Dyadic entropy is $dyadic\n";

