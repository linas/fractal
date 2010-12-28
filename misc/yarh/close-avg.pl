#! /usr/bin/perl
#
# Take a gaussian average of the zero locations as a function of N.
#
# Linas Vepstas December 2010
#

$sigma = 2;
$width = 6 * $sigma;

for ($i=-$width; $i<$width; $i++)
{
	$x = $i/$sigma;
	$gauss[$width+$i] = exp(-$x*$x);
	# print "Gauss $i $gauss[$width+$i]\n";
}

for ($i=-$width; $i<$width; $i++)
{
	$rdata[$width+$i] = 0;
	$idata[$width+$i] = 0;
}

$linecnt = 0;

while (<>)
{
	for ($i=-$width; $i<$width; $i++)
	{
		$rdata[$width+$i] = $rdata[$width+$i+1];
		$idata[$width+$i] = $idata[$width+$i+1];
	}
	($n, $re, $im) = split;

	$rdata[$width+$width] = $re;
	$idata[$width+$width] = $ie;

	$linecnt ++;
	if ($linecnt < 2*$width) next;

	$ravg = 0;
	$iavg = 0;
	for ($i=-$width; $i<$width; $i++)
	{
		$ravg += $rdata[$width+$i] * $gauss[$width+$i];
		$iavg += $idata[$width+$i] * $gauss[$width+$i];
	}

	$i = $n-$width;
	print "$i	$ravg	$iavg\n";
}
