#! /usr/bin/perl
#
# Compute power spectrum of the zero locations as a function of N.
#
# Linas Vepstas December 2010
#

# read in data
$linecnt = 0;
while (<>)
{
	if (/#/) {next; }
	($n, $re, $im) = split;

	$rdata[$linecnt] = $re;
	$idata[$linecnt] = $im;

	$linecnt ++;
}

# compute power 
$freqmin = 0.005;
$freqmax = 1.0;

# for ($freq = $freqmin; $freq < $freqmax; $freq+= 0.005) {
#	$period = 2*3.141592653 / $freq;
for ($period=6; $period<620; $period++)
{
	$freq = 2*3.141592653 / $period;
	$resumcos = 0;
	$resumsin = 0;
	$imsumcos = 0;
	$imsumsin = 0;
	for ($i=0; $i<$linecnt; $i++)
	{
		$resumcos += $rdata[$i] * cos($i*$freq);
		$resumsin += $rdata[$i] * sin($i*$freq);
		$imsumcos += $idata[$i] * cos($i*$freq);
		$imsumsin += $idata[$i] * sin($i*$freq);
	}

	$repow = $resumcos*$resumcos + $resumsin * $resumsin;
	$repow = sqrt($repow);
	$impow = sqrt($imsumcos*$imsumcos + $imsumsin * $imsumsin);

	# $repow /= $linecnt;
	# $impow /= $linecnt;

	print "$freq	$period	$repow	$impow\n";
}
