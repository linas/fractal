#! /usr/bin/env perl
#
# Run processes to bulk convert polylog movie data files to PNG files.
#

# fps == frames per second
# Must match setting for script that initiall generated the frames.
$fps = 30;
$step = 1/ $fps;

$taumax = 300;  # max tau value but also running time in secs

# $fpref == Filename prefix
# $deci == decimal prefix for output.
sub cvt_range($imin, $imax, $deci, $fpref) {

	for ($i=$imin; $i<$imax; $i++)
	{
		$tau = $step * $i;

		# Manually create digits of $tau, for timestamp.
		$itau = int ($tau);
		$t1 = int ($tau * 10 - $itau * 10);
		$t2 = int ($tau * 100 - $itau * 100 - $t1*10);
		print "$i $itau.$t1$t2 \n";

		# system ("ls ptest-$i-*.flo ");
		system ("cat $fpref$i-*.flo | /home/linas/src/fractal/image/flo2mtv |mtvtoppm | pnmtopng > tmp.png");
		system ("convert tmp.png -fill black -draw 'rectangle 5,365,108,393' -pointsize 24 -fill white -gravity SouthWest  -annotate +0+5 ' t=$itau.$t1$t2 ' $fpref$deci$i.png");
	}
}

cvt_range(0,10, "000", "polylog-0.3-");
