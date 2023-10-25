#! /usr/bin/env perl
#
# Run processes to bulk convert polylog movie data files to PNG files.
#

use feature 'signatures';

# fps == frames per second
# Must match setting for script that initiall generated the frames.
$fps = 30;
$step = 1/ $fps;

$taumax = 300;  # max tau value but also running time in secs

# $fpref == Filename prefix
# $deci == decimal prefix for output.
# $sigma == sigma of sigma+i tau
sub cvt_range($imin, $imax, $deci, $sigma, $fpref) {

	for (my $i=$imin; $i<$imax; $i++)
	{
		$tau = $step * $i;

		# Manually create digits of $tau, for timestamp.
		$itau = int ($tau);
		$t1 = int ($tau * 10 - $itau * 10);
		$t2 = int ($tau * 100 - $itau * 100 - $t1*10);
		print "$i $itau.$t1$t2 \n";

		# system ("ls $fpref$i-*.flo ");
		system ("cat $fpref$i-*.flo | /home/linas/src/fractal/image/flo2mtv |mtvtoppm | pnmtopng > tmp.png");
		system ("convert tmp.png -fill black -draw 'rectangle 5,365,183,393' -pointsize 24 -fill white -gravity SouthWest  -annotate +0+5 ' s = $sigma +i $itau.$t1$t2 ' $fpref$deci$i.png");
	}
}

sub cvt($sigma, $fpref) {
	cvt_range(1, 10, "000", $sigma, $fpref);
}

cvt("0.5", "polylog-0.5-");
