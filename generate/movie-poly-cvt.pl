#! /usr/bin/env perl
#
# Run processes to bulk convert polylog movie data files to PNG files.
#
# After this is done, run this:
# ffmpeg -i polylog-0.3-%04d.png -framerate 30 polylog-0.3.mp4
#
# Paste images side by side (from Imagemagick):
# montage a.png b.png c.png -tile 3x1 -geometry +0+0 out.png
# montage a.png b.png c.png d.png -tile 2x2 -geometry +0+0 out.png


use feature 'signatures';

# fps == frames per second
# Must match setting for script that initiall generated the frames.
$fps = 30;
$step = 1/ $fps;

$taumax = 300;  # max tau value but also running time in secs
$nframes = $taumax * $fps; # $taumax seconds at 30 fps

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
		# print "$i $itau.$t1$t2 \n";

		# If we had the full filename, we could have said
		# if (-f $filename) but we don't have that.
		# glob() does the same thing, but with globby filenames.
		# Its a tad slow, but whatever.
		# if (glob("$fpref$i-*.flo"))
		if (-f "$fpref$i-$itau.flo") {
			# print "$i $itau.$t1$t2 \n";

			system ("ls $fpref$i-*.flo ");
			system ("cat $fpref$i-$itau.flo | /home/linas/src/fractal/image/flo2mtv |mtvtoppm | pnmtopng > tmp.png");
			system ("convert tmp.png -fill black -draw 'rectangle 5,365,183,393' -pointsize 24 -fill white -gravity SouthWest  -annotate +0+5 ' s = $sigma +i $itau.$t1$t2 ' $fpref$deci$i.png");
		}
	}
}

sub cvt($sigma, $fpref) {
	cvt_range(1, 10, "000", $sigma, $fpref);
	cvt_range(10, 100, "00", $sigma, $fpref);
	cvt_range(100, 1000, "0", $sigma, $fpref);
	cvt_range(1000, 10000, "", $sigma, $fpref);
}

cvt("0.3", "polylog-0.3-");
# cvt("0.5", "polylog-0.5-");
# cvt("0.7", "polylog-0.7-");
# cvt("0.3", "polylog-0.3-g1m-");
# cvt("0.5", "polylog-0.5-g1m-");
# cvt("0.7", "polylog-0.7-g1m-");
