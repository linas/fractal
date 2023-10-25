#! /usr/bin/env perl
#
# movie-poly-montage.pl
# Bulk montage multiple polylog movie PNG frames.
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

# $deci == decimal prefix for output.
sub monte_range($imin, $imax, $deci) {

	$ap = "polylog-0.3-";
	$bp = "polylog-0.5-";
	$cp = "polylog-0.7-";
	$dp = "polylog-0.3-g1m-";
	$ep = "polylog-0.5-g1m-";
	$fp = "polylog-0.7-g1m-";

	for (my $i=$imin; $i<$imax; $i++)
	{
		$tau = $step * $i;
		$itau = int ($tau);

		$a = "$ap$i-$itau.png";
		$b = "$bp$i-$itau.png";
		$c = "$cp$i-$itau.png";
		$d = "$dp$i-$itau.png";
		$e = "$ep$i-$itau.png";
		$f = "$fp$i-$itau.png";

		$fout = "polylog-montage-$deci$i.png";

		printf "yoooo $a\n";

		# If we had the full filename, we could have said
		# if (-f $filename) but we don't have that.
		# glob() does the same thing, but with globby filenames.
		# Its a tad slow, but whatever.
		# if (glob("$ap$i-*.flo"))
		if (-f $a) {

			system ("ls $a");
			# system ("montage $a $b $c $d $e $f -tile 3x2 -geometry +0+0 $fout");
		}
	}
}

sub monte() {
	monte_range(1, 10, "000");
	# monte_range(10, 100, "00");
	# monte_range(100, 1000, "0");
	# monte_range(1000, 10000, "");
}

monte();
