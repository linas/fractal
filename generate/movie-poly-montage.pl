#! /usr/bin/env perl
#
# movie-poly-montage.pl
# Bulk montage multiple polylog movie PNG frames.
#
# After this is done, run this:
# ffmpeg -i polylog-montage-%04d.png -framerate 30 polylog-montage.mp4
#
# Paste images side by side (from Imagemagick):
# montage a.png b.png c.png -tile 3x1 -geometry +0+0 out.png
# montage a.png b.png c.png d.png -tile 2x2 -geometry +0+0 out.png

use feature 'signatures';

$SIG{'INT'} = sub { print "bye-bye\n";  die "the end"; };

# $deci == decimal prefix for output.
sub monte_range($imin, $imax, $deci) {

	$gp = "polylog-0.3-g1p-";
	$hp = "polylog-0.5-g1p-";
	$ip = "polylog-0.7-g1p-";

	$ap = "polylog-0.3-";
	$bp = "polylog-0.5-";
	$cp = "polylog-0.7-";
	$dp = "polylog-0.3-g1m-";
	$ep = "polylog-0.5-g1m-";
	$fp = "polylog-0.7-g1m-";

	$jp = "polylog-0.3-g10z-";
	$kp = "polylog-0.5-g10z-";
	$lp = "polylog-0.7-g10z-";

	for (my $i=$imin; $i<$imax; $i++)
	{
		$a = "$ap$deci$i.png";
		$b = "$bp$deci$i.png";
		$c = "$cp$deci$i.png";
		$d = "$dp$deci$i.png";
		$e = "$ep$deci$i.png";
		$f = "$fp$deci$i.png";
		$g = "$gp$deci$i.png";
		$h = "$hp$deci$i.png";
		$i = "$ip$deci$i.png";
		$j = "$jp$deci$i.png";
		$k = "$kp$deci$i.png";
		$l = "$lp$deci$i.png";

		$fout = "polylog-montage-$deci$i.png";

		# If we had the full filename, we could have said
		# if (-f $filename) but we don't have that.
		# glob() does the same thing, but with globby filenames.
		# Its a tad slow, but whatever.
		# if (glob("$ap$i-*.flo"))
		if (-f $a) {

			system ("ls $a");
			# system ("montage $a $b $c $d $e $f -tile 3x2 -geometry +0+0 $fout");
			system ("montage $g $h $i $a $b $c $d $e $f $j $k $l -tile 3x4 -geometry +0+0 $fout");
		}
	}
}

sub monte() {
	monte_range(1, 10, "000");
	monte_range(10, 100, "00");
	monte_range(100, 1000, "0");
	monte_range(1000, 10000, "");
}

monte();
