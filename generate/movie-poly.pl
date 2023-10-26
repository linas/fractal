#! /usr/bin/env perl
#
# Perl script to generate frames of the polylog movie
#
$SIG{'INT'} = sub { print "bye-bye\n";  die "the end"; };

$fps = 30;  # frames per second
$taumax = 300;  # max tau value but also running time in secs
$nframes = $taumax * $fps; # $taumax seconds at 30 fps
$tau = 0;
$taudelta = $taumax / $nframes;

# Width, in the complex plane. Since it's centered, this
# will go from -3.5 to +3.5 obn the real axis.
$width = 7.0;

$tmpdir="/home2/linas/tmp";

for ($i=1; $i<$nframes; )
{
	# Manual loop unroll
#	$tau = $i * $taudelta;
#	$itau = int($tau);
#	print ("Start work on $i tau=$tau\n");
#	system ("./polylog-0.3 $tmpdir/polylog-0.3-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.5 $tmpdir/polylog-0.5-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.7 $tmpdir/polylog-0.7-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#
#	system ("./polylog-0.3-g1m $tmpdir/polylog-0.3-g1m-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.5-g1m $tmpdir/polylog-0.5-g1m-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.7-g1m $tmpdir/polylog-0.7-g1m-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#
#	system ("./polylog-0.3-g1p $tmpdir/polylog-0.3-g1p-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.5-g1p $tmpdir/polylog-0.5-g1p-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.7-g1p $tmpdir/polylog-0.7-g1p-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#
#	system ("./polylog-0.3-g10z $tmpdir/polylog-0.3-g10z-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.5-g10z $tmpdir/polylog-0.5-g10z-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.7-g10z $tmpdir/polylog-0.7-g10z-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	$i++;

	$tau = $i * $taudelta;
	$itau = int($tau);
	print ("Start work on $i tau.ext=$tau\n");
#	system ("./polylog-0.3 $tmpdir/polylog-0.3-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.5 $tmpdir/polylog-0.5-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.7 $tmpdir/polylog-0.7-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#
#	system ("./polylog-0.3-g1m $tmpdir/polylog-0.3-g1m-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.5-g1m $tmpdir/polylog-0.5-g1m-$i-$itau.ext 400 400 1 0 0 $width $tau &");
#	system ("./polylog-0.7-g1m $tmpdir/polylog-0.7-g1m-$i-$itau.ext 400 400 1 0 0 $width $tau &");

	system ("./polylog-0.3-g1p $tmpdir/polylog-0.3-g1p-$i-$itau.ext 400 400 1 0 0 $width $tau &");
	system ("./polylog-0.5-g1p $tmpdir/polylog-0.5-g1p-$i-$itau.ext 400 400 1 0 0 $width $tau &");
	system ("./polylog-0.7-g1p $tmpdir/polylog-0.7-g1p-$i-$itau.ext 400 400 1 0 0 $width $tau &");

	system ("./polylog-0.3-g10z $tmpdir/polylog-0.3-g10z-$i-$itau.ext 400 400 1 0 0 $width $tau &");
	system ("./polylog-0.5-g10z $tmpdir/polylog-0.5-g10z-$i-$itau.ext 400 400 1 0 0 $width $tau &");
	system ("./polylog-0.7-g10z $tmpdir/polylog-0.7-g10z-$i-$itau.ext 400 400 1 0 0 $width $tau");
	$i++;
}
