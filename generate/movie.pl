#! /usr/bin/env perl
#
# Perl script to generate frames of a movie
#

$nframes = 5400;
$taumax = 150;
$tau = 0;
$taudelta = $taumax / $nframes;
for ($i=1852; $i<$nframes; $i++)
{
	$tau = $i * $taudelta;
	system ("./polylog /home2/linas/tmp/phase-small-$i-$tau 300 300 3 -0.35 0 3 $tau");
}
