#! /usr/bin/env perl
#
# Perl script to generate frames of a movie
#

$nframes = 30;
$taumax = 100;
$tau = 0;
$taudelta = $taumax / $nframes;
for ($i=0; $i<$nframes, $i++)
{
	system ("./polylog /home2/linas/tmp/ptest-$i-$tau 600 600 3 -0.35 0 3 $tau");
	$tau += $taudelta;
}
