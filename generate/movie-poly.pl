#! /usr/bin/env perl
#
# Perl script to generate frames of the polylog movie
#


$fps = 30;
$taumax = 300;  # max tau value but also running time in secs
$nframes = $taumax * $fps; # $taumax seconds at 30 fps
$tau = 0;
$taudelta = $taumax / $nframes;
for ($i=1; $i<$nframes; )
{
	# Manual loop unroll
	$tau = $i * $taudelta;
	print ("Start work on $i tau=$tau\n");
	system ("./polylog-0.3 /home2/linas/tmp/polylog-0.3-$i-$tau 400 400 1 0 0 7.0 $tau &");
	system ("./polylog-0.5 /home2/linas/tmp/polylog-0.5-$i-$tau 400 400 1 0 0 7.0 $tau &");
	system ("./polylog-0.7 /home2/linas/tmp/polylog-0.7-$i-$tau 400 400 1 0 0 7.0 $tau &");

	system ("./polylog-0.3-g1m /home2/linas/tmp/polylog-0.3-g1m-$i-$tau 400 400 1 0 0 7.0 $tau &");
	system ("./polylog-0.5-g1m /home2/linas/tmp/polylog-0.5-g1m-$i-$tau 400 400 1 0 0 7.0 $tau &");
	system ("./polylog-0.7-g1m /home2/linas/tmp/polylog-0.7-g1m-$i-$tau 400 400 1 0 0 7.0 $tau &");
	$i++;

	$tau = $i * $taudelta;
	print ("Start work on $i tau=$tau\n");
	system ("./polylog-0.3 /home2/linas/tmp/polylog-0.3-$i-$tau 400 400 1 0 0 7.0 $tau &");
	system ("./polylog-0.5 /home2/linas/tmp/polylog-0.5-$i-$tau 400 400 1 0 0 7.0 $tau &");
	system ("./polylog-0.7 /home2/linas/tmp/polylog-0.7-$i-$tau 400 400 1 0 0 7.0 $tau &");

	system ("./polylog-0.3-g1m /home2/linas/tmp/polylog-0.3-g1m-$i-$tau 400 400 1 0 0 7.0 $tau &");
	system ("./polylog-0.5-g1m /home2/linas/tmp/polylog-0.5-g1m-$i-$tau 400 400 1 0 0 7.0 $tau &");
	system ("./polylog-0.7-g1m /home2/linas/tmp/polylog-0.7-g1m-$i-$tau 400 400 1 0 0 7.0 $tau");
	$i++;
}