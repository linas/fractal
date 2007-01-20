#! /usr/bin/env perl

$step = 1/30;

for ($i=0; $i<10; $i++)
{
	$tau = $step * $i;
	$itau = int ($tau);
	$t1 = int ($tau * 10 - $itau * 10);
	$t2 = int ($tau * 100 - $itau * 100 - $t1*10);
	print "$i $itau.$t1$t2 \n";

	# system ("ls ptest-$i-*.flo ");
	system ("cat poly-mid-$i-*.flo | /home/linas/src/fractal/image/flo2mtv |mtvtoppm | pnmtopng > tmp.png");
	system ("convert tmp.png -fill black -draw 'rectangle 5,365,108,393' -pointsize 24 -fill white -gravity SouthWest  -annotate +0+5 ' t=$itau.$t1$t2 ' mid-000$i.png");
}

for ($i=10; $i<87; $i++)
{
	$tau = $step * $i;
	$itau = int ($tau);
	$t1 = int ($tau * 10 - $itau * 10);
	$t2 = int ($tau * 100 - $itau * 100 - $t1*10);
	print "$i $itau.$t1$t2 \n";

	# system ("ls ptest-$i-*.flo ");
	system ("cat poly-mid-$i-*.flo | /home/linas/src/fractal/image/flo2mtv |mtvtoppm | pnmtopng > tmp.png");
	system ("convert tmp.png -fill black -draw 'rectangle 5,365,108,393' -pointsize 24 -fill white -gravity SouthWest  -annotate +0+5 ' t=$itau.$t1$t2 ' mid-00$i.png");
}

