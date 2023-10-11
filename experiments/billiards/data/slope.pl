#! /usr/bin/perl

$f = 0;

while (<>) {
	if (/#/) { next; }
	($r, $mfp) = split /\t/;

	$lmfp = log ($mfp);
	$lr = log ($r);

   if ($f == 0) {
		$x1 = $lr;
		$y1 = $lmfp;
		$f = 1;
		next;
	}
	$slope = ($lmfp-$y1) / ($lr - $x1);
	print "its $r $slope\n";
}
