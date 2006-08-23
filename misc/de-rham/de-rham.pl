#! /usr/bin/perl
#
# script for autolabeling gnuplot graphs
# 

while (<>) {
	if (/d=/) {
		split; 
		$d=$_[1];
		next;
	}
	if (/e=/) {
		split; 
		$e=$_[1];
		next;
	}
	if (/f=/) {
		split; 
		$f=$_[1];
		next;
	}
	if (/g=/) {
		split; 
		$g=$_[1];
		break;
	}
}

print "set term png small\n";
print "set out 'curve_$d_$e_$f_$g.png'\n";
print "set data style lines\n";
print "set key right\n";
print "set size square\n";

print "set title \"General de Rham curve for $d $e $f $g\"\n";
print "set xlabel \"u\"\n";
print "set ylabel \"v\"\n";

print "plot \"curve.dat\" using 3:4 title \"\"\n";

