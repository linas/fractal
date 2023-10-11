#! /usr/bin/perl

$l = log (2);

$x = exp (-2*$l);
$x += exp (-3*$l);
$x += exp (-5*$l);
$x += exp (-7*$l);
$x += exp (-11*$l);
$x += exp (-13*$l);
$x += exp (-17*$l);
$x += exp (-19*$l);
$x += exp (-23*$l);
$x += exp (-29*$l);
$x += exp (-31*$l);
$x += exp (-37*$l);
$x += exp (-41*$l);
$x += exp (-43*$l);
$x += exp (-47*$l);
# $x += exp (-53*$l);
# $x += exp (-59*$l);

print "$x\n";
