#! /usr/bin/env perl
# $x = 1/2 + 1/4 + 1/16;

$x = 0;
$x += 1/(1<<(2-1));
$x += 1/(1<<(3-1));
$x += 1/(1<<(5-1));
$x += 1/(1<<(7-1));
$x += 1/(1<<(11-1));
$x += 1/(1<<(13-1));
$x += 1/(1<<(17-1));
$x += 1/(1<<(19-1));
$x += 1/(1<<(23-1));
$x += 1/(1<<(29-1));
$x += 1/(1<<(31-1));
$x += 1/(1<<(37-1));
print "Its $x\n";