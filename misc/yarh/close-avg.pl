#! /usr/bin/perl
#
# Take a gaussian average of the zero locations as a function of N.
#
# Linas Vepstas December 2010
#

$sigma = 2;
$width = 6 * $sigma;

for ($i=-$width; $i<$width; $i++)
{
	$x = $i/$sigma;
	$gauss[$i] = exp(-$x*$x);
	print "Gauss $i $gauss[$i]\n";
}

while (<>)
{

}
