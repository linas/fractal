#! /usr/bin/env perl
#
# gen-string.pl
#
# Generate strings of random letters.
#
# Linas Vepstas March 2009

use strict;

my $Nletters = 13;
my @letters = ('e', 't', 'o', 'a', 'i', 'n', 's', 'h', 'r', 'd', 'l', 'u');

sub make_random_word
{
	my $word = "";
	my $i;
	for ($i=0; $i<12; $i++)
	{
		my $rn = rand();
		$rn *= $Nletters;
		if ($rn < 1) { last;}
		$rn -= 1;
		my $l = $letters[$rn];
		
		$word = $word . $l;
	}
	$word;
}

my $w = make_random_word();
print "its $w\n";
