#! /usr/bin/env perl
#
# gen-string.pl
#
# Generate strings of random letters.
#
# Linas Vepstas March 2009

use strict;


my $Nletters = 29;
my @letters = ('e', 'e', 'e', 'e', 
't', 't', 't', 
'o', 'o', 'o',
'a', 'a', 'a',
'i', 'i', 'i',
'n', 'n', 's', 's', 'h', 'h', 'r', 'r', 
'd', 'd', 'l', 'l', 'u');

sub make_random_word
{
	my $word = "";
	my $i;
	for ($i=0; $i<12; $i++)
	{
		my $rn = rand();
		$rn *= $Nletters + 3;
		if (($rn < 3) && ($i > 0)) { last;}
		$rn -= 3;
		my $l = $letters[$rn];
		
		$word = $word . $l;
	}
	$word;
}

my $w = make_random_word();
print "its $w\n";
