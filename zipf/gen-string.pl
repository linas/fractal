#! /usr/bin/env perl
#
# gen-string.pl
#
# Generate strings of random letters.
#
# Linas Vepstas March 2009

use strict;
use warnings;


my $Nletters = 29;
my @letters = ('e', 'e', 'e', 'e', 
't', 't', 't', 
'o', 'o', 'o',
'a', 'a', 'a',
'i', 'i', 'i',
'n', 'n', 's', 's', 'h', 'h', 'r', 'r', 
'd', 'd', 'l', 'l', 'u');

my $num_sentences = 0;
my $num_words = 0;
my $num_pairs = 0;

my %word_freq = ();
my %left_word_freq = ();
my %right_word_freq = ();
my %word_pair_freq = ();

my @word_id = ();  # word to index map


sub make_random_word
{
	my $word = "";
	my $i;
	for ($i=0; $i<12; $i++)
	{
		my $rn = rand();
		$rn *= $Nletters + 5;
		if (($rn < 5) && ($i > 0)) { last;}
		$rn -= 5;
		my $l = $letters[$rn];
		
		$word = $word . $l;
	}
	$word;
}

while ($num_sentences < 53123)
{
	my $i;

	# put 13 words into a sentence.
	@word_id = ();
	for ($i=0; $i<13; $i++)
	{
		my $word = make_random_word();
		$word_freq{$word} += 1;
		$num_words ++;
		
		my $right_word = $word;
		my $li = 0;
		for ($li=0; $li<$i; $li++)
		{
			my $left_word = $word_id[$li];
			$left_word_freq{$left_word} += 1;
			$right_word_freq{$right_word} += 1;

			# Count word pairs.
			my $pair = $left_word . "+\@\@\@+" . $right_word;
			$word_pair_freq{$pair} += 1;
			$num_pairs ++;
		}

		$word_id[$i] = $word;
	}

	$num_sentences ++;
	if ($num_sentences %1000 == 0)
	{
		print "done $num_sentences sentences\n";
	}
}

print "Final sentence count = $num_sentences\n";
print "Final num word = $num_words\n";
print "Final num pairs = $num_pairs\n";


# ================================================================
#
# Now, store the results in an SQL database.

use lib '/home/linas/src/novamente/src/cerego/lexical-attr/src/count';

use DBI;
use StoreCount;

my $dbh = DBI->connect('DBI:Pg:dbname=rexat', 'linas', 'asdf')
	or die "Couldn't connect to database: " . DBI->errstr;

StoreCount::set_table_names(
	"RandWordCount", "RandWordPairCount", "RandWords", 
	"RightRandWords", "LeftRandWords", "RandWordPairs");

# store lemma counts
print "Start storing ...\n";
StoreCount::store_counts($dbh, $num_words, \%word_freq,
	 $num_pairs, \%left_word_freq, \%right_word_freq, \%word_pair_freq);

print "... done storing\n";
