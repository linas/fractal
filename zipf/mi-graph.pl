#! /usr/bin/env perl
#
# Perform a bin-count to obtain a distribution of mutal information.
#
# Copyright (C) 2009 Linas Vepstas <linasvepstas@gmail.com>

use DBI;
my $dbh = DBI->connect('DBI:Pg:dbname=lexat', 'linas', 'asdf')
   or die "Couldn't connect to database: " . DBI->errstr;

my $mi_table = "NearbyPairs";
my $start = 90;
my $end = -60;
my $delta = 0.1;

my $mi_table = "Pairs";
my $start = 60;
my $end = -60;
my $delta = 0.5;

# =============================================================
#
my $query = $dbh->prepare('SELECT probability FROM ' . $mi_table . 
	' WHERE mutual_info > ? AND mutual_info <= ? ')
	or die "Couldn't prepare statement: " . $dbh->errstr;

local $| = 1;

print "#\n# Bin count of mutual information for $mi_table\n";
print "#\n#\n";

my $bin;
my $i = 0;
my $tot = 0;
my $wtot = 0;
for ($bin = $start; $bin > $end; $bin -= $delta)
{
	$query->execute($bin, $bin+$delta)
		 or die "Couldn't execute statement: " . $query->errstr;

	my $nr = $query->rows;

	my $wsum = 0;
	for (my $j=0; $j<$nr; $j++)
	{
		my ($prob) = $query->fetchrow_array();
		$wsum += $prob;
	}
	$nr /= $delta;
	$wsum /= $delta;

	print "$i	$bin	$nr	$wsum\n";
	
	$tot += $nr;
	$wtot += $wsum;
	$i ++;
}

$tot *= $delta;
$wtot *= $delta;
print "#\n# total number of entries: $tot  weighted: $wtot\n";
