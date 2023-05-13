#! /usr/bin/env perl
#
# Show distribution of single words.
#
# Copyright (C) 2009 Linas Vepstas <linasvepstas@gmail.com>

use DBI;
my $dbh = DBI->connect('DBI:Pg:dbname=rexat', 'linas', 'asdf')
   or die "Couldn't connect to database: " . DBI->errstr;

my $mi_table = "RandWords";

# =============================================================
#
my $query = $dbh->prepare('SELECT word, probability FROM ' . $mi_table . ' ORDER BY probability DESC')
	or die "Couldn't prepare statement: " . $dbh->errstr;

local $| = 1;

print "#\n# Single words $mi_table\n";
print "#\n#\n";

$query->execute()
		 or die "Couldn't execute statement: " . $query->errstr;
my $nrows = $query->rows;

my $i = 0;
my $tot = 0;
if ($nrows > 10000) { $nrows = 10000; }
while ($i < $nrows)
{
	my ($word, $prob) = $query->fetchrow_array();
	$tot += $prob;
	print "$i	$prob	$tot	$word\n";
	
	$i ++;
}

