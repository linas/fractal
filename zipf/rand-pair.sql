
-- This file describes SQL tables as appropriate for postgres.

-- Table defining word-pair (bigram) counts and statistics.
-- The structure of these tables is (*should be*) identical
-- to that of the LemmaPair, etc. tables. However, the source
-- of data used to populate these tables is different: rather
-- than counting word pairs that are linked to one another by
-- link-grammar, these merely count word pairs that happen to
-- occur in the same sentence.  Thus, this table is more 
-- appropriate for use with Mihalcea Markov-chain type algos.
-- =========================================================

CREATE TABLE NearbyCount (
	count FLOAT
);
INSERT INTO NearbyCount VALUES (0);

CREATE TABLE NearbyPairCount (
	count FLOAT
);
INSERT INTO NearbyPairCount VALUES (0);

CREATE TABLE Nearbys (
	lemma TEXT NOT NULL,
	pos TEXT NOT NULL,
	count FLOAT,
	probability FLOAT,
	PRIMARY KEY (lemma, pos)
);

CREATE TABLE LeftNearbys (
	left_lemma TEXT NOT NULL,
	left_pos TEXT NOT NULL,
	count FLOAT,
	probability FLOAT,
	PRIMARY KEY (left_lemma, left_pos)
);

CREATE TABLE RightNearbys (
	right_lemma TEXT NOT NULL,
	right_pos TEXT NOT NULL,
	count FLOAT,
	probability FLOAT,
	PRIMARY KEY (right_lemma, right_pos)
);

CREATE TABLE NearbyPairs (
	left_lemma TEXT NOT NULL,
	left_pos TEXT NOT NULL,
	right_lemma TEXT NOT NULL,
	right_pos TEXT NOT NULL,
	count FLOAT,
	probability FLOAT,
	mutual_info FLOAT,
 	PRIMARY KEY (left_lemma, left_pos, right_lemma, right_pos)
);

