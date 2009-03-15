
-- This file describes SQL tables as appropriate for postgres.

-- Table defining random word-pair (bigram) counts and statistics.
-- The structure of these tables is (*should be*) identical
-- to that of the Pair, etc. tables. However, the source
-- of data used to populate these tables is different: 
-- the words have been randomly created. All word pairs in
-- a sentence are considered.
-- =========================================================

CREATE TABLE RandWordCount (
	count FLOAT
);
INSERT INTO RandWordCount VALUES (0);

CREATE TABLE RandWordPairCount (
	count FLOAT
);
INSERT INTO RandWordPairCount VALUES (0);

CREATE TABLE RandWords (
	word TEXT PRIMARY KEY,
	count FLOAT,
	probability FLOAT
);

CREATE TABLE LeftRandWords (
	left_word TEXT PRIMARY KEY,
	count FLOAT,
	probability FLOAT
);

CREATE TABLE RightRandWords (
	right_word TEXT PRIMARY KEY,
	count FLOAT,
	probability FLOAT
);

CREATE TABLE RandWordPairs (
	left_word TEXT NOT NULL,
	right_word TEXT NOT NULL,
	count FLOAT,
	probability FLOAT,
	mutual_info FLOAT,
 	PRIMARY KEY (left_word, right_word)
);

