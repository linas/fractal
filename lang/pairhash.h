//
// FILE:
// pairhash.h
//
// FUNCTION:
// The lagWordPairTable class associates a unique ID number with 
// a pair of ints
//
// METHODS:
// The GetWordPairID() method returns the ID number of the pair.
//
// The GetFirstOfPair() method returns the first elt of the pair.
//
// The GetSecondOfPair() method returns the second elt. of the pair.
//
// The GetTopPairContainingWord() returns a pair id that contains the
//    word as the first word.  The returned pairid is the pair 
//    id with the highest count.  
//
// The GetCount() method returns how often the particular word
//    occurs in the text.
//
// The GetTopTen() method returns the n'th most popular pair id.
//
// HISTORY:
// January 1997 Linas Vepstas

#ifndef __LAG_PAIR_TABLE_H__
#define __LAG_PAIR_TABLE_H__

#include "config.h"

class lagWordPairTable {
   public:
      lagWordPairTable (void);
      ~lagWordPairTable ();
      int GetWordPairID (int wordfirst, int wordsecond);
      int GetFirstOfPair (int pairid);
      int GetSecondOfPair (int pairid);
      int GetCount (int pairid);

      int GetTopPairContainingWord (int word);

      int GetTopTen (int n);
      void Dump (void);

   private:
      class Helper {
         public:
            Helper * next;
            unsigned int pair;
            int id;
            int cnt;
      };
      Helper * table [LAG_PAIR_HASH_TABLE_SIZE];
      Helper ** idx;

      Helper * topten [LAG_TOP_TEN];
      int top_ten_count [LAG_TOP_TEN];

      class Concord {
         public:
            Concord * next;
            Helper * where;
      };
      Concord * concordance [LAG_WORD_TABLE_SIZE];

      void AddConcord (int, Helper *);

      int unused_id;
      int num_entries;
      int num_processed;
      int num_collisions;
};
   
#endif // __LAG_PAIR_TABLE_H__

// ================= END OF FILE ==================


