//
// FILE:
// pairhash.h
//
// FUNCTION:
// The lagGenericWordTable class associates a unique ID number with 
// a pair of ints
//
// METHODS:
// The GetGenericWordID() method returns the ID number of the pair.
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

#ifndef LAG_WORD_TUPLE
#define LAG_WORD_TUPLE 2
#endif

class lagGenericWordTable {
   public:
      lagGenericWordTable (void);
      ~lagGenericWordTable ();

#ifdef LAG_TWO_WORD
      int GetID (int wordfirst, int wordsecond);
#endif // LAG_TWO_WORD

#ifdef LAG_THREE_WORD
      int GetID (int, int, int);
#endif // LAG_THREE_WORD

#ifdef LAG_FOUR_WORD
      int GetID (int, int, int, int);
#endif // LAG_FOUR_WORD

#ifdef LAG_FIVE_WORD
      int GetID (int, int, int, int, int);
#endif // LAG_FIVE_WORD

#ifdef LAG_SIX_WORD
      int GetID (int, int, int, int, int, int);
#endif // LAG_SIX_WORD

      int GetElt (int pairid, int elt);
      int GetCount (int pairid);

      int GetTopPairContainingWord (int word);

      int GetTopTen (int n);
      void Dump (void);

   private:
      class Helper {
         public:
            Helper * next;
            unsigned int tuple[LAG_WORD_TUPLE];
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


