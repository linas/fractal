//
// FILE:
// concord.h
//
// FUNCTION:
// The lagConcord class associates builds a concordance
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

#ifndef __LAG_CONCORD_H__
#define __LAG_CONCORD_H__

#include "config.h"

class lagGenericConcordTable :
   public lagGenericWordTable
{
   public:
      lagGenericConcordTable (void);
      ~lagGenericConcordTable ();

      int GetTopTupleContainingWord (int word);

      void Dump (void);

   protected:

      class Concord {
         public:
            Concord * next;
            lagGenericWordTable::Helper * where;
      };
      Concord * concordance [LAG_WORD_TABLE_SIZE];

      void AddID (Helper *);

      int num_concords;
};
   
#endif // __LAG_CONCORD_H__

// ================= END OF FILE ==================


