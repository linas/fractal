//
// FILE:
// concord.h
//
// FUNCTION:
// The lagConcord class builds a concordance.
//
// METHODS:
// The GetTopTupleContainingWord() returns a tuple id that contains the
//    word as the first word.  The returned tupleid is the tuple 
//    id with the highest count.  
//
// The ComputeWeights() method computes the strength of each link in the
//    concordance, normalized in the range 0.0 to 1.0
//
// The ResetToStart() method sets to the top of phrase list.
//
// The GetStart() method provides a cursor-based re-entrant interface
// The GetNext() method provides a cursor-based re-entrant interface
//
// The memebr concordance contains a reverse index.  It is index by the
// first cord of each tuple
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
      virtual ~lagGenericConcordTable ();

      int GetTopTupleContainingWord (int word);
      void ComputeWeights (void);

      void ResetToStart (unsigned int phrase);
      float GetNextLinkWeight (void);
      unsigned int GetNextPhrase (void);

      void * GetStart (unsigned int phrase);
      unsigned int GetPhrase (void *);
      void * GetNext (void *);

      void Dump (void);

   protected:

      class Concord {
         public:
            Concord * next;
            lagGenericWordTable::Helper * where;
#ifdef LAG_USE_OVERLOADED_NEW
            void * operator new (size_t s);
            static int memblocks;
            static int memleft;
            static char * mempool;
#endif // LAG_USE_OVERLOADED_NEW
      };
      Concord ** concordance;

      virtual void AddID (Helper *);

      int num_concords;

      Concord * cursor;
      
};
   
#endif // __LAG_CONCORD_H__

// ================= END OF FILE ==================


