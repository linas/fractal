//
// FILE:
// wordhash.h
//
// FUNCTION:
// The lagWordTable class associates a unique ID number with a string.
//
// The GetWordID() method returns the ID number of the string.
//
// The GetWordFromID() method returns the string associated with the ID.
//
// The GetCount() method returns the number of times the indicated
//    word occurs in the text
//
// The Hash() method returns a fairly evenly distributed short
//    suitable for a hash index.
//
// HISTORY:
// January 1997 Linas Vepstas

#ifndef __LAG_WORD_TABLE_H__
#define __LAG_WORD_TABLE_H__

#include "config.h"

class lagWordTable {
   public:
      lagWordTable (void);
      ~lagWordTable ();
      int GetWordID (char *);
      char * GetWordFromID (int);
      int GetCount (int);
      void Dump (void);

   protected:
      int Hash (char *);

   private:
      class Helper {
         public:
            Helper * next;
            char * theword;
            int id;
            int cnt;
      };
      Helper * table [LAG_HASH_TABLE_SIZE];
      Helper * idx [LAG_WORD_TABLE_SIZE];

      Helper * topten [LAG_TOP_TEN];
      int top_ten_count [LAG_TOP_TEN];

      int unused_id;
      int num_entries;
      int num_processed;
      int num_collisions;
};
   
#endif // __LAG_WORD_TABLE_H__

// ================= END OF FILE ==================



