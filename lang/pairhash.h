//
// FILE:
// pairhash.h
//
// The lagWordPairTable class associates a unique ID number with 
// a pair of ints
// The GetWordPairID() method returns the ID number of the pair.
//
// January 1997 Linas Vepstas

#ifndef __LAG_PAIR_TABLE_H__
#define __LAG_PAIR_TABLE_H__

#define LAG_PAIR_HASH_TABLE_SIZE 65536

class lagWordPairTable {
   public:
      lagWordPairTable (void);
      ~lagWordPairTable ();
      int GetWordPairID (int, int);
      void Dump (void);

   private:
      class Helper {
         public:
            Helper * next;
            int pair;
            int id;
            int cnt;
      };
      Helper * table [LAG_PAIR_HASH_TABLE_SIZE];

      int unused_id;
      int num_entries;
      int num_processed;
      int num_collisions;
};
   
#endif // __LAG_PAIR_TABLE_H__

// ================= END OF FILE ==================

