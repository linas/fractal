//
// FILE:
// pairhash.C
//
// FUNCTION:
// The lagWordPairTable class associates a unique ID number with 
// a pair of ids.
// The GetWordPairID() method returns the ID number of the pair.
//
// HISTORY:
// January 1997 Linas Vepstas

#include <stdio.h>

#include "top.h"
#include "pairhash.h"

// =====================================================

lagWordPairTable :: lagWordPairTable (void) {
   unused_id = 0;
   num_collisions = 0;
   num_processed = 0;
   num_entries = 0;

   int i = 0;
   for (i=0; i<LAG_PAIR_HASH_TABLE_SIZE; i++) {
      table[i] = 0x0;
   }

   idx = new Helper *[LAG_PAIR_TABLE_SIZE];
   for (i=0; i<LAG_PAIR_TABLE_SIZE; i++) {
      idx[i] = 0x0;
   }

   for (i=0; i<LAG_TOP_TEN; i++) {
      topten[i] = 0x0;
      top_ten_count[i] = 0;
   }
}

// =====================================================

lagWordPairTable :: ~lagWordPairTable () {
   unused_id = 0;
   num_collisions = 0;
   num_processed = 0;
   num_entries = 0;

   int i = 0;
   for (i=0; i<LAG_PAIR_HASH_TABLE_SIZE; i++) {
      if (table[i]) {
         Helper * root = table[i];
         while (root) {
            Helper * nxt = root -> next;
            delete root;
            root = nxt;
         }
      }
      table[i] = 0x0;
   }

   for (i=0; i<LAG_PAIR_TABLE_SIZE; i++) {
      idx[i] = 0x0;
   }

   for (i=0; i<LAG_TOP_TEN; i++) {
      topten[i] = 0x0;
      top_ten_count[i] = 0;
   }
}

// =====================================================

int lagWordPairTable :: GetWordPairID (int first, int second) {

   if (!first) return 0;
   if (!second) return 0;
   num_processed ++;

   unsigned short fi = first;
   unsigned short se = second;
   unsigned short hash = fi + se;
   unsigned int cat = fi << 16 | se;

   Helper * root = table[hash];
   while (root) {
      if (cat == root -> pair) {
         root -> cnt ++;
         UPDATE_TOP_TEN ((root->cnt), (root));
         return (root -> id);
      }
      root = root -> next;
   }

   // if we got to here, we don't yet know the word.
   // Add it to the dictionary
   if (table[hash]) num_collisions ++;

   root = new Helper;
   root -> next = table[hash];
   table[hash] = root;

   root -> pair = cat;

   unused_id ++;
   root -> id = unused_id;
   idx[unused_id] = root;

   root -> cnt = 1;
   num_entries ++;

   return (root->id);
}

// =====================================================

int lagWordPairTable :: GetFirstOfPair (int id) {
   if (!idx[id]) return 0;
   unsigned int retval = idx[id] -> pair;
   retval >>= 16;
   return retval;   
}

// =====================================================

int lagWordPairTable :: GetSecondOfPair (int id) {
   if (!idx[id]) return 0;
   unsigned int retval = idx[id] -> pair;
   retval &= 0xffff;
   return retval;   
}

// =====================================================

int lagWordPairTable :: GetCount (int id) {
   if (!idx[id]) return 0;
   unsigned int retval = idx[id] -> cnt;
   return retval;   
}

// =====================================================

int lagWordPairTable :: GetTopTen (int n) {
   if (LAG_TOP_TEN <= n) return 0;
   if (0 > n) return 0;
   if (!topten[n]) return 0;
   unsigned int retval = topten[n] -> id;
   return retval;   
}

// =====================================================

void lagWordPairTable :: Dump (void) {
   printf ("num entries is %d \n", num_entries);
   printf ("num processed is %d \n", num_processed);
   printf ("num collisions is %d \n", num_collisions);

   // for (int i=0; i<LAG_PAIR_HASH_TABLE_SIZE; i++) {
   //   if (table[i]) printf ("its %d %d %d %d \n", i, 
   //        table[i] -> cnt, table[i]->id, table[i] -> pair);
   // }

}

// ===================== END OF FILE ===================

