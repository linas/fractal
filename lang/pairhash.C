//
// FILE:
// pairhash.C
//
// The lagWordPairTable class associates a unique ID number with 
// a pair of ids.
// The GetWordPairID() method returns the ID number of the pair.
//
// January 1997 Linas Vepstas

#include <stdio.h>

#include "pairhash.h"

// =====================================================

lagWordPairTable :: lagWordPairTable (void) {
   unused_id = 0;
   num_collisions = 0;
   num_processed = 0;
   num_entries = 0;
   for (int i=0; i<LAG_PAIR_HASH_TABLE_SIZE; i++) {
      table[i] = 0x0;
   }
}

// =====================================================

lagWordPairTable :: ~lagWordPairTable () {
   unused_id = 0;
   num_collisions = 0;
   num_processed = 0;
   num_entries = 0;
   for (int i=0; i<LAG_PAIR_HASH_TABLE_SIZE; i++) {
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
}

// =====================================================

int lagWordPairTable :: GetWordPairID (int first, int second) {

   if (!first) return 0;
   if (!second) return 0;
   num_processed ++;

   unsigned short fi = first;
   unsigned short se = second;
   unsigned short hash = first + second;

   Helper * root = table[hash];
   while (root) {
      if (hash == root -> pair) {
         root -> cnt ++;
         return (root -> id);
      }
      root = root -> next;
   }

   // if we got to here, we don't yet know the word.
   // Add it to the dictionary
   if (table[hash]) num_collisions ++;
   unused_id ++;

   root = new Helper;
   root -> next = table[hash];
   root -> pair = fi << 16 | se;
   root -> id = unused_id;
   root -> cnt = 1;
   table[hash] = root;
   num_entries ++;

   return (root->id);
}

// =====================================================

void lagWordPairTable :: Dump (void) {
   for (int i=0; i<LAG_PAIR_HASH_TABLE_SIZE; i++) {
      if (table[i]) printf ("its %d %d %d %d \n", i, 
           table[i] -> cnt, table[i]->id, table[i] -> pair);
   }
   printf ("num entries is %d \n", num_entries);
   printf ("num processed is %d \n", num_processed);
   printf ("num collisions is %d \n", num_collisions);
}

// ===================== END OF FILE ===================

