//
// wordhash.C
//
// The lagWordtable class associates a unique ID number with a string.
// The GetWordID() method returns the ID number of the string.
//
// January 1997 Linas Vepstas

#include <stdio.h>

#include "wordhash.h"

// =====================================================

lagWordTable :: lagWordTable (void) {
   unused_id = 0;
   num_collisions = 0;
   num_processed = 0;
   num_entries = 0;
   for (int i=0; i<LAG_WORD_TABLE_SIZE; i++) {
      table[i] = 0x0;
   }
}

// =====================================================

lagWordTable :: ~lagWordTable () {
   unused_id = 0;
   num_collisions = 0;
   num_processed = 0;
   num_entries = 0;
   for (int i=0; i<LAG_WORD_TABLE_SIZE; i++) {
      if (table[i]) {
         Helper * root = table[i];
         while (root) {
            Helper * nxt = root -> next;
            delete [] root -> theword;
            delete root;
            root = nxt;
         }
      }
      table[i] = 0x0;
   }
}

// =====================================================

int lagWordTable :: Hash (char * word) {

   if (!word) return 0;
   
   int i=0;
   short shift = 0;
   short hash = 0;
   while (word[i]) {
      hash += word[i] << shift; 
      shift += 4;
      if (13 < shift) shift = 0;
      i++;      
   }
   return hash;
}

// =====================================================

int lagWordTable :: GetWordID (char * word) {

   if (!word) return 0;
   if (!word[0]) return 0;
   num_processed ++;

   unsigned short hash = Hash (word);

   Helper * root = table[hash];
   while (root) {
      int match = !strcmp (root -> theword, word);
      if (match) {
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
   root -> theword = new char[strlen(word) +1];
   root -> id = unused_id;
   strcpy (root -> theword, word);
   table[hash] = root;
   num_entries ++;

   return (root->id);
}

// =====================================================

void lagWordTable :: Dump (void) {
   for (int i=0; i<LAG_WORD_TABLE_SIZE; i++) {
      if (table[i]) printf ("its %d %d %d %s \n", i, 
           table[i] -> cnt, table[i]->id, table[i] -> theword);
   }
   printf ("num entries is %d \n", num_entries);
   printf ("num processed is %d \n", num_processed);
   printf ("num collisions is %d \n", num_collisions);
}

// ===================== END OF FILE ===================

