//
// FILE:
// wordhash.C
//
// The lagWordtable class associates a unique ID number with a string.
// The GetWordID() method returns the ID number of the string.
//
// January 1997 Linas Vepstas

#include <stdio.h>
#include <string.h>

#include "top.h"
#include "wordhash.h"

// =====================================================

lagWordTable :: lagWordTable (void) {
   unused_id = 0;
   num_collisions = 0;
   num_processed = 0;
   num_entries = 0;
   total_string_bytes = 0;

   int i = 0;
   for (i=0; i<LAG_HASH_TABLE_SIZE; i++) {
      table[i] = 0x0;
   }
   for (i=0; i<LAG_WORD_TABLE_SIZE; i++) {
      idx[i] = 0x0;
   }

   for (i=0; i<LAG_TOP_TEN; i++) {
      topten[i] = 0x0;
      top_ten_count[i] = 0;
   }
}

// =====================================================

lagWordTable :: ~lagWordTable () {
   unused_id = 0;
   num_collisions = 0;
   num_processed = 0;
   num_entries = 0;
   total_string_bytes = 0;

   int i = 0;
   for (i=0; i<LAG_HASH_TABLE_SIZE; i++) {
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

   for (i=0; i<LAG_WORD_TABLE_SIZE; i++) {
      idx[i] = 0x0;
   }

   for (i=0; i<LAG_TOP_TEN; i++) {
      topten[i] = 0x0;
      top_ten_count[i] = 0;
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

   int len = strlen(word) + 1;
   total_string_bytes += len;

   root -> theword = new char[len];
   strcpy (root -> theword, word);

   unused_id ++;
   root -> id = unused_id;
   idx[unused_id] = root;

   root -> cnt = 1;
   num_entries ++;

   return (root->id);
}

// =====================================================

char * lagWordTable :: GetWordFromID (int id) {
   unsigned short uid = id;
   Helper * node = idx[uid];
   if (!node) return 0x0;
   return node -> theword;
}

// =====================================================

int lagWordTable :: GetCount (int id) {
   unsigned short uid = id;
   Helper * node = idx[uid];
   if (!node) return 0;
   return node -> cnt;
}

// =====================================================

void lagWordTable :: Dump (void) {
   // for (int i=0; i<LAG_HASH_TABLE_SIZE; i++) {
   //   if (table[i]) printf ("its %d %d %d %s \n", i, 
   //        table[i] -> cnt, table[i]->id, table[i] -> theword);
   // }
   printf ("Info: lagWordTable::Dump(): \n");
   printf ("num entries is %d \n", num_entries);
   printf ("num processed is %d \n", num_processed);
   printf ("num collisions is %d \n", num_collisions);

   printf ("total string mem usage is %d kbytes \n", total_string_bytes / 1024);
   printf ("helper mem usage is %d kbytes \n", num_entries * sizeof (Helper) / 1024);
   printf ("\n");
}

// ===================== END OF FILE ===================

