//
// FILE:
// pairhash.C
//
// FUNCTION:
// The lagGenericWordTable class associates a unique ID number with 
// a pair of ids.
// The GetGenericWordID() method returns the ID number of the pair.
//
// HISTORY:
// January 1997 Linas Vepstas

#include <stdio.h>

#include "top.h"
#include "pairhash.h"

// =====================================================

lagGenericWordTable :: lagGenericWordTable (void) {
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

   for (i=0; i<LAG_WORD_TABLE_SIZE; i++) {
      concordance[i] = 0x0;
   }

}

// =====================================================

lagGenericWordTable :: ~lagGenericWordTable () {
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

   for (i=0; i<LAG_WORD_TABLE_SIZE; i++) {
      if (concordance[i]) {
         Concord * root = concordance[i];
         while (root) {
            Concord * nxt = root -> next;
            delete root;
            root = nxt;
         }
      }
      concordance[i] = 0x0;
   }
}

// =====================================================

void lagGenericWordTable :: AddConcord (int first, Helper *where) {
   if (!first) return;
   if (!where) return;

   // first, see if the pair is already listed in the concordance
   // If it is, don't add it again.
   // actually, this check is not needed, since the way this method is
   // invoked guarentees that it will ot already bee there.
   Concord *root = concordance [first];
   while (root) {
      if (where == root->where) return;
      root = root -> next;
   }
      
   root = new Concord;
   root -> where = where;
   root -> next = concordance[first];
   concordance[first] = root; 
}

#if (defined LAG_TWO_WORD) | (defined LAG_THREE_WORD) | \
    (defined LAG_FOUR_WORD) | (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)

#ifdef LAG_TWO_WORD
int lagGenericWordTable :: GetID (int first, int second) {
#endif // LAG_TWO_WORD

#ifdef LAG_THREE_WORD
int lagGenericWordTable :: GetID (int first, int second, int third) {
#endif // LAG_THREE_WORD

#ifdef LAG_FOUR_WORD
int lagGenericWordTable :: GetID (int first, int second, int third,
                                  int fourth) {
#endif // LAG_FOUR_WORD

#ifdef LAG_FIVE_WORD
int lagGenericWordTable :: GetID (int first, int second, int third,
                                  int fourth, int fifth) {
#endif // LAG_FIVE_WORD

#ifdef LAG_SIX_WORD
int lagGenericWordTable :: GetID (int first, int second, int third,
                                  int fourth, int fifth, int sixth) {
#endif // LAG_SIX_WORD

   if (!first) return 0;
   if (!second) return 0;
#if (defined LAG_THREE_WORD) | (defined LAG_FOUR_WORD) | (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
   if (!third) return 0;
#endif
#if (defined LAG_FOUR_WORD) | (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
   if (!fourth) return 0;
#endif
#if (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
   if (!fifth) return 0;
#endif
#if (defined LAG_SIX_WORD)
   if (!sixth) return 0;
#endif


   unsigned short hash = (unsigned short) first;
   hash += (unsigned short) second;
#if (defined LAG_THREE_WORD) | (defined LAG_FOUR_WORD) | (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
   hash += (unsigned short) third;
#endif
#if (defined LAG_FOUR_WORD) | (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
   hash += (unsigned short) fourth;
#endif
#if (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
   hash += (unsigned short) fifth;
#endif
#if (defined LAG_SIX_WORD)
   hash += (unsigned short) sixth;
#endif

   num_processed ++;
   Helper * root = table[hash];
   while (root) {
      if ((first == root -> tuple[0])
         && (second == root -> tuple[1])
#if (defined LAG_THREE_WORD) | (defined LAG_FOUR_WORD) | (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
         && (third == root -> tuple[2])
#endif
#if (defined LAG_FOUR_WORD) | (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
         && (fourth == root -> tuple[3])
#endif
#if (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
         && (fifth == root -> tuple[4])
#endif
#if (defined LAG_SIX_WORD)
         && (sixth == root -> tuple[5])
#endif
         ) {
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

   root -> tuple[0] = first;
   root -> tuple[1] = second;
#if (defined LAG_THREE_WORD) | (defined LAG_FOUR_WORD) | (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
   root -> tuple[2] = third;
#endif
#if (defined LAG_FOUR_WORD) | (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
   root -> tuple[3] = fourth;
#endif
#if (defined LAG_FIVE_WORD) | (defined LAG_SIX_WORD)
   root -> tuple[4] = fifth;
#endif
#if (defined LAG_SIX_WORD)
   root -> tuple[5] = sixth;
#endif

   unused_id ++;
   root -> id = unused_id;
   idx[unused_id] = root;

   root -> cnt = 1;
   num_entries ++;

   AddConcord (first, root);

   return (root->id);
}

#endif 

// =====================================================

int lagGenericWordTable :: GetElt (int id, int elt) {
   if (!idx[id]) return 0;
   unsigned int retval = idx[id] -> tuple[elt];
   return retval;   
}

// =====================================================

int lagGenericWordTable :: GetCount (int id) {
   if (!idx[id]) return 0;
   unsigned int retval = idx[id] -> cnt;
   return retval;   
}

// =====================================================

int lagGenericWordTable :: GetTopTen (int n) {
   if (LAG_TOP_TEN <= n) return 0;
   if (0 > n) return 0;
   if (!topten[n]) return 0;
   unsigned int retval = topten[n] -> id;
   return retval;   
}

// =====================================================

int lagGenericWordTable :: GetTopPairContainingWord (int word) {
   if (0 > word) return 0;
   if (LAG_WORD_TABLE_SIZE <= word) return 0;

   Concord * root = concordance[word];
   if (!root) return 0;

   // search for the pair with the highest count
   Helper * top = 0x0;
   int topcnt = 0;
   while (root) {
      if (topcnt < root -> where -> cnt) {
         topcnt = root->where->cnt;
         top = root -> where;
      }
      root = root -> next;
   }

   return (top -> id);
}

// =====================================================

void lagGenericWordTable :: Dump (void) {
   printf ("num entries is %d \n", num_entries);
   printf ("num processed is %d \n", num_processed);
   printf ("num collisions is %d \n", num_collisions);

   // for (int i=0; i<LAG_PAIR_HASH_TABLE_SIZE; i++) {
   //   if (table[i]) printf ("its %d %d %d %d \n", i, 
   //        table[i] -> cnt, table[i]->id, table[i] -> pair);
   // }

}

// ===================== END OF FILE ===================

