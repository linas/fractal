//
// FILE:
// concord.C
//
// FUNCTION:
// The lagGenericConcordTable class associates a unique ID number with 
// a pair of ids.
// The GetGenericConcordID() method returns the ID number of the pair.
//
// HISTORY:
// January 1997 Linas Vepstas

#include <stdio.h>

#include "top.h"
#include "concord.h"

// =====================================================

lagGenericConcordTable :: lagGenericConcordTable (void) {
   num_concords = 0;

   concordance = new Concord* [LAG_PAIR_TABLE_SIZE];
   int i = 0;
   for (i=0; i<LAG_PAIR_TABLE_SIZE; i++) {
      concordance[i] = 0x0;
   }

}

// =====================================================

lagGenericConcordTable :: ~lagGenericConcordTable () {
   num_concords = 0;

   int i = 0;
   for (i=0; i<LAG_PAIR_TABLE_SIZE; i++) {
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

// Say we have the structure    ph1 ------- link ------- ph2
// then (link -> tuple[0] == ph1) is true
// then (link -> tuple[1] == ph2) is true
// then (concordance [ph1] == link) is true
// 
// in fact, cooncordance[ph1] has the property that it points at all
// links that have the exact same ph1 in them.

void lagGenericConcordTable :: AddID (Helper *where) {
   if (!where) return;

   // first, see if the pair is already listed in the concordance
   // If it is, don't add it again.
   // Actually, this check is not needed, since the way this method is
   // invoked guarentees that it will already be there.
   unsigned int first = where -> tuple[0];

   Concord *root = concordance [first];
   while (root) {
      if (where == root->where) return;
      root = root -> next;
   }
      
   root = new Concord;
   root -> where = where;
   root -> next = concordance[first];
   concordance[first] = root; 

   num_concords ++;
}

// =====================================================

void lagGenericConcordTable :: ComputeWeights (void) 
{

   // every node chained to the same slot in the concordance table
   // shares the trait that all of its first words are alike.

   for (int i=0; i<LAG_PAIR_TABLE_SIZE; i++) {
      int total_count = 0;
      Concord * root = concordance [i];
      while (root) {
         total_count += root -> where -> cnt;
         root = root -> next;
      }
      if (total_count) {
         float inv = 1.0 / ((float) total_count);
         Concord * root = concordance [i];
         while (root) {
            root -> where -> activation = inv * (float) root -> where -> cnt;
            root = root -> next;
         }
      }
   }
}

// =====================================================

void lagGenericConcordTable :: ResetToStart (unsigned int phrase) {
   if (!phrase) return;
   cursor = concordance [phrase];
}

// =====================================================

float lagGenericConcordTable :: GetNextLinkWeight (void) {
   if (!cursor) return 0.0;
   return (cursor -> where -> activation);
}

// =====================================================

unsigned int lagGenericConcordTable :: GetNextPhrase (void) {
   if (!cursor) return 0;
   unsigned int ph = cursor -> where -> tuple[1];
   cursor = cursor -> next;
   return ph;
}

// =====================================================

void * lagGenericConcordTable :: GetStart (unsigned int phrase) {
   if (!phrase) return 0x0;
   return ((void *) concordance [phrase]);
}

// =====================================================

unsigned int lagGenericConcordTable :: GetPhrase (void *curse) {
   if (!curse) return 0;
   Concord *c = (Concord *) curse;
   return (c -> where -> tuple[1]);
}

// =====================================================

void *lagGenericConcordTable :: GetNext (void *curse) {
   if (!curse) return 0x0;
   Concord *c = (Concord *) curse;
   return ((void *) c -> next);
}

// =====================================================

int lagGenericConcordTable :: GetTopTupleContainingWord (int word) {
   if (0 > word) return 0;
   if (LAG_PAIR_TABLE_SIZE <= word) return 0;

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

void lagGenericConcordTable :: Dump (void) {

   printf ("Info: lagGenericConcordTable :: Dump(): \n");
   printf ("num concord is %d \n", num_concords);
   printf ("memusage is %d kbytes \n", num_concords * sizeof (Concord) / 1024);
   printf ("sizeof concord is %d bytes \n", sizeof (Concord));
   printf ("\n");

#ifdef LAG_USE_OVERLOADED_NEW
   printf (" num blocks alloced= %d \n", Concord ::memblocks);
   printf (" block mem left = %d \n", Concord ::memleft);
#endif // LAG_USE_OVERLOADED_NEW

   lagGenericWordTable :: Dump();
}

// =====================================================
#ifdef LAG_USE_OVERLOADED_NEW

int lagGenericConcordTable :: Concord ::  memleft = 0;
int lagGenericConcordTable :: Concord ::  memblocks = 0;
char * lagGenericConcordTable :: Concord ::  mempool = 0x0;

void * lagGenericConcordTable :: Concord ::  operator new (size_t s) {

   if (!mempool || memleft < s) {
      mempool = new char [LAG_CHUNK_SIZE];
      memleft = LAG_CHUNK_SIZE;
      memblocks ++;
   }

   memleft -= s;
   void * obj = mempool;
   mempool += s;
   return obj;
}
#endif // LAG_USE_OVERLOADED_NEW

// ===================== END OF FILE ===================

