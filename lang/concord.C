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

