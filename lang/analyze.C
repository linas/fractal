//
// FILE:
// analyze.C
//
// FUNCTION:
// analyze texts
//
// HISTORY:
// January 1997 Linas Vepstas

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "multihash.h"
#include "wordhash.h"

#define LAG_USED_SIZE 130

class lagTextAnalysis {
   public:
      lagTextAnalysis (void);
      void Analyze (FILE *fh);
      void Dump (void);
      void Construct (int);
   private:
      lagWordTable *wt;
      // lagWordQuadTable *pt;
      lagConcordQuadTable *pt;
      int was_used [LAG_USED_SIZE];
};

// =====================================================

lagTextAnalysis :: lagTextAnalysis (void) {
   wt = new lagWordTable;
   // pt = new lagWordTripleTable;
   // pt = new lagWordQuadTable;
   pt = new lagConcordQuadTable;

   int i=0;
   for (i=0; i<LAG_USED_SIZE; i++) {
      was_used [i] = 0;
   }
}

// =====================================================

void lagTextAnalysis :: Analyze (FILE *fh) {

   char buff[5000];

   int this_word_id = 0;
   int last_word_id = 0;
   int blast_word_id = 0;
   int clast_word_id = 0;
   int dlast_word_id = 0;
   int elast_word_id = 0;
   while (!feof(fh)) {
      fgets (buff, 5000, fh);
      buff[4999] = 0x0;

      int i = 0;
      char * word = buff;
      while (buff[i]) {
      
         if (!isalpha (buff[i])) {
            buff[i] = 0x0;
            // printf (" doing %s \n", word);
            if (this_word_id) {
               elast_word_id = dlast_word_id;
               dlast_word_id = clast_word_id;
               clast_word_id = blast_word_id;
               blast_word_id = last_word_id;
               last_word_id = this_word_id;
            }
            this_word_id = wt -> GetWordID (word); 
            // pt -> GetID (blast_word_id, last_word_id, this_word_id);
            pt -> GetID (clast_word_id, blast_word_id, last_word_id, this_word_id);
            word = &buff[i+1];
         } else {
            buff[i] = tolower (buff[i]);
         }
         i++;
      }
   }

}

// =====================================================

void lagTextAnalysis :: Dump (void) {
   wt -> Dump();
   pt -> Dump();

   int i = 0;
   for (i=0; i<LAG_TOP_TEN; i++) {
      int top = pt -> GetTopTen (i);

      int fi = pt -> GetElt (top, 0);
      int se = pt -> GetElt (top, 1);
      int th = pt -> GetElt (top, 2);
      int fo = pt -> GetElt (top, 3);
      int cnt = pt -> GetCount (top);
      char * wfi = wt -> GetWordFromID (fi);
      char * wse = wt -> GetWordFromID (se);
      char * wth = wt -> GetWordFromID (th);
      char * wfo = wt -> GetWordFromID (fo);
      printf (" its id=%d cnt=%d %s %s %s %s \n", top, cnt, wfi, wse, wth, wfo);
   }
}

// =====================================================

void lagTextAnalysis :: Construct (int top) {

   int i = 0;
   for (i=0; i<LAG_USED_SIZE; i++) {
      was_used [i] = 0;
   }

   // construct string
   int pair = pt -> GetTopTen (top);
   int fi = pt -> GetElt (pair, 0);
   char * wfi = wt -> GetWordFromID (fi);
   printf ("%s ", wfi);

   int j=0;
   while (pair) {

      // avoid infinite loops
      for (i=0; i<LAG_USED_SIZE; i++) {
         if (pair == was_used [i]) {
            return;
         }
         if (!was_used[i]) {
            was_used[i] = pair;
            break;
         }
      }
      
      int se = pt -> GetElt (pair, 1);
      int th = pt -> GetElt (pair, 2);
      int fo = pt -> GetElt (pair, 3);
      if (!se) break;
      if (!th) break;
      if (!fo) break;
      char * wse = wt -> GetWordFromID (se);
      char * wth = wt -> GetWordFromID (th);
      char * wfo = wt -> GetWordFromID (fo);
      printf ("%s %s %s ", wse, wth, wfo);
      j++;
      if (3 == j) {
         j = 0;
         printf ("\n");
      }
      pair = pt -> GetTopTupleContainingWord (fo);
   }
}

// =====================================================

main () {

   FILE * fh = stdin;
   lagTextAnalysis * texan = new lagTextAnalysis;
   texan -> Analyze (fh);
   texan -> Dump();

   for (int i=0; i<LAG_TOP_TEN; i++) {
      printf ("\n\n =========== sentance %d =========== \n", i);
      texan -> Construct (i);
   }


}

// ========================== END OF FILE ==================
