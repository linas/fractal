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

#define LAG_USED_SIZE 50

class lagTextAnalysis {
   public:
      lagTextAnalysis (void);
      void Analyze (FILE *fh);
      void Dump (void);
      void Construct (int);
   protected:
      void PrintPhrase (int);

   private:
      lagWordTable *word_table;
      lagWordQuadTable *phrase_table;
      lagConcordPairTable *xref_table;
      int was_used [LAG_USED_SIZE];
};

// =====================================================

lagTextAnalysis :: lagTextAnalysis (void) {
   word_table = new lagWordTable;
   // phrase_table = new lagWordTripleTable;
   phrase_table = new lagWordQuadTable;
   xref_table = new lagConcordPairTable;

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

   int this_phrase_id = 0;
   int last_phrase_id = 0;
   int blast_phrase_id = 0;
   int clast_phrase_id = 0;
   int dlast_phrase_id = 0;
   int elast_phrase_id = 0;

   int lineno = 0;
   while (!feof(fh)) {
      fgets (buff, 5000, fh);
      buff[4999] = 0x0;

      lineno ++;
      if (0 == lineno%1000) { printf ("."); fflush (stdout); }

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
            this_word_id = word_table -> GetWordID (word); 

            if (this_phrase_id) {
               elast_phrase_id = dlast_phrase_id;
               dlast_phrase_id = clast_phrase_id;
               clast_phrase_id = blast_phrase_id;
               blast_phrase_id = last_phrase_id;
               last_phrase_id = this_phrase_id;
            }

            this_phrase_id = phrase_table -> GetID (clast_word_id, blast_word_id, last_word_id, this_word_id);
            xref_table -> GetID (dlast_phrase_id, this_phrase_id);
            word = &buff[i+1];
         } else {
            buff[i] = tolower (buff[i]);
         }
         i++;
      }
   }

   printf ("\n");
}

// =====================================================

void lagTextAnalysis :: PrintPhrase (int phrase_id) {
   if (!phrase_id) return;

   for (int j=0; j<4; j++) {
      int word_id = phrase_table -> GetElt (phrase_id, j);
      char * word = word_table -> GetWordFromID (word_id);
      printf (" %s", word);
   }
}

// =====================================================

void lagTextAnalysis :: Dump (void) {
   word_table -> Dump();
   phrase_table -> Dump();
   xref_table -> Dump();

   printf ("\n");
   printf ("Info: lagTextAnalysis::Dump(): \n");
   printf ("top %d pharses: \n", LAG_TOP_TEN);

   int i = 0;
   for (i=0; i<LAG_TOP_TEN; i++) {
      int top = phrase_table -> GetTopTen (i);

      int cnt = phrase_table -> GetCount (top);
      printf ("top phrase %d cnt=%d ", top, cnt);

      PrintPhrase (top);
      printf ("\n");
   }

   printf ("\n");
   printf ("Info: lagTextAnalysis::Dump(): \n");
   printf ("top %d xrefs: \n", LAG_TOP_TEN);

   for (i=0; i<LAG_TOP_TEN; i++) {
      int top = xref_table -> GetTopTen (i);

      int fi = xref_table -> GetElt (top, 0);
      int se = xref_table -> GetElt (top, 1);
      int cnt = xref_table -> GetCount (top);
      printf ("top xref %d cnt=%d \n", top, cnt);

      PrintPhrase (fi);
      PrintPhrase (se);
      printf ("\n");
   }
}

// =====================================================

void lagTextAnalysis :: Construct (int top) {

   int i = 0;
   for (i=0; i<LAG_USED_SIZE; i++) {
      was_used [i] = 0;
   }

   // construct string
   int phrase = phrase_table -> GetTopTen (top);
   int xref = 0;
   int j=0;
   while (phrase) {

      // avoid infinite loops
      for (i=0; i<LAG_USED_SIZE; i++) {
         if (phrase == was_used [i]) {
            return;
         }
         if (!was_used[i]) {
            was_used[i] = phrase;
            break;
         }
         if (LAG_USED_SIZE-1 == i) return;
      }
      
      PrintPhrase (phrase);
      j++;
      if (3 == j) {
         j = 0;
         printf ("\n");
      }
 
      xref = xref_table -> GetTopTupleContainingWord (phrase);
      if (!xref) break;
      phrase = xref_table -> GetElt (xref, 1); 
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
