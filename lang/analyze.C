//
// FILE:
// analyze.C
//
// FUNCTION:
// analyze texts
//
// METHODS:
// The Analyze() method breaks a text up into a series of words,
// and adds "phrases" of four words to a collection. It also establishes
// links between phrases.
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
      void Dialogue (FILE *fh);
      void Dump (void);
      void Chain (int);
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

            // light up the phrase neuron
            phrase_table -> AccumStrength (this_phrase_id, 1.0);

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

void lagTextAnalysis :: Dialogue (FILE *fh) {

   xref_table -> ComputeWeights ();
   phrase_table -> ResetAllWeights ();

   // light up only this dialogue
   Analyze (fh);

   int num_phrases = phrase_table -> GetTableSize();

   phrase_table -> FlipAllWeights ();
   for (int i=1; i<= num_phrases; i++) {
      
      float phrase_activation = phrase_table -> GetWeight (i);

      // a negative activation implies that the current neuron is "hot"
      // and should be only a source, not a sink.
      if (0.0 > phrase_activation) {

         phrase_activation = -phrase_activation;
         // activate first layer of neurons
         xref_table -> ResetToStart (i);
   
         float link_strength = xref_table -> GetNextLinkWeight ();
         unsigned int ph = xref_table -> GetNextPhrase();
         while (ph) {
            float not_hot = phrase_table -> GetWeight (ph);

            // light up the neuron only if its not a hot neuron
            if (-1.0e-8 <= not_hot) {
               phrase_table -> AccumStrength (ph, link_strength*phrase_activation);
            }
            link_strength = xref_table -> GetNextLinkWeight ();
            ph = xref_table -> GetNextPhrase();
         }
      }
   }

   // put all non-hot neurons through the activation (squashing)
   // function
   phrase_table -> ActivateAll ();
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

void lagTextAnalysis :: Chain (int top) {

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

main (int argc, char * argv[]) {

   if (3 > argc) {
      printf ("Usage: %s <textfile> <dialogue> \n", argv[0]);
      exit (1);
   }

   printf ("Info: %s: Opening %s for text analysis\n", argv[0], argv[1]); 
   FILE * text_fh = fopen (argv[1], "r");
   if (!text_fh) {
      printf ("Error: %s: no such text file %s \n", argv[0], argv[1]);
      exit (1);
   }

   printf ("Info: %s: Opening %s for dialogue\n", argv[0], argv[2]); 
   FILE * dialogue_fh = fopen (argv[2], "r");
   if (!dialogue_fh) {
      printf ("Error: %s: no such dialogue file %s \n", argv[0], argv[2]);
      exit (1);
   }

   lagTextAnalysis * texan = new lagTextAnalysis;
   texan -> Analyze (text_fh);

   texan -> Dialogue (dialogue_fh);

   texan -> Dump();

   for (int i=0; i<LAG_TOP_TEN; i++) {
      printf ("\n\n =========== sentance %d =========== \n", i);
      texan -> Chain (i);
   }
}

// ========================== END OF FILE ==================
