//
// FILE:
// analyze.C
//
// FUNCTION:
// analyze texts
//
// METHODS:
// The Analyze() method breaks a text up into a series of words,
//    and adds "phrases" of four words to a collection. It also establishes
//    links between phrases.
//
// The Chain() method chains together phrases based on last word of last
//    phrase, and first word of next phrase. Only strongest links are
//    used in the chain.
//
// The Dialogue() method uses a short snippet of dialogue to light
//    up a neural net of nodes.
//
// The PathWalk() method walks hill tops.
//
// HISTORY:
// January 1997 Linas Vepstas

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "multihash.h"
#include "wordhash.h"

#define LAG_THREE_WORD_PHRASE

#define LAG_USED_SIZE 50

class lagTextAnalysis {
   public:
      lagTextAnalysis (void);
      void Analyze (FILE *fh);
      void Dialogue (FILE *fh);
      void Dump (void);
      void HillClimb (void);
      void Chain (int);

      void PathWalk (void);
   protected:
      float TreeWalk (unsigned int lead_off_phrase,
                      unsigned int was_visited_arrray [],
                      short depth);
      void PrintPhrase (unsigned int);

   private:
      lagWordTable *word_table;
#ifdef LAG_THREE_WORD_PHRASE
      lagWordTripleTable *phrase_table;
#endif // LAG_THREE_WORD_PHRASE
#ifdef LAG_FOUR_WORD_PHRASE
      lagWordQuadTable *phrase_table;
#endif // LAG_FOUR_WORD_PHRASE
      lagConcordPairTable *xref_table;
      unsigned int was_used [LAG_USED_SIZE];
};

// =====================================================

lagTextAnalysis :: lagTextAnalysis (void) {
   word_table = new lagWordTable;
#ifdef LAG_THREE_WORD_PHRASE
   phrase_table = new lagWordTripleTable;
#endif // LAG_THREE_WORD_PHRASE
#ifdef LAG_FOUR_WORD_PHRASE
   phrase_table = new lagWordQuadTable;
#endif // LAG_FOUR_WORD_PHRASE
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

#ifdef LAG_THREE_WORD_PHRASE
            this_phrase_id = phrase_table -> GetID (blast_word_id, last_word_id, this_word_id);
#endif // LAG_THREE_WORD_PHRASE
#ifdef LAG_FOUR_WORD_PHRASE
            this_phrase_id = phrase_table -> GetID (clast_word_id, blast_word_id, last_word_id, this_word_id);
#endif // LAG_FOUR_WORD_PHRASE

            // light up the phrase neuron
            phrase_table -> AccumStrength (this_phrase_id, 1.0);

#ifdef LAG_THREE_WORD_PHRASE
            xref_table -> GetID (clast_phrase_id, this_phrase_id);
#endif // LAG_THREE_WORD_PHRASE
#ifdef LAG_FOUR_WORD_PHRASE
            xref_table -> GetID (dlast_phrase_id, this_phrase_id);
#endif // LAG_FOUR_WORD_PHRASE
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

#define LAG_NUM_CYCLES 5

void lagTextAnalysis :: Dialogue (FILE *fh) {

   xref_table -> ComputeWeights ();
   phrase_table -> ResetAllWeights ();

   // light up only this dialogue
   Analyze (fh);

   int num_phrases = phrase_table -> GetTableSize();

   for (int repeat=0; repeat < LAG_NUM_CYCLES; repeat ++) {
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
}

// =====================================================

void lagTextAnalysis :: HillClimb (void) {

   unsigned int hottest_phrase = 0;
   float heat = 0.0;

   // find the hottest phrase
   int num_phrases = phrase_table -> GetTableSize();
   int i = 0;
   for (i=1; i<= num_phrases; i++) {
      float activity = phrase_table -> GetWeight (i);
      if (activity > heat) {
         heat = activity;
         hottest_phrase = i;
      }
   }

   printf ("\nInfo: lagTextAnalysis :: HillClimb(): \n");
   printf ("Hotest phrase id = %d weight = %f \n", hottest_phrase, heat);
   PrintPhrase (hottest_phrase);
   printf ("\n");
   printf ("\nInfo: lagTextAnalysis :: HillClimb(): response text is \n");
   PrintPhrase (hottest_phrase);
   printf ("			id = %d weight = %f \n", hottest_phrase, heat);

   phrase_table -> AbsAllWeights();

   // setup to avoid infinite loops
   for (i=0; i<LAG_USED_SIZE; i++) {
      was_used [i] = 0;
   }

   // now follow the highest ridge
   while (0.05 < heat) {
      heat = 0.0;
      int found_one = 0;
      xref_table -> ResetToStart (hottest_phrase);
   
      unsigned int phrase = xref_table -> GetNextPhrase();
      while (phrase) {
         float activity = phrase_table -> GetWeight (phrase);

/*
printf ("next ph is %d %f \n", phrase, activity);
PrintPhrase (phrase);
printf ("\n");
*/

         if (activity > heat) {
            heat = activity;
            hottest_phrase = phrase;
            found_one = 1;
         }
         phrase = xref_table -> GetNextPhrase();
      }

      // avoid infinite loops
      for (i=0; i<LAG_USED_SIZE; i++) {
         if (hottest_phrase == was_used [i]) {
            printf ("\n\n");
            return;
         }
         if (!was_used[i]) {
            was_used[i] = hottest_phrase;
            // i = LAG_USED_SIZE + 50;  // break out of for loop
            break;
         }
         if (LAG_USED_SIZE-1 == i) return;
      }

      if (found_one) {
         PrintPhrase (hottest_phrase);
         printf ("			id = %d weight = %f \n", hottest_phrase, heat);
      } else {
         break;
      }
   }

   printf ("\n\n");
}

// =====================================================
#define LAG_TREE_DEPTH 7

float lagTextAnalysis :: TreeWalk 
   (unsigned int lead_off_phrase,
   unsigned int was_visited_array [],
   short depth) 
{

   float heat = 0.0;
   unsigned int hottest_phrase = 0;

   heat = phrase_table -> GetWeight (lead_off_phrase);
   if (LAG_TREE_DEPTH < depth) return heat;
   if (0.0005 > heat) return 0.0;

   heat = 0.00001;

   unsigned int best_path [LAG_USED_SIZE];
   int i = 0;
   for (i=0; i<LAG_USED_SIZE; i++) {
      best_path[i] = was_visited_array[i];
   }

   // now follow the highest ridge
   void * cursor = xref_table -> GetStart (lead_off_phrase);

   // walk precisely one level down the tree
   unsigned int phrase = xref_table -> GetPhrase(cursor);
   while (phrase) {

      // avoid infinite loops -- check to see if we've been here already
      short dont_do_it = 0;
      for (i=0; i<LAG_USED_SIZE; i++) {
         if (!was_visited_array [i]) break;
         if (phrase == was_visited_array [i]) {
            dont_do_it = 1;
            break;
         }
      }

      if (dont_do_it) {
         if (LAG_USED_SIZE-1 == i) break;
         cursor = xref_table -> GetNext (cursor);
         phrase = xref_table -> GetPhrase(cursor);
         continue;
      }

      // recursive inifinite loop trimmer.
      unsigned int stopper [LAG_USED_SIZE];
      for (i=0; i<LAG_USED_SIZE-1; i++) {
         stopper [i] = was_visited_array[i];
         if (!stopper[i]) {
            stopper[i] = phrase;
            stopper[i+1] = 0;
            break;
         }
      }

      // OK, now recurse down a level
      float path_sum = TreeWalk (phrase, stopper, depth+1);

      // locate hottest tree branch
      if (path_sum > heat) {
         heat = path_sum;
         hottest_phrase = phrase;
         for (i=0; i<LAG_USED_SIZE; i++) {
            best_path [i] = stopper[i];
            if (!stopper[i]) break;
         }
      }

      cursor = xref_table -> GetNext (cursor);
      phrase = xref_table -> GetPhrase(cursor);
   }

   // record the actual path taken
   for (i=0; i<LAG_USED_SIZE; i++) {
      was_visited_array[i] = best_path[i];
      if (!best_path[i]) break;

   }

   // add activation of this node to the hottest path.
   heat += phrase_table -> GetWeight (lead_off_phrase);
   return (heat);
} 

// =====================================================

void lagTextAnalysis :: PathWalk (void) 
{
   unsigned int hottest_phrase = 0;
   float heat = 0.0;

   // find the hottest phrase
   int num_phrases = phrase_table -> GetTableSize();
   int i = 0;
   for (i=1; i<= num_phrases; i++) {
      float activity = phrase_table -> GetWeight (i);
      if (activity > heat) {
         heat = activity;
         hottest_phrase = i;
      }
   }

   printf ("\nInfo: lagTextAnalysis :: PathWalk(): \n");
   printf ("Hotest phrase id = %d heat=%f \n", hottest_phrase, heat);
   PrintPhrase (hottest_phrase);
   printf ("\n");

   phrase_table -> AbsAllWeights();

   // setup to avoid infinite loops
   for (i=0; i<LAG_USED_SIZE; i++) {
      was_used [i] = 0;
   }

   was_used[0] = hottest_phrase;
   TreeWalk (hottest_phrase, was_used, 0);

   printf ("\nInfo: lagTextAnalysis :: PathWalk(): response text is: \n");
   for (i=0; i<LAG_USED_SIZE; i++) {
      if (was_used[i]) {
         PrintPhrase (was_used [i]);
         float heat = phrase_table -> GetWeight (was_used[i]);
         printf ("			id = %d weight = %f \n", was_used[i], heat);
      } else {
         break;
      }
   }
   printf ("\n\n");
}


// =====================================================

void lagTextAnalysis :: PrintPhrase (unsigned int phrase_id) {
   if (!phrase_id) return;

#ifdef LAG_THREE_WORD_PHRASE
   for (int j=0; j<3; j++) {
#endif // LAG_THREE_WORD_PHRASE
#ifdef LAG_FOUR_WORD_PHRASE
   for (int j=0; j<4; j++) {
#endif // LAG_FOUR_WORD_PHRASE
      unsigned int word_id = phrase_table -> GetElt (phrase_id, j);
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

   // follow least steepest descent
   texan -> HillClimb ();

   // follow highest ridge.
   texan -> PathWalk();

   texan -> Dump();

/*
   for (int i=0; i<LAG_TOP_TEN; i++) {
      printf ("\n\n =========== sentance %d =========== \n", i);
      texan -> Chain (i);
   }
*/

}

// ========================== END OF FILE ==================
