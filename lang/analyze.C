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

#include "pairhash.h"
#include "wordhash.h"

main () {

   FILE * fh = stdin;
   
   lagWordTable *wt = new lagWordTable;
   lagWordPairTable *pt = new lagWordPairTable;

   char buff[500];

   int this_word_id = 0;
   int last_word_id = 0;
   while (!feof(fh)) {
      fgets (buff, 500, fh);
      buff[499] = 0x0;

      // reject the mail headers
      if (!strncmp (buff, "From", 4)) continue;
      if (strstr (buff, "From:")) continue;
      if (strstr (buff, "Date:")) continue;
      if (strstr (buff, "Subject:")) continue;
      if (strstr (buff, "To:")) continue;
      if (strstr (buff, "Re:")) continue;
      if (strstr (buff, "Status:")) continue;
      if (strstr (buff, "MIME-Version:")) continue;
      if (strstr (buff, "Content-Type:")) continue;

      int i = 0;
      char * word = buff;
      while (buff[i]) {
      
         if (!isalpha (buff[i])) {
            buff[i] = 0x0;
            // printf (" doing %s \n", word);
            if (this_word_id) last_word_id = this_word_id;
            this_word_id = wt -> GetWordID (word); 
            pt -> GetWordPairID (last_word_id, this_word_id);
            word = &buff[i+1];
         } else {
            buff[i] = tolower (buff[i]);
         }
         i++;
      }
   }

   pt -> Dump();

   int i = 0;
   for (i=0; i<LAG_TOP_TEN; i++) {
      int top = pt -> GetTopTen (i);

      int fi = pt -> GetFirstOfPair (top);
      int se = pt -> GetSecondOfPair (top);
      int cnt = pt -> GetCount (top);
      char * wfi = wt -> GetWordFromID (fi);
      char * wse = wt -> GetWordFromID (se);
      printf (" its id=%d cnt=%d %s %s \n", top, cnt, wfi, wse);
   }

   // cnstruct string
   i = 0;
   int pair = pt -> GetTopTen (0);
   int fi = pt -> GetFirstOfPair (pair);
   char * wfi = wt -> GetWordFromID (fi);
   printf ("%s ", wfi);
   while (pair) {
      int se = pt -> GetSecondOfPair (pair);
      if (!se) break;
      char * wse = wt -> GetWordFromID (se);
      printf ("%s ", wse);
      i++;
      if (10 == i) {
         i = 0;
         printf ("\n");
      }
      pair = pt -> GetTopPairContainingWord (se);
   }

}

// ========================== END OF FILE ==================
