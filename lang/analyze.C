//
// analyze.C
//
// The lagWordtable class associates a unique ID number with a string.
// The GetWordID() method returns the ID number of the string.
//
// January 1997 Linas Vepstas

#include <ctype.h>
#include <stdio.h>

#include "wordhash.h"

main () {

   FILE * fh = stdin;
   
   lagWordTable *wt = new lagWordTable;

   char buff[500];

   while (!feof(fh)) {
      fgets (buff, 500, fh);
      buff[499] = 0x0;

      int i = 0;
      char * word = buff;
      while (buff[i]) {
      
         if (!isalpha (buff[i])) {
            buff[i] = 0x0;
            wt -> GetWordID (word); 
            word = &buff[i+1];
         } else {
            buff[i] = tolower (buff[i]);
         }
         i++;
      }
   }

   wt -> Dump();
}

// ========================== END OF FILE ==================
